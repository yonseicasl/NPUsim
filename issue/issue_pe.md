# PE 내부 논리 구현 진단

> **재판정 (2026-09-02): 아래 판정의 P0/P1 대부분이 이후 작업으로 해소됐다.**
>
> - **P0 activation → 해결.** ownership 계약이 확정·구현됨(plan/plan_sfu.md):
>   FUNCTIONAL 값의 유일한 owner는 Nebula `forward()`(activation 포함 적용, npu.cc의
>   `#ifdef FUNCTIONAL` 경로)이고, activation의 hardware cost는 SFU가 최종 output
>   commit 이벤트에서 모델링한다. `pe_t::activation()`(하드코딩 ReLU)은 **선언째
>   제거**됐고 재도입은 컴파일 에러 + `unittest/run_validation.sh` grep 가드.
>   SFU pure evaluator가 실제 Nebula 구현과 직접 링크 비교된다(validation/sfu).
> - **P1 mac_width → 해결.** lane→accumulator 계약 구현(PE1) + active-lane 과금(L9,
>   `utils/pe_lane.h`) + `mac_energy_<in>_<wt>` precision family(E3) —
>   [timing/pe.md](timing/pe.md) 재판정 참조.
> - **P1 절삭/시간축 → 해결.** cycle 통계 double 유지 + B12 경고(PE5), static
>   energy는 layer critical-path window(SE1–4).
> - **잔존:** sparse(functional 선행 조건, 전 계층 명시적 거부로 차단 — A7 정책),
>   PE-local shared/bypass(명시적 미지원 유지).
>
> "검증 상태"의 빌드 실패 언급도 stale — 현재 전체 빌드·validation 스위트 통과.

## 판정

현재 PE 구현은 기본 dense tile reduction의 메모리 안전성과 reset 안정성은 개선됐지만, **완전한 PE 모델로는 판단할 수 없다.** FUNCTIONAL 경로에는 layer activation이 연결되어 있지 않고, `num_macs`와 `mac_width`의 하드웨어 의미가 정의·구현되지 않았으며, sparse/shared/bypass 및 timing/energy 모델도 실제 하드웨어 동작과 일치하지 않는다.

- 심각도: **P0 / FUNCTIONAL 결과 신뢰성**, **P1 / 성능·에너지 결과 신뢰성**
- 분석 방법: 현재 작업 트리의 정적 코드·상태 전이·수식 분석 및 제공 validation 실행
- 영향 범위: ReLU convolution/FC의 FUNCTIONAL 빌드, `mac_width > 1` 구성(Simba 포함), sparse 실행 및 nonzero static-energy 구성

## 현재 해결된 항목

다음은 더 이상 미해결 이슈로 추적하지 않는다.

| 해결 항목 | 현재 근거 |
|---|---|
| output accumulator 범위 밖 접근 | 모든 MAC-side 배열이 `num_macs * mac_width` 용량으로 할당되고, tile/active-MAC 용량을 검증한다. [`pe.cc`](../components/pe.cc#L215-L221), [`pe.cc`](../components/pe.cc#L465-L478) |
| 새 output accumulator 초기화 누락 | 기존 partial sum이 없을 때 accumulator를 zero-fill한다. [`pe.cc`](../components/pe.cc#L1720-L1725) |
| PE reset 제어 상태 누락 | index, flush counter, 존재/request signal, offsets, buffer가 reset된다. [`pe.cc`](../components/pe.cc#L2444-L2501) |
| output load의 batch/group destination index | batch와 group 항의 덧셈 및 destination group divisor가 수정됐다. [`scheduler.cc`](../scheduler/scheduler.cc#L1596-L1603) |
| 잘못된 PE-local shared/bypass 동작 | 구현되지 않은 PE-local bypass와 shared buffer를 명시적으로 실패 처리해, 지원하는 것처럼 잘못 계산하지 않는다. [`pe.cc`](../components/pe.cc#L242-L246), [`pe.cc`](../components/pe.cc#L314-L325) |
| stationary/config 기본 검증 일부 | stationary, line size, bitwidth, MAC 수의 기본 유효성 검사가 추가됐고 output-stationary 기본 order가 수정됐다. [`pe.cc`](../components/pe.cc#L211-L214), [`pe.cc`](../components/pe.cc#L265-L325) |

`shared`/`bypass`는 "모델링 완료"가 아니라 **명시적 미지원** 상태다.

## 남은 문제

| 우선순위 | 문제 | 결과 |
|---|---|---|
| P0 | activation 함수가 정의만 되고 호출되지 않음 | `activation=relu`인 convolution/FC의 FUNCTIONAL output이 activation 없이 전달됨 |
| P1 | `num_macs`와 `mac_width`가 scalar register capacity 외의 의미를 갖지 않음 | vector-MAC lane reduction, accumulator 수, utilization 및 cycle의 하드웨어 의미가 불명확 |
| P1 | sparse functional path와 cost path가 분리됨 | compressed storage/decoder/irregular issue를 실행하지 않은 sparse 성능 추정 |
| P1 | cycle 계산에 `double`→`unsigned` 절삭 및 0-access 경계 위험이 남음 | fractional latency와 overlap cycle이 왜곡될 수 있음 |
| P1 | static energy가 elapsed time 기반이 아님 | leakage/static-energy 결과를 물리적 시간 기반 전력으로 해석할 수 없음 |
| P2 | PE-local shared/bypass는 미구현 | 이 기능을 요구하는 PE 구성은 실행할 수 없음 |

### 1. P0: activation이 reduction 완료 경계에도 적용되지 않는다

`activation()`은 accumulator에 ReLU를 적용하도록 구현돼 있지만 호출처가 없다. [`pe.cc`](../components/pe.cc#L2334-L2338) stationary 구현은 MAC과 output store를 수행할 뿐 activation을 호출하지 않는다. 예를 들어 input-stationary의 compute/store 경로는 [`pe.cc`](../components/pe.cc#L2586-L2633)에 있다.

따라서 기존의 "partial sum마다 ReLU" 오류는 제거됐지만, 현재는 반대로 ReLU가 전혀 적용되지 않는다. network config에는 실제로 `activation=relu`가 다수 존재한다. [`resnet50.cfg`](../models/networks/resnet50.cfg#L22-L52)

필요한 동작은 MAC의 곱셈·누산과 최종 activation을 분리하는 것이다. output identity별 C/R/S 및 array-level reduction 완료를 추적하고, 최종 output commit 직전에 layer policy(`relu`, `linear` 등)를 한 번 적용해야 한다. activation 비용을 PE에 포함할지 별도 operator로 모델링할지도 명시해야 한다.

### 2. P1: MAC-width 마이크로아키텍처가 정의되지 않았다

현재 `mac_register_capacity = num_macs * mac_width`이며 세 MAC-side storage를 그 크기로 할당한다. [`pe.cc`](../components/pe.cc#L215-L221) `mac_operation()`은 tensor index를 따라 output accumulator에 reduction을 수행한다. [`pe.cc`](../components/pe.cc#L2306-L2332) 기존 memory-safety 문제는 해결됐지만, `active_mac_width`는 여전히 사용되지 않는다. [`pe.h`](../components/pe.h#L186-L191)

따라서 Simba의 8 MAC × 8-wide가 64 scalar FMA인지, 8 accumulator와 8 multiplier-lane인지 불명확하다. 두 모델은 accumulator routing, register-port 수, lane reduction latency, utilization, energy가 다르므로 하나를 명시하고 mapping의 active component 수 및 cycle 계산과 일치시켜야 한다.

### 3. P1: sparse는 compressed execution이 아니다

FUNCTIONAL input transfer는 dense tensor를 MAC register로 복사한 뒤 compression type에 따라 비용을 분기한다. [`pe.cc`](../components/pe.cc#L633-L639), [`pe.cc`](../components/pe.cc#L743-L770) COO는 명시적으로 미지원 종료한다. [`pe.cc`](../components/pe.cc#L743-L746)

현재 sparse 결과는 **dense functional execution에 zero-gating 및 compressed traffic 추정을 덧붙인 모델**로 해석해야 한다. 압축 value/index/pointer storage, decoder, nonzero matching, irregular issue, accumulator routing을 구현하거나 sparse를 명시적으로 미지원 처리해야 한다.

### 4. P1: cycle/energy 모델의 단위와 시간축이 일관되지 않다

일부 pipeline stage는 unit cycle이 `double`인데도 `unsigned`에 저장한다. [`pe.cc`](../components/pe.cc#L1186-L1209) static energy는 elapsed cycle에 비례하지 않고 inactive PE에 transfer event마다 한 번 더한다. [`spatial_arch.cc`](../components/spatial_arch.cc#L1444-L1447)

cycle stage는 `double`로 유지하고 transaction count 0/1/2 이상의 공통 overlap helper를 사용해야 한다. static energy는 명시적 leakage-power 단위와 elapsed time/cycle을 곱해 계산해야 한다.

### 5. P2: PE-local shared/bypass는 지원하지 않는다

현재는 잘못된 데이터 경로 대신 초기화 단계에서 실패하므로 기능적 오계산은 방지한다. 다만 shared local buffer 또는 MAC forwarding/bypass를 연구 대상 하드웨어로 사용하려면 실제 구현이 필요하다.

- shared: 하나의 allocation, base offset, capacity/live-range allocator, bank/port contention
- bypass: 상위 memory/NoC에서 MAC register로 향하는 forwarding path와 LB access/cycle/energy 제거

## 검증 상태

- `./unittest/run_validation.sh`는 통과했다. 이 테스트는 config/parser와 mapping validation 중심이며 PE functional 수학 결과, activation, lane reduction을 검증하지 않는다. [`unittest/validation_test.cc`](../tests/validation_test.cc#L47-L102)
- 전체 library build는 현재 환경에 Nebula의 `convolutional.h`가 없어 완료하지 못했다. 따라서 end-to-end FUNCTIONAL 실행 검증은 보류 상태다.

## 완료 조건

- 음수·양수 partial product가 섞인 convolution이 `activation(sum(products))` reference와 일치한다.
- `linear` layer에는 ReLU가 적용되지 않는다.
- input/weight/output-stationary가 동일 layer에 대해 동일한 최종 tensor를 만든다.
- Simba 8×8의 lane/accumulator mapping, cycle, utilization, energy 정의가 문서와 테스트로 고정된다.
- dense와 지원 sparse format의 functional event와 cost event가 동일한 compressed data path를 사용한다. 미지원 format은 초기화 단계에서 거부한다.
- fractional unit cycle 및 0/1/2 transaction 경계에서 절삭·underflow 없이 예상 cycle이 나온다.
- nonzero static-energy 설정에서 elapsed simulation time에 비례한 leakage가 보고된다.

## 최종 결론

현재 PE는 **기본 dense tile reduction의 메모리 안전성과 reset 안정성은 개선됐지만**, activation 및 lane 마이크로아키텍처의 핵심 계약이 빠져 있다. 따라서 ReLU layer의 FUNCTIONAL 결과와 `mac_width > 1` PE의 성능/에너지 결과는 아직 신뢰할 수 없다.

```text
reduction 완료 추적 + layer-aware activation
→ num_macs/mac_width 의미와 cycle 모델 확정
→ sparse/shared/bypass 지원 범위 정리
→ cycle/static-energy 모델과 functional regression test
```
