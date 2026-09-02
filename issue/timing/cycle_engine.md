# Timing engine cycle 정합성 이슈

## 범위

이 문서는 2026-08-17 현재 작업 트리의 **timing simulation cycle 경로만** 재감사한 결과다.

- 포함: layer critical path, stage busy cycle, repetition scaling, overlap, fold/setup cycle, pipeline 수식, cycle 통계 및 golden validation 계약
- 제외: functional simulation, energy 값·단가·누적, traffic byte/transaction 정확도, sparse timing
- 기존 문서와 중복되는 일반적인 contention/back-pressure 미지원 자체는 반복하지 않고, 현재 코드에서 수치 불변식을 깨는 구체 구현만 기록한다.

최신 코드에는 mc→GLB fill을 `global_buffer_t::modeled_elapsed_cycles()`에 포함하는 보정과 output load-event 라벨 정정이 이미 반영되어 있으므로, 이 두 항목은 새 이슈에서 제외했다.

## 판정 (2026-08-18 갱신: CE1~CE7 전체 해결)

> **해결 후 상태**: 공식 metric("Compute-schedule latency")이 시뮬레이터에서 직접
> 출력되고 golden gate와 일치하며, critical path·busy·bottleneck은 최종 스케일된
> vector에서 재계산되어 T7/T8 invariant로 고정된다. Critical-path의 절대값은
> config 단가(GLB 1cyc/8b-line 등)를 반영하는 informational 축으로, 외부 검증은
> config 단가 캘리브레이션(후속) 이후 가능하다. 이하 원문은 해결 전 감사 기록.

Gemmini RTL과 Eyeriss silicon에 대해 현재 통과하는 latency 수치는 simulator가 출력하는 `Critical-path latency`가 아니다. 검증기는 별도 schedule proxy인 `Computation cycle + Fold fill cycle`을 사용하며, global timeline은 repetition 전의 stage busy를 일괄 배수로 확대한다. 그 결과 headline latency, bottleneck, MAC utilization과 golden validation 결과가 서로 다른 cycle 계약을 사용한다.

| ID | 우선순위 | 상태 | 요약 / 해결 (2026-08-18) |
| --- | --- | --- | --- |
| CE1 | P1 | **해결됨** | `stats_t::finalize_layer_timeline()` 신설 — scale_serial_repetitions 종단(양쪽 분기)에서 **최종 per-type vector로 stage busy·critical path·bottleneck·mac_available을 재계산**; leakage는 최종/uniform latency 비율로 선형 재조정. 검증: conv3 critical 644M→121.7M, bottleneck DRAM→GLB(문서의 예측대로 stale busy가 가리던 축 노출); gemm64 55,560 = 36,896+1,536+17,128 세그먼트 산술 정확 |
| CE2 | P1 | **해결됨** | compute 축 = `computation + fold_fill(+setup)`로 PE stage busy에 배치(`stage_axis_compute`), 결과에 **"Compute-schedule latency"로 직접 출력**. A/B 검증: `layer_setup_cycle` +10,000 → compute-schedule 정확히 +10,000, critical path는 세그먼트 max 계약대로 반응 |
| CE3 | P1 | **해결됨** | 공식 metric 계약 명문화: **"Compute-schedule latency" = 검증된 공식 metric**(RTL/silicon golden), "Critical-path latency" = memory-inclusive 분석 timeline(informational — config 단가 반영, 외부 미검증). traffic 하네스에 **T8**(critical ≥ compute-schedule) gate 추가. golden gate는 기존대로 compute-schedule 검사(이제 시뮬레이터가 그 값을 직접 출력) |
| CE4 | P1 | **해결됨** | 결합 규칙 확정: **단일 물리 link/port(shared) = 타입 합산, per-type partition(separate) = max** — DRAM device/link, NoP, GLB shared SRAM, 배포 fabric은 합산; separate GLB/LB는 max. 부수 발견·수정: **pass-through(temporal buffer 부재) array에 buffer access cycle이 무조건 과금**되던 3개 경로(`global_buffer.cc` 배포 write, `pe_array.cc` 배포 read·writeback read/write)를 `exist_temporal_buffer`로 게이트 |
| CE5 | P2 | **해결됨** | 표준 makespan `src+link+dst+(N−1)·max`로 교체(비대칭 1-cycle 과다 제거) + **count-aware 오버로드**(stage별 상이한 transaction 수) 신설, 6개 call site 전부 전환. unittest: (4,5,1,2)==23, bottleneck 3방향 대칭, width-ratio(1,8,8)==22 assert |
| CE6 | P3 | **해결됨** | `first_active_pe` 플래그로 MIN 3종(computation/mac/lb)을 전 chip에 걸쳐 1회만 초기화. 1-chip 회귀 무영향 확인; 2-chip live 테스트는 다중칩 config 부재로 미실행(코드 검토 검증) — 다중칩 활성화 시 테스트 추가 항목 |
| CE7 | P3 | **해결됨** | "Busy-cycle axes (access/link/overlap)" 5-stage 출력 신설 + traffic 하네스 **T7**(busy == max(printed axes); PE는 compute 축 포함) invariant — 8케이스 전부 통과 |

## CE1 — repetition 이후 stage busy와 critical path가 재계산되지 않음

### 코드 근거

`stats_t::update_stats()`는 한 번 실행한 live pass의 component state에서 다음을 먼저 계산한다.

```text
busy_cycle_{dram,multi_chip,global_buffer,pe_array,pe}
layer_latency
```

그 뒤 `scale_serial_repetitions()`는 다음처럼 처리한다.

1. `layer_latency`와 모든 `busy_cycle_*`를 global `m_repetitions`로 일괄 곱한다.
2. DRAM/multi-chip/fill의 실제 cycle vector는 `m_datatype_repetitions[type]`으로 서로 다르게 곱한다.
3. 바뀐 component cycle vector에서 stage busy와 critical path를 다시 계산하지 않는다.

관련 위치:

- `scheduler/stats.cc:449-479`
- `scheduler/stats.cc:779-869`

### 실행 근거

최신 Eyeriss AlexNet conv3(`result/eyeriss/alexnet/silicon/layer_4.txt`)에서:

| 항목 | cycle |
| --- | ---: |
| Headline critical path | 644,167,680 |
| Headline DRAM busy | 276,117,504 |
| Headline GLB busy | 184,025,088 |
| 최종 출력된 DRAM total-cycle 최대 | 5,309,952 |
| 최종 출력된 GLB access-cycle 최대 | 93,782,016 |

최종 세부 cycle만 보면 GLB 축이 DRAM 축보다 훨씬 크지만, stale busy 값 때문에 결과는 DRAM을 bottleneck으로 표시한다. 같은 레이어에서 DRAM busy는 최종 DRAM cycle 축의 약 52배다.

### 영향

- `Critical-path latency`가 최종 출력된 component cycle과 화해되지 않는다.
- bottleneck level이 바뀔 수 있다.
- `mac_available_cycle = physical_scalar_macs × layer_latency`이므로 MAC utilization도 함께 왜곡된다.
- GLB fill을 busy에 포함한 최신 보정은 `busy >= final access` 불변식을 복구하지만, uniform repetition 때문에 반대로 큰 과대값을 남긴다.

### 완료 조건

- 모든 repetition scaling과 one-time cycle을 적용한 후 stage busy를 최종 cycle vector에서 재계산한다.
- 같은 최종 stage busy로 critical path, bottleneck, MAC available cycle을 한 번만 계산한다.
- 각 stage에 대해 `busy`와 출력된 구성 axis 사이의 계약을 regression test로 고정한다.

## CE2 — Fold fill과 layer setup cycle이 critical path에서 누락됨

### 코드 근거

`fold_fill_cycle`은 PE-array WEIGHT distribution에서 누적되지만 `pe_array_t::modeled_elapsed_cycles()`는 다음만 읽는다.

- `write_back_cycle`, `overlapped_transfer_cycle`
- `access_cycle[type]`
- `transfer_cycle[type]`
- `cycle_temporal_pe[type]`

`fold_fill_cycle`과 `u_layer_setup_cycle`은 읽지 않는다. 더구나 `stats_t::update_stats()`가 PE-array busy와 layer latency를 계산한 뒤에야 `fold_fill_cycle_pe_array`를 수집하고, `scale_serial_repetitions()`에서 setup을 추가한다.

관련 위치:

- `components/pe_array.cc:136-143`
- `components/pe_array.cc:181-190`
- `scheduler/stats.cc:449-479`
- `scheduler/stats.cc:620-622`
- `scheduler/stats.cc:772-829`

### 실행 근거

Gemmini `gemm_64x64x64` 결과:

```text
Computation cycle = 1,024
Fold fill cycle   = 2,494  (224 fold + 2,270 one-time setup)
```

Fold/setup이 computation보다 2.44배 큰데도 PE-array busy와 `Critical-path latency` 계산에는 이 값이 들어가지 않는다. 현재 RTL 비교가 이 값을 사용하는 것은 validation 스크립트가 결과 파일에서 별도로 더하기 때문이다.

### 영향

- simulator의 headline latency가 RTL calibration에 사용한 cycle을 누락한다.
- fold/setup knob을 바꿔도 다른 stage가 지배하는 한 critical path와 bottleneck이 반응하지 않는다.
- `MAC utilization (time)`의 분모가 fold/setup stall을 포함하지 않는다.

### 완료 조건

- fold-fill을 PE-array/compute schedule의 명시적 serial 또는 overlap axis로 배치한다.
- one-time setup은 repetition과 분리해 layer timeline에 정확히 한 번 포함한다.
- `Critical-path latency`가 fold/setup knob 변화량만큼 예측 가능하게 반응하는 테스트를 추가한다.

## CE3 — Golden validation과 headline latency가 서로 다른 지표를 사용함

### 코드 근거

`validation/check_timing.py`는:

- Gemmini: `Computation cycle + Fold fill cycle`
- Eyeriss: `Computation cycle + Fold fill cycle`

만 golden과 비교한다. `Critical-path latency`는 파싱하거나 assert하지 않는다. Phase-1은 다시 `max(computation, PE-array input IC, weight IC)`라는 세 번째 proxy를 사용한다.

관련 위치:

- `validation/check_timing.py:28-44`
- `validation/check_timing.py:47-79`
- `validation/phase1/compare.py:38-53, 81-96`

### 실행 근거

현재 `validation/check_timing.py --check-baseline`은 통과한다.

```text
Gemmini RTL: MAPE=4.40%, max=7.86%
Eyeriss silicon latency: MAPE=4.26%, max=6.39%
```

그러나 Eyeriss의 headline critical path를 동일 파일의 검증 schedule과 비교하면:

| 레이어 | Critical / (Computation + Fold) |
| --- | ---: |
| conv1 | 277.5× |
| conv2 | 145.7× |
| conv3 | 151.4× |
| conv4 | 168.4× |
| conv5 | 190.0× |
| fc1 | 1,188.8× |
| fc2 | 1,610.2× |
| fc3 | 2,734.6× |

따라서 현재의 4~8% 정확도는 compute-schedule proxy의 정확도이며 global critical-path 정확도로 해석할 수 없다.

### 영향

- regression이 통과해도 headline latency는 임의로 변하거나 크게 틀릴 수 있다.
- 사용자에게 어떤 latency를 공식 결과로 사용해야 하는지 명확한 계약이 없다.
- Phase-1/2/3이 서로 다른 cycle proxy를 사용해 하나의 timing 정확도 수치로 합치기 어렵다.

### 완료 조건

- 공식 latency metric을 하나로 정의한다.
- global timeline을 공식 metric으로 선택하면 golden gate가 `Critical-path latency`를 직접 검증해야 한다.
- compute-only/array-schedule/memory-inclusive latency가 모두 필요하면 서로 다른 이름과 acceptance threshold로 분리한다.

## CE4 — Shared resource의 datatype cycle을 `max`로 완전 overlap함

### 코드 근거

다음 `modeled_elapsed_cycles()`는 INPUT, WEIGHT, OUTPUT의 cycle을 모두 `max`로 결합한다.

- `global_buffer_t::modeled_elapsed_cycles()`
- `multi_chip_t::modeled_elapsed_cycles()`
- `dram_t::modeled_elapsed_cycles()`
- `pe_array_t::modeled_elapsed_cycles()`
- `pe_t::modeled_elapsed_cycles()`

이 계산은 component의 `memory_type=SHARED/SEPARATE`, 포트 수, link 공유 여부를 사용하지 않는다. 특히 DRAM link와 NoP/NoC link는 datatype별로 별도 물리 링크가 정의되어 있지 않으므로 세 stream을 무조건 동시에 처리할 수 없다.

관련 위치:

- `components/global_buffer.cc:1586-1597`
- `components/multi_chip.cc:1838-1845`
- `components/dram.cc:1457-1464`
- `components/pe_array.cc:136-143`
- `components/pe.cc:2683-2694`

### 실행 근거

Gemmini square GEMM에서 최종 DRAM transfer cycle은 INPUT/WEIGHT/OUTPUT이 동일하다. 예를 들어 `gemm_64x64x64`은 각 1,536 cycles로:

```text
max(type transfer) = 1,536
sum(type transfer) = 4,608
```

현재 component elapsed 계약은 공유 링크에서도 max를 택한다. CE1의 stale uniform busy가 이 차이를 가리지만, CE1을 고쳐 최종 vector로 busy를 재계산하면 shared-resource 과소 계상이 바로 드러난다.

### 영향

- shared SRAM의 input/weight/output read-write 충돌이 숨겨진다.
- DRAM/NoP/NoC의 datatype arbitration 없이 최대 3-way full overlap을 가정한다.
- `SEPARATE`와 `SHARED` memory type의 cycle 차이가 용량 외에는 충분히 반영되지 않는다.

### 완료 조건

- separate banks/ports와 shared port/link에 서로 다른 결합 규칙을 적용한다.
- 최소 analytical 모델에서는 shared single-port/link의 service cycle을 합산하고, multi-port는 port 수만큼 스케줄한다.
- 세 datatype stream을 하나씩 활성화한 결과의 합과 동시 활성 결과를 비교하는 contention micro-test를 추가한다.

## CE5 — `pipelined_transfer_cycles()`의 수식과 API가 일반적이지 않음

### 동일 transaction 수에서도 발생하는 수식 오류

동일한 `N`개 token이 source/link/destination 3-stage pipeline을 통과하는 ideal makespan은 다음이다.

```text
source + link + destination + (N - 1) × max(source, link, destination)
```

현재 구현은 양 끝 stage에 별도 `max`를 적용한다. 비대칭 stage 비용에서 1 cycle 과다다.

| N, source, link, destination | 현재 | ideal |
| --- | ---: | ---: |
| 2, 5, 1, 2 | 14 | 13 |
| 4, 5, 1, 2 | 24 | 23 |
| 4, 1, 5, 2 | 23 | 23 |
| 4, 2, 1, 5 | 24 | 23 |

`unittest/validation_test.cc`는 `(4,5,1,2)==24`를 정답으로 assert하고 있어 현재 오류를 regression contract로 고정한다.

관련 위치:

- `utils/interconnect_timing.cc:70-80`
- `unittest/validation_test.cc:232-237`

### Endpoint transaction 수가 다를 때의 API 한계

`datatype_transfer_timing()`은 source access, destination access, link transaction 수가 서로 다를 수 있지만 `pipeline_transactions = max(세 수)` 하나만 만든다. 그 뒤 helper는 이 최대 횟수만큼 세 stage가 모두 일한다고 가정한다.

예를 들어 256-bit source line, 32-bit link, 32-bit destination처럼 width가 다르면 한 source access가 여러 link/destination token을 만든다. 올바른 stage workload는 각각:

```text
source_accesses × source_cycle
link_transactions × link_cycle
destination_accesses × destination_cycle
```

인데 현재 API는 세 개의 count를 보존하지 않아 큰 cycle 차이를 표현할 수 없다.

### 완료 조건

- 동일 token 3-stage pipeline 수식을 표준 makespan으로 교체한다.
- 서로 다른 endpoint/link count를 받는 workload 기반 API를 도입하거나 width-conversion stage를 명시한다.
- source/link/destination 중 각각 하나가 bottleneck인 비대칭 테스트와 width-ratio 테스트를 추가한다.

## CE6 — PE MIN cycle 통계가 chip마다 재초기화됨

### 코드 근거

`stats_t::update_stats()`는 active chip `i` 바깥 루프, active PE `j` 안쪽 루프 구조다. 그러나 MIN 초기화 조건은 `j == 0`만 사용한다.

```cpp
min_computation_cycle = (j == 0) ? current : min(min_computation_cycle, current);
min_access_cycle_mac[k] = (j == 0) ? current : min(...);
min_access_cycle_lb[k] = (j == 0) ? current : min(...);
```

두 번째 chip의 첫 PE에서 이전 chip의 global minimum이 사라진다.

관련 위치:

- `scheduler/stats.cc:493-536`

예를 들어 chip0 PE cycle이 `[1, 100]`, chip1이 `[50, 60]`이면 전체 minimum은 1이지만 현재 결과는 50이다.

### 영향

- multi-chip 결과의 `MIN` 및 `MIN access cycle`이 마지막 active chip의 minimum만 나타낸다.
- layer latency는 max 기반이라 직접 변하지 않지만 load imbalance와 tail-chip 진단이 틀린다.

### 완료 조건

- 최초 active PE인지 나타내는 전체 counter 또는 `i==0 && j==0` 조건으로 한 번만 초기화한다.
- 서로 다른 chip별 cycle을 가진 2-chip 통계 테스트를 추가한다.

## CE7 — Stage busy를 구성하는 cycle axis가 결과에 모두 노출되지 않음

`modeled_elapsed_cycles()`는 PE-array와 multi-chip의 access, transfer, overlap/writeback axis를 모두 사용한다. 그러나 결과 파일은:

- PE-array: interconnection과 temporal overlap은 출력하지만 `access_cycle_pe_array`는 출력하지 않음
- multi-chip: interconnection은 출력하지만 `access_cycle_multi_chip`과 overlap total을 충분히 출력하지 않음

따라서 `Busy cycles` 값이 어떤 세부 축에서 나왔는지 결과 파일만으로 재구성할 수 없다. CE1처럼 busy와 최종 axis가 어긋나도 자동·수동 감사가 어렵다.

관련 위치:

- `scheduler/stats.cc:1047-1065`
- `scheduler/stats.cc:1300-1360` 부근 PE-array 출력
- `scheduler/stats.cc:1498-1555` 부근 multi-chip 출력

### 완료 조건

- 각 stage에 대해 busy 계산에 사용한 모든 axis 또는 최소한 `dominant_axis`와 그 값을 출력한다.
- 결과 생성 시 `busy == recompute(printed axes)` invariant를 검사한다.

## 권장 해결 순서

1. **CE2**: fold/setup을 timeline 입력으로 연결해 누락 cycle을 제거한다.
2. **CE1**: repetition 이후 최종 vector로 stage busy와 critical path를 재계산한다.
3. **CE4**: shared/separate 및 link serialization 계약을 확정한다.
4. **CE5**: pipeline helper를 workload/count-aware 수식으로 교체한다.
5. **CE3**: 확정된 공식 latency를 golden gate가 직접 검증하도록 전환한다.
6. **CE6/CE7**: multi-chip MIN과 cycle 관측성을 보강한다.

CE1~CE5를 먼저 해결하지 않고 현재 `Critical-path latency`를 baseline으로 고정하면, 잘못된 uniform scaling과 hidden overlap을 golden 값으로 굳힐 위험이 있다.

## 재현에 사용한 명령

```bash
./npusim.sh build npusim
./unittest/run_timing_validation.sh
python3 validation/check_timing.py --check-baseline
```

검증 결과 자체는 통과했지만, 위에서 설명한 대로 그 gate는 `Critical-path latency`를 검사하지 않는다.
