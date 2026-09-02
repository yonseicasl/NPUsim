# NPUsim timing simulation issue 재감사

- 감사일: 2026-08-18
- 대상: `issue/timing/*.md` 13개 문서와 현재 timing 코드
- 포함: dense timing, cycle 집계, overlap, buffer/NoC/NoP/DRAM timing, runtime datatype, static/timing energy, 검증 계약
- 제외: functional-only 실행 정확성

## 최종 판정

Dense single-chip의 레이어 단위 timing 경로는 상당 부분 구현되어 있다. Gemmini RTL과 Eyeriss silicon 기준의 compute-schedule 검증도 통과한다. 그러나 `issue/timing/cycle_engine.md`의 “CE1~CE7 전체 해결”을 포함해 timing issue 전체가 완료된 상태는 아니다.

핵심 미완료 영역은 다음과 같다.

1. network timing 총계의 신규 axis 누락
2. 공식 compute-schedule 출력과 golden checker의 직접 연결 부재
3. 병렬 PE/chip과 shared resource의 `max`/`sum` 집계 순서
4. pass-through PE-array와 GLB bypass 의미론
5. width-conversion pipeline 및 single-buffer LB 직렬화
6. multi-chip multicast, contention, back-pressure
7. stage-granular 모델을 넘어선 per-tile global timeline
8. 상세 DRAM, format-IP, PE reduction, systolic pipeline 모델과 calibration

문서의 이슈 ID에는 중복 별칭이 있다. 예를 들어 `SP3`, `LB5`, `GB7`, `MC5`는 같은 `pipelined_transfer_cycles()` 문제를 참조한다. 따라서 단순 ID 개수 대신 문서와 실행 경로별로 분류한다.

## 해결된 이슈

| 문서 | 해결된 항목 | 현재 판정 근거 |
|---|---|---|
| [adder_tree](../issue/timing/adder_tree.md) | AT1, AT2, AT3, AT4, AT6, AT7 | reduction 비용이 OUTPUT write-back override에서 fan-in과 data volume에 따라 과금되고 overlap 통계에 연결됨 |
| [cycle_engine](../issue/timing/cycle_engine.md) | CE1의 레이어 단위 동작 | repetition 이후 최종 vector로 레이어 stage busy와 critical path를 재계산함 |
| [dram](../issue/timing/dram.md) | DR1, memory-controller DR4, Phase-3 DR4, DR5, DR6 | dense guard, queue UB, datatype repetition, grouped-conv 의미론, final output flush가 구현됨 |
| [global_buffer](../issue/timing/global_buffer.md) | GB1~GB5의 dense/write-back 수정, legacy weight 이중 계상, fill per-type scaling, GB8, GB9 | descriptor 기반 write-back, double-buffer gating, utilization, fill 분리, 출력 및 busy 축 연결이 존재함 |
| [local_buffer](../issue/timing/local_buffer.md) | LB1, LB2, LB4, LB6, LB5의 1-transaction 항 | utilization, leakage, zero-active guard와 단일 transaction 수식이 구현됨 |
| [multi_chip](../issue/timing/multi_chip.md) | MC1~MC6의 기본 dense 경로 | OUTPUT write-back, leakage, 요청 guard, NoP 단가 초기화, 1-transaction 오류가 수정됨 |
| [pe](../issue/timing/pe.md) | PE3, PE4, PE6, PE7, PE8 | peak utilization hardening, leakage window, zero-active guard, 방어 코드, MAC available window가 연결됨 |
| [pe_array](../issue/timing/pe_array.md) | PA1~PA8 | shared-stream multicast, leakage, elapsed window, overlap report, topology, occupancy, fold fill, edge accumulation 구현 |
| [runtime_datatype](../issue/timing/runtime_datatype.md) | dense packing, MXFP metadata, hierarchy descriptor와 report | dense runtime datatype의 storage/access/link transaction 경로가 공통 descriptor를 사용함 |
| [spatial_architecture](../issue/timing/spatial_architecture.md) | SP1, SP2, SP3의 기존 단일 transaction 문제 | analytical hop fill, topology validation, 1-transaction 수식이 구현됨 |
| [static_energy](../issue/timing/static_energy.md) | SE1~SE4 | 레이어 critical-path leakage window, always-on domain, multi-chip leakage, layer/network report 연결 |
| [systolic_array](../issue/timing/systolic_array.md) | SY1~SY5의 analytical 범위, V2, V3 | active array diameter fill, topology guard, fold-fill calibration, edge accumulator가 구현됨 |

## 부분 해결된 이슈

| 항목 | 재판정 | 남은 내용 |
|---|---|---|
| CE2 | 부분 해결 | 레이어 compute schedule은 정상이나 network 총계의 `stage_axis_compute`가 0 |
| CE3 | 부분 해결 | 공식 metric을 출력하지만 golden checker는 출력값을 직접 검사하지 않음 |
| CE4 | 부분 해결 | datatype별 sum/max는 추가됐지만 병렬 entity와 datatype 결합 순서가 부정확함 |
| CE5 / SP3 / LB5 / GB7 / MC5 | 부분 해결 | 동일 count와 단일 transaction은 수정됐으나 width-conversion packet dependency가 없음 |
| CE6 | 구현 해결, 검증 미완료 | `first_active_pe`는 수정됐지만 2-chip live regression이 없음 |
| CE7 | 부분 해결 | 레이어 T7은 통과하나 network axis가 전부 0 |
| DR2 | 부분 해결 | 순차 dense stream의 row activation만 모델링. bank conflict, random/strided access, tRC/tRAS/tRP 없음 |
| DR3 / PE5 | 부분 해결 | zero/frequency guard와 fractional 경고는 있으나 절삭 및 non-power-of-two 허용은 유지 |
| Global timeline A1 | 부분 해결 | 5-stage analytical max/sum 모델만 있음. tile ready/available cycle과 stall propagation 없음 |
| Multi-chip route | 부분 해결 | Manhattan hop은 있으나 shared-link arbitration과 back-pressure 없음 |
| PE1 / PE2 | 부분 해결 | `mac_width`가 reduction cycle에 연결됐지만 reduction energy와 active-lane 구조 모델은 미완성 |
| Runtime format-IP | 부분 해결 | configurable cycle/energy hook은 있으나 throughput 및 다른 stage와의 overlap 없음 |
| Systolic SY2 | analytical 해결 | diameter fill은 있으나 operand skew, pipeline registers, partial-sum propagation, tile-boundary bubble 없음 |
| PA9 | 외부 정보 대기 | Eyeriss mapper parameter 부재로 GLB point-match는 완료할 수 없고 경계 검증만 종결됨 |

## 명확히 미해결된 이슈

### 1. Network timing 총계 axis 누락

`scheduler/stats.cc::update_network_stats()`는 `layer_latency`와 `busy_cycle_*`는 합산하지만 다음 필드를 합산하지 않는다.

- `stage_axis_compute`
- `stage_axis_access[5]`
- `stage_axis_link[5]`
- `stage_axis_overlap[5]`

최신 `result/eyeriss/alexnet/silicon/network.txt` 재현값:

```text
Compute-schedule latency :        0.0 cycles
Critical-path latency    : 2408877220.0 cycles
Computation cycle        :   21621248.0 cycles
Fold fill cycle          :    2378337.3 cycles
```

올바른 network compute schedule 합은 `23,999,585.3` cycle이다. 모든 printed busy-cycle axis도 0이어서 network-level CE2/CE7은 실패한다.

### 2. Golden checker가 공식 metric을 직접 검증하지 않음

`validation/check_timing.py`는 `Compute-schedule latency`를 읽지 않고 다음 두 값을 따로 파싱해 합산한다.

- `Computation cycle`
- `Fold fill cycle`

따라서 공식 출력이 깨져도 golden validation이 통과할 수 있다. 위 network의 `0.0` 문제가 기존 gate에서 발견되지 않는 이유다.

필요한 검증 계약:

1. golden 비교량은 `Compute-schedule latency`를 직접 파싱한다.
2. 별도 invariant로 `Compute-schedule latency == Computation cycle + Fold fill cycle`을 검사한다.
3. layer뿐 아니라 `network.txt`에서도 축 합산을 검사한다.

### 3. 병렬 entity와 datatype의 집계 순서

현재 stats는 datatype별로 먼저 PE/chip의 최대치를 구한 뒤 datatype을 합산한다.

```text
현재: sum_type(max_entity(cycle[entity][type]))
필요: max_entity(sum_type(cycle[entity][type]))
```

예를 들어 chip 0이 input 100 cycle, chip 1이 weight 100 cycle이면 실제 병렬 elapsed는 100 cycle이지만 현재 계산은 200 cycle이 된다. shared GLB, PE-array fabric, PE local link 등에서 비대칭 mapping의 critical path를 과대평가할 수 있다.

### 4. Pass-through PE-array에 존재하지 않는 destination stage 과금

`global_buffer_t::account_descriptor_dense_transfer()`는 `pe_array->exist_temporal_buffer == false`일 때 PE-array write access/energy를 생략한다. 그러나 `cycle_pe_array_global_buffer`의 `pipelined_transfer_cycles()`에는 여전히 다음이 들어간다.

- `timing.destination_accesses`
- `pe_array->u_write_cycle[type]`

즉 존재하지 않는 temporal-buffer write stage가 overlap cycle에는 남는다. 현재 배포된 다수 config는 PE-array temporal buffer가 기본 off이므로 live 결과에 영향을 줄 수 있다.

### 5. GLB bypass 의미론 불일치

GLB `bypass`는 capacity/utilization 검사에만 사용되고 descriptor access/transfer accounting에는 적용되지 않는다. `eyeriss.cfg` 등의 `bypass=0:1:0`은 weight GLB 접근을 실제로 제거하지 않는다.

필요한 결정:

- bypass를 direct stream으로 정의해 GLB access를 생략하거나
- 현재 물리 모델에 맞게 config에서 bypass를 제거한다.

### 6. Single-buffer LB serialization 부재

`pe_t::modeled_elapsed_cycles()`는 다음을 모두 `max`로 결합한다.

- computation
- LB access
- LB↔MAC transfer
- write-back

LB는 single-buffered인데 tile load와 compute가 full overlap된 것으로 처리된다. LB7은 여전히 미해결이며 single/double-buffer 설정을 cycle상 구분하는 모델이 필요하다.

### 7. Multi-chip multicast와 contention

현재 NoP distribution은 active chip마다 source read와 link stream을 반복 과금한다. 다음 기능은 없다.

- 공유 input/weight multicast
- 동일 physical link를 공유하는 route 간 contention
- arbitration
- destination buffer pressure
- producer/consumer stall propagation
- 신뢰 가능한 2-chip 이상 회귀

### 8. Width mismatch pipeline dependency

Count-aware `pipelined_transfer_cycles()`는 각 stage의 남은 작업이 이상적으로 overlap된다고 가정한다. `8-bit source line -> 256-bit link -> 8-bit destination line`처럼 여러 source access를 모아 하나의 link transaction을 구성하는 경우의 assemble/disassemble dependency는 표현하지 못한다.

현재 unit test에는 `(1,8,8)` 방향만 있고 실제 config에서 흔한 `(32,1,32)`와 tail packet 사례가 없다.

## 의도적으로 보류된 항목

| 항목 | 보류 사유 |
|---|---|
| AT5 adder cycle/energy calibration | 모델 hook은 있으나 측정·문헌 기반 단가가 없음 |
| PA9 Eyeriss GLB point-match | 칩 mapper parameter가 공개되지 않음 |
| Sparse CSC/COO timing | functional sparse execution이 선행 조건이며 현재 top-level에서 명시적으로 거부됨 |
| 상세 DRAM timing | DRAMSim3 또는 별도 bank/timing model과 calibration 필요 |
| Format-IP 절대 단가·throughput | 하드웨어별 calibration 필요 |
| Router queue/VC/adaptive routing | per-tile/global event timeline 위의 별도 topology model 필요 |
| Functional-only GB6/LB3 | 이번 timing-only 감사 범위 밖 |

## 문서 상태 불일치

현재 issue 문서에는 해결 전 본문과 해결 후 표기가 함께 남아 있어 구현 계획의 기준으로 바로 사용하기 어렵다.

- `cycle_engine.md`: CE1~CE7 전체 해결 표기가 부정확
- `systolic_array.md`: 상단 판정과 표의 해결 상태가 충돌
- `global_cycle_overlap.md`: “timeline 없음” 설명과 stage timeline 구현 판정이 혼재
- `pe_array.md`: PA1 본문은 해결됐으나 완료 체크박스는 미완료
- `local_buffer.md`: LB2/LB4 본문은 해결됐으나 완료 체크박스는 미완료
- `pe.md`: PE1/PE3 본문과 완료 조건 불일치
- `adder_tree.md`: 해결된 PA2~PA4를 미해결로 가리키는 참고문이 남음
- `dram.md`: memory-controller 수정 후에도 완료 조건이 미완료로 남음

## 검증 실행 결과

### 실행 명령

```bash
bash unittest/run_validation.sh
bash unittest/run_timing_validation.sh
python3 validation/traffic/check.py
./npusim.sh run eyeriss_energy alexnet silicon
python3 validation/energy/check.py
bash unittest/run_full_sanitizers.sh
```

### 결과

| 검증 | 결과 |
|---|---|
| Unit/config validation | 통과, 203개 config 검증 |
| Full ASan/UBSan smoke | 통과, Gemmini GEMM64와 Eyeriss AlexNet |
| Gemmini RTL timing | MAPE 4.40%, max 7.86%, 기준 8% 이내 |
| Eyeriss silicon latency | MAPE 4.26%, max 6.39%, 기준 이내 |
| Traffic T1~T8 | 8개 대상 전부 통과 |
| Energy E1~E5b | 8개 레이어 전부 통과 |
| Network-level official metric | 실패: compute schedule 0.0 |
| Network-level busy axes | 실패: 5-stage axis 전부 0.0 |
| 2-chip CE6 live regression | 미완료: 기존 TPUv3/MobileNetV3 config는 있으나 runtime fixture가 없음 |
| 전체 AlexNet timing | 부분 결과: pooling/softmax 등 5개 non-convolution/connected layer 제외 |

기존 gate 통과는 dense single-chip의 선택된 레이어 경로가 안정적이라는 뜻이다. Network aggregate, multi-chip, bypass, single-buffer, packetization의 정확성을 증명하지는 않는다.

## 해결 우선순위

### P0 — 결과 계약과 회귀 gate 복구

1. `update_network_stats()`에 `stage_axis_compute/access/link/overlap` 합산 추가
2. `network.txt`에서 compute schedule과 busy axis invariant 추가
3. `check_timing.py`가 `Compute-schedule latency`를 직접 golden과 비교하도록 변경
4. `Compute-schedule == computation + fold_fill` identity gate 추가

이 단계가 끝나야 이후 cycle 변경이 공식 metric과 network 결과에 일관되게 반영된다.

### P1 — 현재 critical path를 직접 왜곡하는 문제

1. CE4의 entity-local 합산 후 global max로 집계 순서 교체
2. pass-through PE-array에서 nonexistent write stage 제거
3. GLB bypass 계약 확정 및 descriptor timing에 연결
4. single-buffer LB의 load/compute serialization 구현
5. width-conversion packet dependency를 반영하는 pipeline API와 micro-test 추가

### P2 — Multi-chip timing 완성

1. 다운로드 없는 최소 2-chip synthetic fixture 추가
2. CE6 min 통계와 entity aggregation 회귀 추가
3. NoP multicast/source sharing 계약 구현
4. route별 shared-link contention과 arbitration 구현
5. buffer pressure와 stall propagation을 timeline에 연결

### P3 — Per-tile global timing 정밀화

1. 계층별 `ready_cycle`과 `available_cycle` 도입
2. tile 단위 producer/consumer interval 기록
3. double-buffer depth에 따른 overlap 허용량 적용
4. fill/drain과 back-pressure를 layer critical path에 반영
5. leakage와 MAC availability를 새 wall-clock에 재연결

### P4 — 컴포넌트 상세 모델과 calibration

1. DRAM bank/row state와 tRC/tRAS/tRP 또는 DRAMSim3 필수 경로
2. PE lane reduction energy와 active-lane 구조
3. systolic skew, partial-sum propagation, tile-boundary bubble
4. format-IP throughput/overlap
5. adder-tree와 component unit cost calibration

### Deferred — 선행 조건 이후

1. functional sparse execution 구현
2. CSC/COO value/index/pointer 및 decoder timing 구현
3. 외부 mapper·silicon 세부 데이터 확보 후 PA9 point calibration

## 완료 판정 기준

Timing simulation 전체 완료를 선언하려면 최소한 다음 조건이 필요하다.

1. layer와 network에서 공식 metric 및 stage axis가 모두 자기 일관적이다.
2. single-chip과 multi-chip 모두 shared resource aggregation 수기 사례를 통과한다.
3. buffer 존재/bypass/single-double 설정 변경이 예측 가능한 cycle delta를 만든다.
4. equal-width와 width-conversion pipeline 사례가 packet-level 수기 계산과 일치한다.
5. per-tile timeline 또는 명시적으로 제한된 analytical 계약이 back-pressure와 overlap 의미를 일관되게 정의한다.
6. 모든 지원 workload가 부분 결과 여부를 명확히 표시하고, 공식 gate가 simulator의 직접 출력값을 검증한다.

