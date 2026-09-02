# NPUsim Energy·Traffic 정합성 검증

## 검증 개요

- 검증일: 2026-08-17
- 기준 커밋: `a9d76b7` + 현재 작업 트리의 energy/traffic 수정분
- 범위: timing simulation의 dense 경로
- 대상: Eyeriss AlexNet(CONV 5개, FC 3개), Gemmini GEMM 3개
- 제외: functional simulation 정확도, sparse/compressed timing 경로, 절대 전력 캘리브레이션

## 결론

Traffic의 핵심 내부 항등식과 GLB fill energy의 datatype별 반복 스케일링은 정상 동작한다. 하지만 mc→GLB fill cycle이 layer busy-cycle 계산에서 빠지는 코드 정합성 결함이 남아 있어, 현 상태를 energy/traffic 수정 완료로 판정할 수는 없다. 또한 Eyeriss 외부 traffic 기준 오차와 검증 자동화의 누락도 남아 있다.

## 수행한 검증

```bash
./npusim.sh build npusim
./npusim.sh run eyeriss alexnet silicon
./npusim.sh run gemmini gemm_64x64x64 ws
./npusim.sh run gemmini gemm_256x256x256 ws
./npusim.sh run gemmini gemm_512x512x512 ws
python3 validation/traffic/check.py

./npusim.sh run eyeriss_energy alexnet silicon
python3 validation/energy/check.py

./unittest/run_validation.sh
./unittest/run_full_sanitizers.sh
python3 validation/check_timing.py
python3 validation/check_timing.py --check-baseline
python3 validation/check_timing.py --check-traffic
python3 validation/phase3/compare.py
git diff --check
```

## 통과 결과

### Traffic 내부 정합성

최신 바이너리로 결과를 다시 생성한 후 `validation/traffic/check.py`의 T1~T6을 실행했다.

| 검사 | 결과 |
| --- | --- |
| T1 output DRAM bytes == output volume | 8/8 케이스 `1.000` |
| T2 weight DRAM bytes >= weight volume | 8/8 케이스 `1.00` |
| T3 input DRAM bytes >= input volume | 통과; CONV `2.60~4.90`, GEMM `1.00` |
| T4 GLB bytes >= DRAM bytes | 8/8 통과 |
| T5 GLB weight bytes == event × array tile | 8/8 케이스 `1.000` |
| T6 computation cycles == MACs / active PEs | 8/8 케이스 `1.000` |

DR6 final output flush는 검증 대상 8개 workload에서 output DRAM volume을 정확히 1회 기록한다. GB5의 legacy weight double-count 제거도 T5 `1.000`으로 확인됐다.

### GLB fill energy 스케일링

Eyeriss conv3 weight에서 다음 항등식이 성립한다.

- GLB→array weight read: `92,012,544 B`
- DRAM/multi-chip→GLB weight fill: `1,769,472 B`
- GLB weight access 합계: `93,782,016` byte accesses
- 설정 단가: `3 pJ/byte access`
- 출력 energy: `93,782,016 × 3 = 281,346,048 pJ`

따라서 fill energy를 일반 GLB access와 분리하고 datatype repetition으로 스케일한 뒤 합치는 방향은 정상이다.

### Timing 및 일반 회귀

| 항목 | 결과 |
| --- | --- |
| Gemmini RTL latency | MAPE `4.40%`, max `7.86%` |
| Eyeriss silicon latency | MAPE `4.26%`, max `6.39%` |
| Config validation | 203개 파일 통과 |
| Python unit tests | 5/5 통과 |
| ASan/UBSan smoke | Gemmini 및 Eyeriss 통과 |
| Patch whitespace 검사 | `git diff --check` 통과 |

## 잔여 문제

### P1 — GLB fill cycle이 layer timeline에서 누락됨

`multi_chip_t::account_descriptor_dense_distribution()`은 mc→GLB write cycle을 `global_buffer_t::fill_access_cycle`에 기록한다. 그러나 `global_buffer_t::modeled_elapsed_cycles()`는 `access_cycle`, `transfer_cycle`, `cycle_pe_array_global_buffer`만 확인하고 `fill_access_cycle`은 확인하지 않는다.

또한 `stats_t::update_stats()`가 GLB busy cycle, layer latency, static energy를 먼저 계산하고, `scale_serial_repetitions()`가 실행된 뒤에야 `merge_global_buffer_fill()`이 fill cycle을 최종 GLB access cycle에 합친다.

그 결과 Eyeriss의 지원 레이어 8개 모두 다음 불변식을 위반했다.

```text
busy_cycle_global_buffer >= max(access_cycle_global_buffer)
```

대표 사례:

| 레이어 | GLB busy cycle | 최종 max GLB access cycle | 차이 |
| --- | ---: | ---: | ---: |
| conv1 | 122,696,640 | 122,734,656 | -38,016 |
| conv3 | 92,052,480 | 93,782,016 | -1,729,536 |
| fc1 | 306,708,480 | 377,487,360 | -70,778,880 |
| fc3 | 33,280,000 | 40,960,000 | -7,680,000 |

현재 Eyeriss 검증 설정은 static energy가 0이고 DRAM이 전체 critical-path bottleneck이라 최종 latency/energy 표에서 영향이 숨겨진다. 하지만 GLB-bound workload 또는 non-zero leakage 설정에서는 latency와 static energy가 과소 계상될 수 있다.

관련 코드:

- `components/multi_chip.cc:414-445`
- `components/global_buffer.cc:1583-1593`
- `scheduler/stats.cc:418-491`
- `scheduler/stats.cc:760-885`

### P1 — 외부 Eyeriss traffic 기준 미달

`python3 validation/check_timing.py --check-traffic`은 exit code 1로 실패했다.

```text
Eyeriss traffic: max=669.91% (milestone-3 limit 15%)
```

레이어별 최신 비교:

| 레이어 | DRAM ratio | GLB ratio |
| --- | ---: | ---: |
| conv1 | 1.3× | 7.7× |
| conv2 | 1.5× | 2.6× |
| conv3 | 1.2× | 2.0× |
| conv4 | 1.7× | 2.6× |
| conv5 | 2.3× | 2.6× |

이는 기존 PA9로 기록된 batch filter sharing 및 psum spill/reload 모델 부재와 일치한다. T1~T6 통과는 내부 카운터 항등식을 보장하지만 실제 칩 traffic 정확도를 보장하지 않는다.

### P2 — Energy checker가 판정하지 않음

`validation/energy/check.py`는 E1~E5를 설명하고 결과표를 출력하지만 threshold assertion이나 실패 return code가 없다. 따라서 잘못된 energy 값도 파일 파싱만 가능하면 성공으로 종료된다.

독립 검산 결과 E1(ALU energy == MAC count)은 Eyeriss 8개 레이어 모두 성립했다. 그러나 문서의 E5인 `DRAM energy == serialized traffic × 200/word`는 packet padding이 존재하는 레이어에서 정확히 성립하지 않는다.

| 레이어 | DRAM energy / serialized-traffic energy 차이 |
| --- | ---: |
| conv1 | -0.485% |
| conv2 | 0.000% |
| conv3 | 0.000% |
| conv4 | -0.364% |
| conv5 | -0.438% |
| fc1 | -4.410% |
| fc2 | -4.408% |
| fc3 | -10.336% |

현재 모델은 DRAM access energy에는 logical data bytes를, link traffic에는 packet padding이 포함된 serialized bytes를 사용한다. 이 계약을 유지한다면 E5를 logical DRAM accesses 기준으로 고치고, link padding energy는 별도 interconnect 항목으로 검증해야 한다.

### P2 — DR6 반영 후 timing baseline 미갱신

`python3 validation/check_timing.py --check-baseline`은 다음 차이로 실패했다.

```text
conv1 measured DRAM = 6.532416 MB
conv1 baseline DRAM = 4.166976 MB
```

DR6 final output flush가 의도적으로 traffic을 변경했지만 `validation/phase3/npusim_baseline.csv`는 이전 값이다. P1 GLB timeline 문제를 먼저 해결한 후 최신 결과로 baseline을 갱신해야 한다.

### P3 — Output event counter의 방향성이 출력에서 불명확함

Output writeback transaction과 energy는 정상적으로 기록되지만 DRAM과 multi-chip의 `num_data_transfer[OUTPUT]`은 0이다. 헤더 계약상 이 카운터는 DRAM→multi-chip 및 multi-chip→GLB 방향만 세므로 동작 자체는 일관되지만, 결과 파일의 `# of data transfer` 라벨은 양방향 통계처럼 보일 수 있다. `read/load events`로 라벨을 명확히 하거나 read/write event를 분리하는 것이 좋다.

## 권장 수정 순서

1. GLB fill cycle을 datatype repetition 적용 후 busy-cycle 및 layer critical path에 포함하고, leakage 시간창도 같은 latency를 사용한다.
2. `validation/energy/check.py`에 E1~E5 assertion과 비정상 종료 코드를 추가한다.
3. E5를 logical memory access와 serialized link traffic으로 분리한다.
4. 수정 완료 후 Phase-3 NPUsim baseline을 갱신한다.
5. PA9의 batch filter reuse와 psum GLB spill/reload를 구현한 후 15% external traffic gate를 활성화한다.

## 생성된 최신 결과 위치

- Eyeriss timing/traffic: `result/eyeriss/alexnet/silicon/`
- Eyeriss energy: `result/eyeriss_energy/alexnet/silicon/`
- Gemmini timing/traffic: `result/gemmini/`

이 보고서를 작성하는 과정에서는 simulator 소스 코드를 변경하지 않았으며, 검증 실행으로 결과 파일만 재생성했다.
