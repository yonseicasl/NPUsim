# 2026-08-20 energy 재평가 (E20) 수정 기록

기준: [../../assessment/TIMING_SIMULATION_ENERGY_REEVALUATION_2026-08-20.md](../../assessment/TIMING_SIMULATION_ENERGY_REEVALUATION_2026-08-20.md)

권장 순서(Phase 1 → 6)를 따라 **E20-1, E20-2, E20-4, E20-6**을 해결했다. E20-3과 E20-5는 성격이
다르므로 아래 "남은 것"에 규모와 이유를 적었다.

## 0. 먼저 — 평가 §7이 지적한 오염 제거

문서 §7이 "같은 workspace에서 mutation campaign이 source를 변이하며 build를 반복해 binary와 result가
오염됐다"고 기록했다. **그 campaign이 이 세션의 것이었고 아직 돌고 있었다.** 결과를 수확한 뒤
(13 mutant: 4 killed / 9 survived) 중단하고, 소스를 stable backup에서 복원해 `diff -rq`로 동일함을
확인하고, 알려진 두 지점(`num_additions = leaves - 1`, `column_height = std::max`)이 원본인지
확인한 뒤에야 작업을 시작했다. 이후 build와 gate 실행은 직렬로만 했다.

## E20-1 + E20-2 (P0) — 미교정 결과의 absolute/wattage 자격 차단

### 먼저 필요했던 것: 활동량 counter 3개

`UNPRICED_ACTIVE`는 **event가 일어났는가**와 **가격이 선언됐는가**를 교차해야 판정된다. 그런데
energy 축은 두 경우 모두 0을 낸다. 세 event에 count가 아예 없어 먼저 만들었다:

| event | 추가한 counter |
| --- | --- |
| DRAM row activation | `dram_t::row_activation_events` |
| Format-IP payload / metadata | `pe_t::format_payload_events` / `format_metadata_events` |
| 배열 reduction addition | `pe_array_t::reduction_additions` |

### 판정 규칙

event는 **자기 counter가 0이 아니면 ACTIVE**, config가 그 key를 **어떤 값으로든(0 포함) 선언하면
PRICED**다. `key = 0`은 "이 설계는 여기에 비용을 두지 않는다"는 **진술**이고, key 부재는 진술이
아니다. PE local-buffer leakage만 조건부다 — config가 `static_energy`를 선언해 leakage를 scope에
넣었을 때만 active로 센다.

### 결과

`gemmini_cacti22`가 이 문제의 표본이었다. SRAM은 CACTI, DRAM은 DRAMsim3에서 왔으므로 이전의 모든
absolute 조건을 만족해 **418.685 mW를 보고**했다. 같은 report에 layer setup, accumulator reload
1,048,560 byte, spill 1,048,576 byte, final cast 4,096 byte가 **전부 0.00 pJ**로 찍혀 있었다.
산술은 맞았고 분자에 미선언 항이 6개 있었다.

| 완료조건 | 결과 |
| --- | --- |
| `gemmini_cacti22`가 wattage 거부 | ✅ 미과금 event 6개를 이름과 함께 나열 |
| 완전한 fixture만 wattage | ✅ `gemmini_power_calibrated` 9.394 mW (미과금 0) |
| component cost 누락 fixture 거부 | ✅ 신규 `gemmini_power_missing_dram` |

power 자격은 이제 다섯 조건을 모두 요구한다: single clock domain, absolute unit + provenance,
현재 precision의 calibrated compute cost, `undercounted() == 0`, `UNPRICED_ACTIVE` 없음.

새 gate [validation/unpriced](../../validation/unpriced/check.py) UP1~UP5. **UP4**가 핵심 구별을
못 박는다 — explicit 0은 플래그를 지우고 에너지에 기여하지 않으며, 부재는 검출된다. **UP2**는
priced key 6개를 하나씩 제거해 전부 검출되는지 본다. **UP3**은 비활성 feature(dense int8의 metadata,
bus array의 adder tree, row buffer 없는 DRAM)를 잡지 않는지 본다. network rollup 전파도 확인했다.

### Phase 4 gate의 방향 전환

P4-4는 `gemmini_cacti22`가 wattage를 **받아야** 한다고 주장했다. E20-1이 그것을 거부하므로 회귀가
아니라 **이전 판정이 틀렸던 것**이다. 이제 P4-4는 거부를 확인하고, **거부 이유가 올바른지**(unit
미선언이나 precision 미교정이 아니라 분자 undercount) 까지 검사한다.

## E20-4 (P1) — retained-output accumulator reload 경계

지적대로 `account_accumulator_reload()`가 `skip_transfer[OUTPUT]` 가드 **앞에서** 호출됐다. 물리적
transfer 분기 안으로 옮기고, retention을 `Accumulator retained events`로 보이게 했다. spill은 이미
물리적 write-back과 짝지어져 있어 경계가 맞았다.

**크기를 측정했다.** 기존 fixture에서는 pre/post가 1,048,560으로 동일했다 — 버그가 잠재적이었고,
그래서 accumulator gate에 pre-fix 형태를 주입해도 실패하지 않았다. E20-4가 요구한 retained mapping
(PE temporal level에 reduction 차원 C만) 을 만들자:

| | pre-fix | post-fix |
| --- | ---: | ---: |
| reload bytes | 1,048,512 | **16,384** |
| retained events | 16,127 | 16,127 |

**64배 과대 과금** — 일어나지 않은 read-back 1,032,128 byte에 대한 traffic과 energy다.

AC14~AC16 추가. **AC15**가 특히 구별력이 있다: retained/create는 fp32와 fp16에서 **동일**하고
(retention은 mapping 속성) reload와 spill만 정확히 2배다. retained pass를 reload로 청구하는 모델은
reload도 retention을 따라 움직이므로 걸린다.

## E20-6 (P2) — 2D mesh multicast

### 전제를 먼저 검증했다

인용된 생존자 문맥에는 `tallest_row = std::max(...)`와 `column_height[x] = std::max(...)`가 둘 다
있었고, `tallest_row`는 `(void)tallest_row;`로 폐기되는 **죽은 변수**다. 어느 쪽이 변이됐는지에 따라
결론이 달라지므로 각각 직접 시험했다 — **둘 다 생존했다.** 전자는 진짜 equivalent mutant이고
후자는 실제 검출력 공백이다. 즉 E20-6의 지적은 옳다.

원인도 확인했다: multicast fixture가 2×1과 1×4뿐이고, **1×4는 모든 chip이 row 0이라
`column_height`가 항상 0**이어서 max와 min이 같은 값을 낸다.

### 작업 중 잡은 제 실수

첫 2×2 mapping은 K를 CHIPS_X에, B를 CHIPS_Y에 두었다. 그러면 **모든 datatype이 chip에 분할되어
multicast 분기가 아예 실행되지 않는다** — 덮으려던 경로를 조용히 비활성화한 것이다. 탐침
(`tree_links = 100`)으로 발견했고, weight가 broadcast되도록 chip 차원에 batch만 두어 고쳤다.

### 결과

| fixture | baseline | `column_height` max→min | tree links |
| --- | ---: | ---: | --- |
| 1×4 (기존) | 5120 | 5120 (불변) | 4 → 4, column 항이 0 |
| **2×2** | 5120 | **2560** | 4 → 2 |
| **3-wide, 4 active** | 5120 | **3840** | 4 → 3 |

손 계산과 정확히 일치한다: 2×2는 `1 + 1 + (1+1) = 4`, 3-wide 4-active는 `1 + 2 + 1 = 4`.
NP5~NP7 추가. **NP5**는 "N개 router의 spanning tree + ingress = N link, grid shape 무관"이라는
성질을 고정한다(세 shape 모두 4). 변이 시 NP5·NP6이 정확한 link 수와 함께 실패한다.

## E20-5 (P1) — 실행 가능한 절반

전체는 외부 reference가 필요해 닫을 수 없지만(아래 참조), 두 항목은 지금 처리했다.

**항목 4의 두 번째 선택지** — shipped config는 `uncalibrated`로 두고 wattage를 금지한다. 8개 config
전부 이미 그렇게 동작하지만 **검사되지 않아** 조용히 회귀할 수 있었다. **PW8**로 고정했다.

**항목 5** — `dram_config` 계약을 하나로 정하고 report에 명시했다:

```
DRAM cost provenance  : the [dram] unit costs, verbatim. `dram_config` selects a DRAMsim3
                        device model in a -DDRAMSIM3 build only and derives no cost here
```

device 이름은 비용의 provenance가 아니다. 그래서 HBM2를 적은 config와 DDR3를 적은 config가 key가
같으면 동일하게 과금된다 — Phase 4가 기록한 그 사실이다.

## 회귀

Gate 24종 전부 통과, unit test 통과, 외부 reference 전부 불변:
Gemmini RTL 4.4%/7.9%, Eyeriss latency 4.3%/6.4%, SCALE-Sim PASSED, `check_timing` exit 0.

## 남은 것

**E20-3 (P1) — 외부 traffic 정합성.** psum GLB spill/reload source와 batch filter reuse를 **모델에
구현**해야 한다. 새 traffic source를 넣는 작업이고 Eyeriss DRAM/GLB 수치를 의도적으로 움직이므로,
외부 baseline 3종에 대한 영향 검토가 선행돼야 한다. 규모상 별도 착수가 맞다.

**E20-5 (P1) — 나머지 calibration 축.** MAC compute는 명시된 node의 합성 또는 검증 가능한 공개
수치가, NoP link는 pJ/traversal reference가 필요하다. Simba는 유료(HTTP 403)라 이 환경에서 얻을 수
없고, 기억에 의존한 수치를 reference로 쓰지 않는다는 원칙은 유지한다. 필요한 항목은
[validation/phase4/PROVENANCE.md](../../validation/phase4/PROVENANCE.md)에 적혀 있다.

---

# E20-3 (P1) — 외부 traffic 정합성

착수해 보니 평가 문서가 지목한 **두 원인이 모두 측정으로 반박됐다.** 실제 원인을 찾아 GLB 절반을
구현했고, DRAM 절반은 원인을 확정한 상태에서 멈췄다.

## 반박 1 — "psum GLB spill/reload traffic source 부재"

문서는 GLB 부족분의 원인을 source 부재로 지목했고, 내가 `phase3/compare.py`에 써둔 주석도
"모델 내부의 어떤 mapping으로도 칩 수치에 도달할 수 없다"고 단언하고 있었다. conv3의 mapping에서
output 인자를 array level → GLB level로 옮겨봤다 (**곱은 동일, 작업량 동일**, level만 이동):

| conv3 mapping | GLB SRAM traffic | computation cycle |
| --- | ---: | ---: |
| shipped | 9.2 MB | 3,833,856 |
| K를 GLB로 | 32.2 MB | 15,335,424 |
| K+P를 GLB로 | 101.2 MB | 199,360,512 |
| **칩 실측** | **50.2 MB** | **4,360,000** |

칩 수치는 모델 안의 mapping으로 **도달 가능하다.** psum GLB 경로는 존재하고 동작한다. 정확한
진술은 더 좁다: **traffic과 latency를 동시에 맞출 수 없다.** 두 주석 모두 고쳤고, 반박 근거를
`eyeriss_psumprobe*.cfg`에 재현 가능하게 남겼다.

## 진짜 원인과 구현

`equal_output_tile`은 배열의 output **tile 하나**를 GLB의 tile과 비교해 "같으니 움직일 필요 없다"고
판단했다. **그 위 loop가 여러 tile을 순회한다는 사실이 비교에 들어 있지 않았다.** conv3은 DRAM
level C=64이고 GLB가 312개 output tile을 순회하므로, 배열은 tile 하나를 들고 나머지를 전부 돌고
돌아온다 — 재방문마다 물리적 read-back과 write-out이다.

`mapping_table_t::reduction_tiled_above_array()`를 추가했다: C/R/S 중 하나라도 GLB·CHIPS·DRAM
level에서 factor > 1이면 retention은 성립하지 않는다. `pe_array_t::psum_retention_valid`로 layer마다
한 번 기록하고(`request_data()`에는 scheduler가 없다) 세 stationary 분기 **전부**에 적용한다.

### 헛짚은 두 번

1. **용량 기반 spill을 먼저 시도했다가 폐기했다.** conv3의 배열 psum tile은 용량의 20.6%로 충분히
   들어간다 — 용량은 구속 조건이 아니었다.
2. **output-stationary 분기에만 guard를 걸었다.** eyeriss는 `pe_stationary = weight_stationary`라
   else 분기를 탄다. 논리는 stationary type과 무관하므로 세 분기 전체로 옮겼다.

두 번 다 측정이 가설을 반박했다.

### 결과

| conv3 | before | after | 칩 |
| --- | ---: | ---: | ---: |
| GLB SRAM traffic | 9.2 MB | **40.9 MB** | 50.2 MB |
| 부족 배수 | 5.5× | **1.23×** | — |
| compute schedule 오차 | −2.4% | **−2.4%** | — |

**latency는 전혀 움직이지 않았다.** loop nest이 함의하는 66.5 MB와도 같은 규모다. 변경은
`equal_output_tile`이 잘못 참이 되던 경우에만 발화하므로 gemmini(이미 false, edge accumulation
경로)를 포함해 다른 config는 영향이 없다. Gate 24종·unit test 전부 통과, 외부 baseline 3종 불변.

## input halo 해결 — "batch filter reuse 미모델링" 가설의 종결

수정 전 DRAM input만 과대였고 weight는 이미 fetch-once였다. 따라서 batch/filter
reuse 문제가 아니라 P/Q sliding-window의 중첩 영역을 매 tile 재인출한 것이 원인이었다.

| layer | 수정 전 input | 수정 후 unique elements | 수정 후 DRAM txn (64b) | T3 |
| --- | ---: | ---: | ---: | ---: |
| conv1 | 3.31× | 629,244 | 157,311 | **1.00×** |
| conv2 | 5.06× | 380,928 | 95,232 | **1.00×** |
| conv3 | 2.60× | 230,400 | 57,600 | **1.00×** |
| conv4 | 2.62× | 345,600 | 86,400 | **1.00×** |
| conv5 | 2.62× | 345,600 | 86,400 | **1.00×** |

`mapping_table_t::input_halo_reuse()`가 padded P/Q·stride·filter를 포함한 전체
input union과 sliding ring working set을 계산한다. `global_buffer_t`는 separate/shared
partition, resident weight/output tile, double buffering을 포함해 실제 수용 가능 여부를
판정한다. 수용 가능할 때만 logical request를 유지한 채 중복 payload를 병합한다. DRAM
source access는 exact 정수 access 수에서 cycle/energy를 재구성하고, multi-chip과 GLB-fill은
같은 coalescing factor로 traffic/cycle/energy를 줄인다.

결과 파일은 replicated→unique element, required working-set bytes, DRAM 전후 transaction을
출력한다. `validation/traffic/check.py` T3/T11은 다섯 convolution에서 이 값이 exact union과
일치하고 halo 계약이 실제로 applied였는지를 함께 고정한다. conv2처럼 P와 Q를 동시에
tile하거나 conv1처럼 마지막 transaction tail이 있는 경우에도 단순 반복비가 아니라 exact
transaction 수를 사용한다.

legacy R/S factor가 1이 아닌 mapping은 P/Q sliding과 loop-order 의미가 다르므로 아직 이
계약을 적용하지 않는다. 검증되지 않은 reuse를 추정하지 않는 의도적인 보수 범위다.

## 남은 것

**외부 DRAM point-match.** input halo 이후 NPUsim은 모든 convolution에서 ideal dense union을
정확히 한 번 읽는다. 반면 RLC dense-equivalent chip 비교의 오차는 더 작은 방향(0.58~0.95×,
MAPE 23.8%)으로 남는다. 이는 halo 중복이 아니라 chip의 비공개 mapper/refetch scheduling과
Table V의 RLC dense-equivalent 변환에서 생기는 차이이며, 공개 데이터만으로 둘을 분리할 수 없다.
따라서 외부 DRAM 축은 계속 informational로 유지한다.

**GLB 나머지.** chip이 현재 mapping의 [완전보존 하한, literal 상한] 안에 있는 것은 1/5층이다.
나머지는 psum spilling이 GLB capacity가 아니라 tiling hierarchy를 따라 발생하는 문제가 남는다.

## 후속 마무리

- FUNCTIONAL 빌드의 retained-output 분기도 analytical 경로와 동일하게
  `accumulator_retained_events`를 기록한다. 실제 read-back이 없으므로 reload는 과금하지 않는다.
- `energy_units_t::describe()`는 layer 실행 전에 출력되므로 run-level unpriced event를 알 수 없다.
  따라서 기존의 `totals and power are meaningful` 단정을 제거하고
  `declared absolute-unit candidate; run-level completeness required`로 바꿨다. 최종 absolute/전력
  자격은 계속 `stats_t`의 UNDERCOUNT 및 UNPRICED_ACTIVE gate가 결정한다. `UP5`가 incomplete
  fixture에서 과도한 문구가 다시 나타나지 않는지 검사한다.
