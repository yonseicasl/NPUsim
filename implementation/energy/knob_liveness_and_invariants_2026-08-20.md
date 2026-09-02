# 2026-08-20 dead-knob sweep + metamorphic invariant 검증 추가

기존 gate 17종은 모두 **동작하는 축**을 검사한다: 양을 지목하고, 손으로 유도하고, 모델이 일치하는지
확인한다. 이 형태로는 이 프로젝트가 평가 회차마다 한 번씩 맞아온 결함을 **원리적으로** 잡을 수 없다 —
파싱되고 출력되지만 **아무도 소비하지 않는 설정**. 그런 knob은 어떻게 설정해도 report가 bit-identical
이므로, 동작하는 축을 아무리 많이 검사해도 드러나지 않는다.

지금까지 발견된 사례는 전부 손으로, 하나씩 찾은 것이다:

| knob | 증상 | 발견 |
| --- | --- | --- |
| `accumulator_format` | 파싱·출력되지만 spill은 OUTPUT datatype으로 크기가 정해져서 fp32/fp16이 동일한 traffic·energy | E4 |
| `row_miss_energy` | 내부 상태에만 누적되고 어디에도 출력되지 않음 | E1a |
| `noc_energy` | 모든 shipped config가 0이라 multicast multiplier(~17배 오차)가 한 번도 실행되지 않음 | E2/RE6 |

새 gate [validation/knobs/check.py](../../validation/knobs/check.py)가 이 부류를 **기계적으로**
찾는다: 선언된 설정을 하나씩 흔들고 재실행해 report에서 최소 한 숫자가 움직이는지 요구한다.

## 검사 (KN1~KN8)

| | 내용 |
| --- | --- |
| KN1 | **determinism** — 같은 config 두 번 실행이 동일한 report. 아래 모든 비교의 전제 |
| KN2 | **energy scale invariance** — 모든 energy unit cost ×10 → 보고된 모든 energy가 정확히 ×10, cycle·latency는 **하나도** 안 움직임. energy 회계는 unit cost에 선형·동차이므로, 배율이 안 맞는 energy는 하드코딩 상수에서 나온 것이고, 움직이는 latency는 energy cost가 timing path로 새는 것 |
| KN3 | **zero collapse** — 모든 energy cost를 0으로 두면 보고된 모든 energy가 **정확히** 0 (tolerance 없음) |
| KN4 | **dead-knob sweep** — 위 참조. 면제는 사유와, 해당되면 실제로 그 knob을 검사하는 gate 이름을 함께 선언 |
| KN5 | **energy/timing 분리** (knob 단위) — energy cost 하나만 올렸을 때 어떤 cycle·latency 필드도 움직이지 않음. KN2는 일괄 scaling에 대해 주장하고, 이건 하나씩 주장하므로 잘못 배선된 cost 하나가 다른 것들 뒤에 숨을 수 없다 |
| KN6 | **energy monotonicity** — energy unit cost를 올려서 어떤 energy도 내려가지 않음 |
| KN7 | **latency monotonicity** — cycle cost를 올려서 critical-path latency가 내려가지 않음 |
| KN8 | **면제 부패 방지** — 면제로 선언된 모든 knob이 실제로 주장대로 동작해야 함. 살아난 면제는 아무것도 주장하지 않는 checker이므로, 손으로 관리하는 allowlist가 흔히 겪는 부패를 막는다 |

Fixture는 `gemmini_all_knobs` — optional knob이 전부 non-zero인 그 config다. feature가 꺼진
knob은 재미없는 이유로 죽어 있으므로, 한 번의 sweep이 의미를 가지려면 이 fixture여야 한다.

## 결과

```
fixture gemmini_all_knobs: 251 reported numeric fields, 108 declared settings
  live knobs      : 75
  dead (declared) : 12
  invalid when perturbed (declared): 4
  non-numeric (covered by A/B fixtures): 17
```

**KN2가 통과한 것이 이번 라운드의 가장 안심되는 결과다** — 61개(현재 251개 필드 중 energy 축)
전부 정확히 10배, latency 이동 0건. RE7에서 하드코딩 `pJ` 15곳을 뜯어낸 직후라 특히 그렇다.

## sweep이 찾아낸 것

### 1. `t_wtr_cycle`이 어디에서도 검증되지 않고 있었다 → DR8로 해결

DRAM write→read bus turnaround는 `gemmini_dram_detail`과 `gemmini_dram_serial` **양쪽에서 죽어
있었다**. 모델 버그는 아니었다 — 구조는 맞다(read는 `false`, output write-back은 `true`로
`account_bus_turnaround`를 호출한다). 원인은 **DRAM gate가 `gemm_64x64x64`만 돌린다**는 것이었다.
그 schedule은 모든 load를 먼저, store를 나중에 내보내므로 bus가 read→write로 한 번 뒤집히고
write→read로는 **한 번도** 뒤집히지 않는다.

`gemm_512x512x512`에서는 tile이 여러 개라 store가 나중 load보다 앞선다. 항등식은 정확히 선형이다:

| `t_wtr_cycle` | Weight stream DRAM transfer cycles | delta |
| --- | --- | --- |
| 0 | 99,584 | — |
| 2 | 99,776 | +192 |
| 4 | 99,968 | +384 |
| 8 | 100,352 | +768 |

**96회의 write→read 뒤집힘**이고, 전부 뒤따르는 read를 하는 stream(여기서는 weight)에 과금된다.
[validation/dram/check.py](../../validation/dram/check.py) **DR8**이 이 항등식을 두 값에서
검사하고, 다른 stream이 움직이지 않는 것도 함께 확인한다.

### 2. multi-chip temporal buffer만 occupancy가 보고되지 않았다 → 메꿈

`multi_chip.memory_size`를 256 → **16** (16배 축소)으로 줄여도 실행되고 아무것도 바뀌지 않았다.
PE local buffer는 초과하면 에러를 낸다. 원인: `memory_size`는 capacity 검사에서 **소비되지만**
(`multi_chip.cc:302`) 그 occupancy가 어디에도 출력되지 않았다 — PE local buffer, PE-array temporal
buffer, GLB는 모두 출력하는데 multi-chip 계층만 없었다. E1a와 같은 부류(누적되지만 출력 안 됨)다.

capacity 검사가 이미 계산하는 값을 `buffer_utilization`으로 기록하고 report에 추가했다. 이제
`memory_size`가 live이므로 영구 면제가 실제 coverage로 바뀌었다.

### 3. 선언되었지만 죽어 있는 12개 — 사유 문서화

| knob | 사유 |
| --- | --- |
| `num_threads` | 시뮬레이션 병렬성. **죽어 있어야 정상**이며, 움직이면 비결정성이므로 `MUST_BE_INVARIANT`로 반대 방향을 주장한다 |
| `computation_energy` | 이 fixture에서 `mac_energy_int8_int8_fp32`가 우선한다(RE4 precedence). fallback 경로 자체는 PR5가 검사 |
| `bandwidth` (4개 section) | **bandwidth는 이 모델에서 rate limit이 아니다.** bitwidth가 없을 때 그것을 유도하는 source일 뿐이고(`derived_link_bitwidth`), 모든 section이 bitwidth를 명시한다. 사용자가 오해할 수 있는 지점이라 명시적으로 적었다 |
| `line_size`, `pe_array_read/write_cycle`, `pe_array_read/write_energy` | **PE-array temporal buffer를 실행하는 config가 repo에 하나도 없다** (access energy가 전부 0). 즉 이 5개는 어디에도 coverage가 없다 — 남은 실제 공백이며, 암묵적으로 두지 않고 표에 기록했다 |
| `t_wtr_cycle` | 위 참조. DR8이 다른 workload에서 검사 |

## sweep을 만들면서 잡은 자체 버그 3개

정직하게 적어둔다 — 셋 다 처음에는 "dead knob 발견"으로 보고됐다.

1. **`double_buffer`를 ×2로 흔들었다.** 이 knob은 depth가 아니라 `0 = single, 1 = double` 플래그다.
   1→2는 범위 밖이고 1처럼 동작한다. 올바른 perturbation은 toggle이며, A/B는 buffering gate가
   이미 `gemmini` vs `gemmini_single_buffer`로 검사하고 있었다.
2. **정수 설정에 `"2.0"`을 넣었다.** `get_setting<unsigned>`는 trailing `.0`에서 파싱 실패하고
   **조용히 기본값을 유지한다**. `number_of_macs`가 이 때문에 죽은 것처럼 보였다. (부수적으로:
   정수 knob에 소수를 쓴 config는 경고 없이 기본값으로 돌아간다. energy key는 RE5 schema가 막지만
   비-energy key는 막지 않는다 — 작지만 실재하는 robustness 공백이다.)
3. **필드 소멸을 변화로 세지 않았다.** `frequency`를 한 component만 올리면 single-clock-domain
   전제가 깨져 Power summary가 사라지는데(정상 동작), 숫자 이동만 보던 diff는 이를 "변화 없음"으로
   읽었다. 구조 변화도 liveness 증거다.

`gemmini_all_knobs`에 DRAM timing knob(`t_ras`/`t_rp`/`t_wtr`/`t_rtw`)을 추가했다. 그것들이 없으면
bank 병렬성이 감출 시간이 없어 `num_banks`가 완전히 죽는데, sweep이 처음 그렇게 보고했다.

## 검출력 확인

| 주입 | 결과 |
| --- | --- |
| energy를 하드코딩 상수로 계산 (`cast_bytes*0.25`) | **KN2 실패** — `Layer dynamic energy 14091076.4, expected 14100292.4` 등 4개 필드에서 지목 |
| live knob의 `get_setting` 제거 (`t_rtw_cycle`) | **KN4 실패** — `dram.t_rtw_cycle: 2 -> 4 changed NOTHING in the report` |

KN8은 구현 중 실제로 제 면제 표의 오류 7건을 잡아냈다(정수 보존 수정 후 buffer/array 크기 확대가
정상 실행되게 되자 `EXEMPT_INVALID` 항목들이 무효화됐다).

## 정정

작성 중 KN3 설명을 과장했다가 고쳤다. 원래 "additive term은 ratio test를 통과하지만 KN3만
잡는다"고 썼는데, KN2는 ratio가 아니라 **10×base와의 절대 비교**이므로 additive term도 잡는다.
KN3의 실제 고유 가치는 KN2의 출력 반올림 여유가 흡수할 만큼 작은 항(9배해도 여유 안에 남는 항)과,
tolerance 없이 성질을 못 박는 것이다 — 한 번의 추가 실행 값은 한다.

## 남긴 공백

- **PE-array temporal buffer를 실행하는 fixture가 없다** — knob 5개 무보증. mapping이 그 buffer를
  경유하도록 config를 하나 쓰면 닫힌다.
- **비-energy 정수 knob의 소수 값이 조용히 무시된다** (위 자체 버그 2). RE5 schema를 energy 밖으로
  확장하면 닫힌다.
- shipped architecture config 6개(`tpu`, `tpuv3`, `simba`, `maeri`, `eyerissv2`, `fsd`)는 여전히
  gate가 없다. 기본 workload가 repo에 없는 `weights/*.wgh`를 요구하고, gemmini의 mapping을 그대로
  복사할 수는 없다(simba는 chip당 4×4). 각 architecture용 합성 mapping을 써야 한다.

---

# 후속 1: PE-array temporal buffer coverage (2026-08-20)

sweep이 남긴 첫 숙제를 닫았다. 원인은 한 줄이었다.

```
#exist_temporal_buffer = 1      <- 48개 shipped config 전부, 두 section 모두 주석 처리
```

`exist_temporal_buffer`는 PE-array temporal buffer 경로 전체를 gate하고 `pe_array_t`에서 기본값이
**false**다. 모든 config가 주석 처리해 배포하므로 buffer는 어디서나 pass-through였고, 결과 파일의
access energy는 항상 0.00이었으며, cost knob 5개(`line_size`, `pe_array_read/write_cycle`,
`pe_array_read/write_energy`)는 무엇에 연결돼 있어도 무방했다.

## sweep 방법론의 공백 — 그래서 KN9를 추가했다

이건 dead-knob sweep에도 **안 보였다**. sweep은 config가 **선언한** 설정을 흔들기 때문에, 어떤
config도 선언하지 않는 knob은 사정거리 밖이다. 찾아낸 것은 보완적 스캔이었다: 코드가
`get_setting()`을 호출하는 모든 key **−** 어떤 config든 선언하는 모든 key.

**KN9**가 그 스캔을 영구화한다. 19개가 남고, 전부 분류를 요구한다:

| 분류 | 개수 | 예 |
| --- | --- | --- |
| alias (선언된 key로 fallback) | 4 | `nop`→`noc`, `num_processors`→`num_chips` |
| optional alternate (기본값 있음) | 6 | `parameter_order`, `data_format` |
| **UNCOVERED cost knob** | **9** | 아래 |

미분류 key가 나타나면 KN9가 실패한다. `undocumented_knob`을 하나 심어 확인했다.

### 기록된 backlog — coverage가 전무한 cost knob 9개

| knob | 내용 |
| --- | --- |
| `accumulator_spill_cycle` | RE1의 spill/reload **energy**는 AC1~AC9이 고정하지만 **time**은 어디에도 선언되지 않아, 모든 fixture에서 spill이 시간상 무료다 |
| `output_cast_cycle` | RE1 최종 cast도 동일 — energy는 AC6이 고정, time은 미선언 |
| `row_miss_cycle` | activation latency 직접 override. DR3~DR5는 `t_ras_cycle + t_rp_cycle`로 얻으므로 이 경로는 미실행 |
| `format_payload_energy`, `format_metadata_cycle`, `format_metadata_energy` | format fixture는 `format_payload_cycle`만 선언 → Format-IP energy가 모든 run에서 0이고 validation/format은 latency만 검사 |
| `lb_static_energy` | `static_energy`와 별개인 PE local buffer leakage. validation/leakage는 `static_energy`만 실행 |
| `adder_cycle`, `adder_energy` | `[adder_tree]` 전용 = `maeri.cfg`, gate가 전혀 없는 shipped config 6개 중 하나 |

## 새 gate: validation/array_buffer/check.py (AB1~AB6)

`gemmini_all_knobs`에서 양쪽 section의 flag를 켜고 read(0.4)/write(0.5) 단가를 구분했다.
`gemmini.cfg`는 pass-through 사례로 남는다.

단가 하나를 0으로 두어 **측정으로** 두 항을 분리하니 access 수가 그대로 나온다:

| | Input | Weight | Output |
| --- | --- | --- | --- |
| read accesses | 4,096 | 4,096 | **0** |
| write accesses | 4,096 | 4,096 | 65,536 |

**output의 read가 0인 것은 계약이다** — `edge_accumulation`에서는 array가 partial sum을 temporal
buffer에서 되읽지 않고 edge accumulator로 흘려보낸다. 쓰기는 한다.

| 검사 | 내용 |
| --- | --- |
| AB1 | 세 datatype 모두 non-zero — 아래 항등식이 공허하지 않음을 보장 |
| AB2 | 측정으로 분리한 read/write 항이 보고된 총합을 정확히 재구성 |
| AB3 | 두 방향이 독립 축. read 단가를 2배해도 **OUTPUT은 불변** |
| AB4 | `line_size` 2배 → 모든 access 수가 정확히 절반. array level에서 `line_size`를 실행하는 유일한 검사 |
| AB5 | flag를 지우면 모든 축이 **정확히 0** — knob들을 안 보이게 만들었던 그 기본값을 못 박는다 |
| AB6 | RE8 containment — 40,140.8이 PE array subtotal과 layer total에 각각 한 번, 다른 component는 불변 |

검출력: `pe_array.cc`의 `if(exist_temporal_buffer)` guard를 제거하니 AB5가 실패했다.

## 부수 발견 2개

**1. 두 `exist_temporal_buffer`의 기본값이 다르다.** `pe_array_t`는 `false`,
`multi_chip_t`는 **`true`**(`multi_chip.cc:18`). 그래서 multi-chip temporal buffer는 줄곧 켜져
있었고, flag를 1로 명시해도 아무 변화가 없어 처음엔 "consumed by nothing"으로 보였다. 명시적으로
0을 넣자 critical-path latency가 98,198 → 44,806으로 움직였다 — 살아 있다. 두 컴포넌트의 서로 다른
기본값은 버그는 아니지만 혼란의 소지가 있어 TOGGLE 주석에 적어뒀다.

**2. `stringstream >> bool`이 `"2"`에서 파싱 실패하고 조용히 기본값을 남긴다.** ×2 perturbation이
bool knob에서 dead로 오보하는 원인이었다. 앞서 정수 knob의 `"2.0"` 문제와 같은 부류다 — 두 경우 모두
config 값이 조용히 무시된다.

## 결과

| | 이전 | 이후 |
| --- | --- | --- |
| live knob | 75 | **82** |
| dead (사유 문서화) | 12 | **7** |

남은 7개: `num_threads`(죽어 있어야 정상), `computation_energy`(RE4 precedence),
`bandwidth` ×4(이 모델에서 rate limit이 아님), `t_wtr_cycle`(DR8이 다른 workload에서 검사).

---

# 후속 2: RE1의 latency 절반 (2026-08-20)

KN9 backlog의 첫 두 항목을 닫았다. RE1은 accumulator spill/reload와 최종 cast를 각자의 boundary로
분리하고 **energy**를 AC1~AC9로 고정했다. 그 **time**은 모델에 있었지만 —
`accumulator_spill_cycle`이 Format-IP cycle로, `output_cast_cycle`이 multi-chip write-back으로 —
**어떤 config도 두 key를 선언하지 않아** 모든 fixture에서 spill과 cast가 wall-clock상 무료였다.
값 검사로도, dead-knob sweep으로도 안 보인다(sweep은 선언된 것만 흔든다). KN9가 찾아냈다.

## 먼저 정리한 것 — RE1이 남긴 잔재 세 벌

`output_cast_cycle`/`output_cast_energy`를 **세 컴포넌트가 읽고 있었다**. RE1의 cast 경계는 세 번
시도했고(PE `write_output` → 262,144, GLB readout → 16,384, off-chip store → 4,096), 각 시도가
자기 컴포넌트에 unit-cost 쌍을 남겼다:

| 위치 | 상태 |
| --- | --- |
| `multi_chip_t` | **살아 있는 유일한 경로** |
| `pe_t` | `account_output_cast()`가 정의만 되고 **한 번도 호출되지 않음**. 그 함수만이 `u_output_cast_*`를 쓰고, `output_cast_bytes`/`output_cast_energy` 멤버도 아무도 읽지 않음 |
| `global_buffer_t` | 두 vector를 config에서 **파싱만 하고 어디서도 사용하지 않음** (두 생성자 모두) |

`accumulator_format`이 파싱·출력되고 아무도 안 쓴 것(E4)과 정확히 같은 부류이고, 이번엔 **내가**
남긴 것이다. 죽은 쪽 전부 제거했다.

## 관측성 — cast time은 report에 없었다

`output_cast_cycle`은 `write_back_cycle`에 누적되고, 그것은 `modeled_elapsed_cycles()`에서
access/transfer 축과 **max**로만 합쳐진다. 즉 그 축들보다 싼 cast는 critical path에서 전혀 보이지
않는다. 검증할 수가 없으므로 `Output cast cycle`을 별도 counter로 분리해 report에 추가했고, 계약도
함께 적었다:

```
Output cast cycle     :         2048.0 cycles  [enters the fabric's busy time as a MAX
                                                against its access/transfer axes, not a
                                                serial addition]
```

## 항등식 — 측정으로 확인

`gemmini_accum_fp32`/`fp16`에 `accumulator_spill_cycle = 0:0:0.25`,
`output_cast_cycle = 0:0:0.5`를 선언했다.

| 실험 | Format-IP cycle (Output) | Output cast cycle | edge acc energy |
| --- | --- | --- | --- |
| fp32 baseline | 2048.0 | 2048.0 | 1,048,568 |
| fp16 baseline | **1024.0** (정확히 절반) | **2048.0** (불변) | 524,284 |
| spill_cycle ×2 | 4096.0 | 2048.0 | 1,048,568 |
| cast_cycle ×2 | 2048.0 | 4096.0 | 1,048,568 |
| spill_cycle 미선언 | **0.0** | 2048.0 | 1,048,568 |
| cast_cycle 미선언 | 2048.0 | **0.0** | 1,048,568 |
| spill_**energy** ×2 | 2048.0 | 2048.0 | **2,097,136** |

`Format-IP cycle`은 PE 간 **max** × repetitions이고, cast cycle은 `4096 bytes × 0.5 = 2048`로 정확하다.

| 검사 | 내용 |
| --- | --- |
| AC10 | spill/reload TIME이 accumulator width를 따름 — fp32가 fp16의 정확히 2배. AC1/AC2의 시간 거울 |
| AC11 | cast TIME = cast bytes × `output_cast_cycle`, accumulator format에 **불변**. AC3/AC7의 시간 거울 |
| AC12 | 각 cycle knob이 **자기 축만** 선형으로 움직이고, 미선언이면 정확히 0 — 둘을 숨기던 그 상태 |
| AC13 | **energy/time 독립** — energy 단가를 2배해도 cycle 축이 안 움직이고, cycle 단가를 2배해도 energy 축이 안 움직인다. RE1이 네 event를 각자의 boundary로 나눈 것의 두 번째 차원 |

검출력:

| 주입 | 결과 |
| --- | --- |
| spill을 OUTPUT datatype으로 sizing (pre-E4 버그) | **AC10 실패** (`512.0 (fp32) != 2x 512.0 (fp16)`), AC1도 동시에 실패 |
| cast time을 MAC issue당 과금 (64배) | **AC11 실패** (`131072.0 != hand 2048.0`) |

## KN9 backlog: 9 → 7

닫힌 것: `accumulator_spill_cycle`, `output_cast_cycle`.
남은 것: `row_miss_cycle`, `format_payload_energy`, `format_metadata_cycle`,
`format_metadata_energy`, `lb_static_energy`, `adder_cycle`, `adder_energy`.

---

# 후속 3: Format-IP metadata와 energy (2026-08-20)

KN9 backlog의 format 3개를 닫았다. `format_payload_energy`, `format_metadata_cycle`,
`format_metadata_energy`가 코드에 읽히지만 어떤 config도 선언하지 않아 **Format-IP energy가 모든
run에서 0**이었고 block-scale metadata stream은 한 번도 과금되지 않았다. FM1~FM6은 payload
**latency**만 검증한다.

## 먼저 막힌 것 — 내 RE5 schema가 불완전했다

새 fixture를 쓰자마자 거부당했다:

```
Error: invalid energy unit cost: [systolic_array] format_metadata_energy is not a declared
energy key; an unrecognized energy key silently leaves its component at zero cost
```

RE5 schema 표를 **shipped config가 선언한 key 조사**로 만들었기 때문에, 코드가 읽지만 아무 config도
설정한 적 없는 key가 전부 빠져 있었다 — `format_payload_energy`, `format_metadata_energy`,
`lb_static_energy`, `adder_energy`, 그리고 `nop_energy` alias. 새 fixture가 그중 하나를 선언하려
하자 **오타로 거부**했다. schema는 config가 쓰는 것이 아니라 **코드가 읽는 것**에 대해 완전해야 한다.

**KN10**이 그 교차 검사를 영구화한다: `utils/energy_units.cc`의 표를 파싱해, 코드가 읽는 key 중
이름에 `energy`가 든 것이 전부 표에 있는지 확인한다. 표에서 한 항목을 지워 확인했다.

## fixture 설계 — payload를 고정하고 metadata만 켠다

**mxfp8의 payload는 8비트로 int8과 동일하고** block scale metadata만 추가된다. 그래서
`gemmini_format_mxfp`(operand만 mxfp8, output은 int8 유지)는 `gemmini_format`과 **metadata stream
하나만** 다르다 — payload traffic이 그대로이므로 두 항이 깔끔히 분리된다.

단가 하나씩 0으로 두어 분리 측정:

| | Input | Weight | Output |
| --- | --- | --- | --- |
| cycle: payload 항 | 4096 | 64 | 4096 |
| cycle: metadata 항 | 3072 | 48 | **0** |
| cycle: baseline | **7168** | **112** | **4096** |
| energy: payload 항 | 131072 | 2048 | 131072 |
| energy: metadata 항 | 65536 | 1024 | **0** |
| energy: baseline | **196608** | **3072** | **131072** |

int8 output은 block scale이 없으므로 metadata가 정확히 0이다. 네 단가를 전부 0으로 두면 두 축 모두 0.

| 검사 | 내용 |
| --- | --- |
| FM7 | mxfp8 operand는 metadata 항을 갖고, int8 output은 갖지 않으며, all-int8 fixture는 energy가 0 |
| FM8 | payload-only + metadata-only가 baseline을 **정확히** 재구성 (cycle·energy 양쪽, 세 stream 모두) |
| FM9 | 네 단가가 각자 자기 항에만 선형 |
| FM10 | cycle 단가는 energy 축을 안 움직이고 그 역도 성립. Format-IP energy 330,752이 MAC subtotal에 정확히 한 번 (RE8) |
| FM11 | **granularity 계약을 명시** — 아래 참조 |

## FM11 — 발견한 granularity 사실을 계약으로 못 박았다

block이 32 element인데 **metadata transaction 수가 payload와 같다**(input 262,144 / 262,144).
MAC tile 단위에서 둘 다 1 transaction으로 올림되기 때문이다.

이것을 버그로 보지 않았다. format IP는 block의 일부라도 디코드하려면 그 block의 scale을 가져와야
하므로 **tile당 scale 1회 접근**이 물리적으로 맞다 — 즉 이 축은 volume 비율이 아니라 **access
횟수**다. 다만 report만 보면 metadata 용량이 payload와 같다고 오해할 수 있으므로, FM11이 그 등식을
**선언된 계약**으로 고정하고 한계도 기록한다: 이 축은 payload:metadata의 32:1 용량비를 표현할 수
없다. RE1이 accumulator에서 같은 문제를 만나 per-byte 과금으로 옮긴 것과 동일한 구조다.

검출력:

| 주입 | 결과 |
| --- | --- |
| metadata를 payload 단가로 과금 (하나의 blended 양) | **FM7·FM11 실패** (metadata events 0 != payload events 524288) |
| RE5 schema에서 key 하나 삭제 | **KN10 실패** |

## KN9 backlog: 7 → 4

남은 것: `row_miss_cycle`, `lb_static_energy`(둘 다 기존 경로의 **대체 경로**라 "무엇을 하기로
되어 있는가"를 먼저 정해야 한다), `adder_cycle`, `adder_energy`(`maeri.cfg` 전용 — gate가 전혀 없는
shipped config 6개 문제와 묶여 있다).

---

# 후속 4~7: shipped architecture coverage, 엄격한 parser, hierarchy 보존, 대체 경로 knob

네 항목을 순서대로 진행했다. 결과부터: **KN9 backlog가 9 → 0이 되었다.** 남은 10개는 alias 4개와
기본값 있는 optional alternate 6개뿐이다.

## 4-1. shipped architecture 6개에 coverage 확보

`tpu`, `tpuv3`, `simba`, `maeri`, `eyerissv2`, `fsd`가 validation 파일에서 참조 0건이었다. 원인은
데이터 의존성이다 — 기본 mapping이 실제 network(resnet50, bert…)를 쓰고 그 weight blob이 repo에
없다. 각 architecture의 array·chip grid에 맞춘 **합성 GEMM mapping**을 써서 해결했다.

> 이 mapping들은 **coverage용이고 성능용이 아니다.** spatial factor만 array에 올리고 나머지는 DRAM에
> 두는데, 1-deep systolic PE buffer에 맞는 유일한 형태다. 그래서 예컨대 tpu의 critical path는 compute
> schedule의 약 26만 배다(batch마다 전부 refetch). "이 architecture가 자기모순 없이 맞는 계산을
> 하는가"로 읽어야 하고, mapping 품질 결과로 읽어선 안 된다.

실행시키는 과정에서 **실제 결함 2개**가 나왔다:

1. **multi-chip capacity 검사가 GLB bypass를 무시했다.** `global_buffer_t`는 bypass된 datatype을
   자기 검사에서 면제하고, `multi_chip_t`도 delivery traffic 회계에서는 이미
   `chips[i]->bypass[type]`를 참조한다 — capacity 검사만 안 했다. 그래서 `eyerissv2.cfg`(weight를
   ingress에서 spad로 직접 스트리밍, 양쪽 레벨 `weight_size = 0`)는 **실행 자체가 불가능**했다.
2. **`multi_chip_t`는 buffer 크기를 `unsigned`로, `global_buffer_t`는 같은 key를 `double`로**
   저장했다. eyerissv2는 `input_size = 4.5`(KB)를 선언하는데 unsigned 파싱이 실패해 **조용히 0**이
   되었고, 모든 tile이 그 0을 "초과"했다. 아무 gate도 이 config를 돌리지 않아 아무도 몰랐다.

capacity 오류 메시지도 어느 datatype이 얼마나 초과했는지 말하지 않아서 함께 고쳤다.

새 gate [validation/architectures/check.py](../../validation/architectures/check.py) SM1~SM10.
핵심은 **SM2 — 7개 architecture 전부 정확히 262,144 MAC**이다. mapping은 일을 어떻게 펼칠지만
정하고 얼마나 할지는 정하지 않으므로, workload를 정확히 덮지 못한 mapping이 여기서만 드러난다.

| architecture | MACs | sched | critical | PE% | chip% |
| --- | --- | --- | --- | --- | --- |
| tpu | 262144 | 127 | 33,185,950 | 6.2 | 100.0 |
| tpuv3 | 262144 | 127 | 12,468,606 | 6.2 | 100.0 |
| simba | 262144 | 256 | 239,621 | 100.0 | 11.1 |
| maeri | 262144 | 4096 | 3,173,378 | 38.1 | 100.0 |
| eyerissv2 | 262144 | 4096 | 1,306,635 | 33.3 | 50.0 |
| fsd | 262144 | 64 | 159,189 | 22.2 | 100.0 |
| gemmini | 262144 | 3518 | 59,010 | 100.0 | 100.0 |

**SM8~SM10 — MAERI의 adder tree.** maeri의 mapping에서 C를 **spatial**로 올렸다(C가 temporal이면
fan-in이 1이고 경로 전체가 early-return한다). 항등식이 정확하다:

- energy = distinct outputs(4096) × (fan_in−1 = 63) × 단가 = **154,828.8**
- cycle delta = 262,144 write-back × 6 tree level × 단가

`adder_tree_reduction_cost`의 N−1은 **검증이 전무했다**. `validation/reduction` RD2는 `utils/pe_lane.cc`의
per-MAC lane tree를 담당하는 **별개 구현**이라, N을 주입해도 RD2는 통과하고 SM9만 실패한다.
`transfer_cycle`에 tree fill이 **element마다** 더해지는 것(amortize 아님)은 상한으로 기록했다 —
systolic fold bubble 때와 같은 처리다.

## 4-2. 선언되었는데 쓸 수 없는 값은 거부한다

`section_config_t::get_setting()`은 값이 요청된 C++ 타입으로 파싱되지 않으면 `false`를 반환했고,
모든 호출자는 `false`를 "미선언"으로 읽는다 — 즉 **선언한 값이 조용히 사라지고 기본값이 남았다.**
세 경로 모두 하드 오류로 바꿨고, section·key·값을 이름으로 말한다.

바꾸자마자 **또 하나의 실제 버그**가 나왔다: `maeri.cfg`의 `[dram] transfer_cycle = 2.75:2.75:2.75`.
모델은 이 key를 scalar로 읽으므로 값이 통째로 무시되어 그 architecture의 DRAM transfer 비용이
0이었다. scalar `2.75`로 고쳤다.

`std::string` 설정도 전체 값(공백 포함)을 받도록 오버로드를 추가했다 — 이전에는 첫 토큰에서 멈추고
eof 검사에 실패해 여러 단어짜리 `energy_reference`가 조용히 버려졌을 것이다.

`get_vector_setting`은 실패 시 vector를 **절반만 갱신한 채** 포기할 수 있었다(선언값과 기본값이 섞인
상태). 이제 완성된 임시 vector를 만들어 한 번에 대입한다.

hard error는 in-process 단위 테스트로 관측할 수 없어, 세 거부 경우를 subprocess 기반인 knob gate로
옮겼다 — **KN11**. 단위 테스트에는 이유를 남겼다.

**KN10**도 이 라운드에서 추가했다: RE5 schema 표를 *shipped config가 선언한 key* 조사로 만들었기
때문에 코드가 읽지만 아무 config도 설정한 적 없는 energy key 5개가 빠져 있었고, 새 fixture가 그중
하나를 선언하려 하자 **오타로 거부**했다. schema는 config가 쓰는 것이 아니라 코드가 읽는 것에 대해
완전해야 한다.

## 4-3. hierarchy 보존 — 말할 수 있는 것과 없는 것

새 gate [validation/hierarchy/check.py](../../validation/hierarchy/check.py) HC1~HC4, 7개
architecture 전부에 적용.

**먼저 유도되지 않는 것을 분명히 했다.** 5개 경계에 걸친 **byte 단위 보존 법칙은 현재 report에서
나오지 않는다.** transaction이 경계별 link 폭 단위로 세어지고, 세밀한 granularity에서는 한
transaction이 link word보다 훨씬 적게 나른다 — MAC↔LB는 128비트 link에서 262,144 transaction으로
262,144 byte를 나르므로 폭을 곱하면 byte를 16배 부풀린다. FM11이 metadata에 대해, RE1이
accumulator에 대해 기록한 것과 같은 한계다.

**output request 사슬도 단순한 단조성이 없다.** 하나가 아니라 7개 architecture를 측정해서 가설 세
개가 연달아 반박됐다:

| 가설 | 반박 |
| --- | --- |
| array→GLB out == GLB→chip out − 1 | eyerissv2 (32,760 vs 2,048) |
| 바깥으로 단조 비증가 | tpu (63 < 64) |
| edge_accumulation이 켜졌을 때 정확히 0 | simba, fsd (꺼져 있는데 0) |

그래서 output 사슬은 **미해결 문제로 기록**했고 주장하지 않았다. 기여 요인 하나(`equal_output_tile`)는
`Output tile residency`로 report에 노출했다.

| 검사 | 내용 |
| --- | --- |
| HC1 | **load request 단조성** — input/weight 요청 수가 바깥으로 비증가. load는 miss일 때만 바깥으로 전파하므로, 바깥이 더 많으면 traffic이 무에서 생긴 것 (missing-source의 형태) |
| HC2 | DRAM traffic ≥ distinct tensor volume (모든 architecture) |
| HC3 | DRAM traffic이 volume의 **정수배** — mapping factor가 차원을 나누므로(SM2가 고정) 부분 tensor가 건너가면 안 된다 |
| HC4 | output tile residency가 report에 명시됨 |

측정된 load 사슬(input, 안쪽→바깥): gemmini `[2048, 4, 4, 4]`, tpu `[266240, 64, 64, 64]`,
simba `[256, 8, 8, 2]`, maeri `[262208, 4096, 4096, 4096]`, eyerissv2 `[131104, 32768, 32768, 4096]`,
fsd `[8192, 2, 2, 1]`.

## 4-4. 대체 경로 knob 2개 — 제 추측이 틀렸다

이전 보고에서 "`row_miss_cycle`은 `t_ras+t_rp`와 중복이라 제거가 맞아 보인다"고 적었다. **코드를
읽으니 틀렸다.**

- **`row_miss_cycle`**은 중복이 아니라 **대안 flat 모델**이다:
  `row_activation_cycle = (t_ras > 0 && t_rp > 0) ? t_ras + t_rp : u_row_miss_cycle`.
  `describe_timing_model()`이 이미 두 모델을 구분해 출력한다("JEDEC tRC=tRAS+tRP" vs "flat row-miss
  cost"). 모든 config가 row buffer를 끄거나 tRAS/tRP를 선언해서 flat 분기만 실행되지 않았다.
  → `gemmini_dram_flat` + **DR9**(분기가 선택·명시되고 비용에 선형, busiest bank 4 activation).
- **`lb_static_energy`**도 대안이 아니라 **가산 항**이다:
  `static_energy[type] = (u_static_energy[type] + u_lb_static_energy[type]) * elapsed_cycles`.
  → `gemmini_lb_leakage` + **LK6**(0.03을 0.01 위에 더하면 PE leakage가 정확히 4배, latency와
  dynamic energy는 불변).

LK6은 처음에 3.954배로 실패했다. layer static total에는 `pe_array_static_energy`처럼 이 두 key와
무관한 leakage가 섞여 있어서다 — **PE row**를 비교해야 정확히 4배가 된다.

## Gate 수

| | 시작 | 현재 |
| --- | --- | --- |
| Gate | 17종 | **21종** (`knobs`, `array_buffer`, `architectures`, `hierarchy` 신규) |
| KN9 UNCOVERED backlog | 9 | **0** |
