# 2026-08-19 timing simulation energy·power 재평가 수정 기록 (RE1~RE8)

기준 문서: [../../assessment/TIMING_SIMULATION_ENERGY_REEVALUATION_2026-08-19.md](../../assessment/TIMING_SIMULATION_ENERGY_REEVALUATION_2026-08-19.md)

권장 순서(RE1 → RE2 → RE3+RE5+RE8 → RE4+RE6 → RE7)를 그대로 따랐다. **Phase 4(외부 traffic 및
absolute energy calibration)만 제외**했다 — measured component energy reference가 선행되어야
하므로 코드로 닫을 수 없다.

각 issue를 닫은 직후 검증을 붙였고, **모두 negative test로 검출력을 확인**했다: 수정 전 형태를
코드에 다시 주입해 새 checker가 실제로 실패하는지 확인한 뒤 되돌렸다.

## RE1 — accumulator spill / reload / final cast event 분리

`account_accumulator_*`가 세 event를 각자의 boundary에서 만들도록 바꿨고, **최종 cast를
`multi_chip_t::account_output_writeback_to_dram()`에서 과금**한다. Ownership은
`edge_accumulation`을 따른다(PE 소유 vs PE array edge 소유).

Scaling을 uniform `m_repetitions`에서 **OUTPUT datatype repetition**으로 바꾼 것이 핵심이다.
그 전에는 cast가 최종 output element 수가 아니라 GLB pass 수만큼 부풀었다.

검증 — [validation/accumulator/check.py](../../validation/accumulator/check.py) AC1~AC9. 절대값을
고정하는 것이 **AC6**이다: final cast = 4,096 byte = 최종 output element 수. 잘못된 두 형태가
각각 다른 값을 주므로 구분된다.

| 형태 | cast bytes |
| --- | --- |
| PE `write_output` (reduction pass마다) | 262,144 |
| GLB readout (GLB pass마다) | 16,384 |
| **최종 off-chip store (정답)** | **4,096** |

Negative test: uniform `m_repetitions` scaling을 다시 주입하니 AC6이 fp32/fp16 양쪽에서
`16384 != 4096`으로 실패했다.

## RE2 — uncalibrated config가 absolute pJ / mW로 인정되던 문제

이것이 이번 라운드에서 가장 실질적인 결함이었다. `energy_unit` 기본값이 pJ였기 때문에 **아무
것도 선언하지 않은 `gemmini` config가 4.273 mW를 출력**했다. 같은 report가 바로 위에서
`MAC energy basis: UNCALIBRATED`, `provenance: not declared`라고 말하면서도 그랬다. 산술은
맞았고 **자격(qualification)이 틀렸다.**

- `energy_unit_t::UNSPECIFIED`를 기본값으로 두었다 (pJ 아님).
- `is_absolute()`가 **`energy_unit = pJ`와 non-empty `energy_reference`를 모두** 요구한다.
  단위만 선언한 pJ는 검증 불가능한 주장이다.
- Absolute power는 추가로 **현재 precision의 compute cost가 calibrated**일 것을 요구한다.
- 조건 미달이면 total에 `(ESTIMATED: not calibrated for absolute comparison)`을 붙이고 power는
  이유를 명시한 `unsupported`로 둔다.

Fixture 4종을 추가했고 power checker의 absolute fixture를 명시적 calibrated-pJ로 옮겼다:

| fixture | 상태 | power |
| --- | --- | --- |
| `gemmini_power_calibrated` | pJ + provenance + per-precision MAC energy | 5.833 mW |
| `gemmini_power_no_provenance` | pJ, reference 없음 | unsupported |
| `gemmini_power_precision_fallback` | pJ + provenance, precision fallback | unsupported |
| `eyeriss_energy` | normalized | unsupported |
| `gemmini` | 미선언 | unsupported |

검증 — [validation/power/check.py](../../validation/power/check.py) **PW7**. Negative test로
pre-RE2 형태(기본값 pJ + unit-only 판정 + precision gate 없음)를 주입하니 7건이 실패했고, 그
중 첫 줄이 정확히 이 결함이었다: `gemmini_power_no_provenance: reports 4.273 mW without
qualifying for absolute power`.

## RE3 — layer setup event가 unit energy에 종속됐던 문제

Setup event count를 `u_layer_setup_energy > 0`이 아니라 **`u_layer_setup_cycle > 0`(실제 실행)**
에서 만든다. 역방향(energy만 양수, cycle 0)은 event source가 없으므로 config load에서 거부한다.

수정 전 gemmini는 2270 cycle을 지불하면서 `0.00 pJ over 0 setup event(s)`를 출력했다 — "setup이
없었다"로 읽힌다. 이제:

```text
Layer-setup energy    :           0.00 uncalibrated over 1 setup event(s)
                        [unit cost UNCALIBRATED: the setup executes for 2270 cycle(s)
                         but no layer_setup_energy is declared]
```

검증 — [validation/energy_schema/check.py](../../validation/energy_schema/check.py) ES6/ES7,
[validation/fold_energy/check.py](../../validation/fold_energy/check.py) FS6, 그리고
`unittest/run_validation.sh`의 `setup-energy-no-event.cfg` 거부.

## RE5 — energy config 검증을 schema 수준으로

`validate_energy_settings()`가 값 위생만 보던 것을 **schema 검증**으로 올렸다. 모든 energy key를
arity와 함께 선언하고(`energy_key_schema()`), 밖에 있는 것은 **가장 가까운 known key를 알려주는
error**로 만든다.

새로 잡히는 것: key 오타, per-datatype vector 길이 오류, 빈 중간 field, scalar 자리에 vector.
275개 config 전부 통과한다.

`Uncosted components`를 **결과값이 아니라 선언 상태**에서 계산한다. 0이라는 숫자 하나가 서로
전혀 다른 네 상황을 뭉뚱그리고 있었다:

| 상태 | 의미 | fixture |
| --- | --- | --- |
| `NOT MODELED` | energy key를 하나도 선언하지 않음 | `gemmini_cost_not_modeled` |
| `PARTIAL (UNDERCOUNT)` | 일부만 선언 — 오타가 이 모습이다 | `gemmini_cost_partial` |
| `modeled zero` | 전부 선언하고 전부 0 — 의도된 free | `gemmini` (PE array, NoP) |
| `NO ACTIVITY` | 전부 선언되고 non-zero인데 layer가 안 건드림 | `gemmini_cost_no_activity` |

Optional-feature cost(`layer_setup_energy`, `weight_fold_fill_energy`,
`accumulator_spill_energy`, `output_cast_energy`, `row_miss_energy`, `mac_energy_*`)는 부재가
gap이 아니지만 **선언된 non-zero는 component를 costed로 만든다.** 이 구분이 없으면 setup만
과금한 config가 7.50을 출력하면서 "modeled zero"로 표시됐다.

검증 — ES1~ES5. Negative test로 결과 기반 판정을 되돌리니 12건이 실패했다.

## RE8 — non-zero checker에 event 의미와 containment 추가

| 공백 | 대응 |
| --- | --- |
| accumulator checker가 최종 output 수를 대조 안 함 | **AC6** (RE1) |
| reduction energy의 layer/network containment 미확인 | **RD5** |
| energy summary checker 주 fixture에서 optional cost가 전부 0 | **E11** + `gemmini_all_knobs` |
| power checker가 default-pJ fixture를 absolute로 사용 | **PW7** (RE2) |
| fold/setup checker가 uncalibrated event count 미확인 | **FS6** (RE3) |

**E11**이 특히 중요하다. `eyeriss_energy`에서는 reduction/fold/setup/accumulator/cast가 모두
0이므로, 어떤 축을 subtotal에서 빼먹어도 총합 항등식이 통과한다 — 0을 맞는 row에 더하는 것과
틀린 row에 더하는 것이 구별되지 않는다. `gemmini_all_knobs`는 이 축들을 **동시에 non-zero**로
만들고, checker가 fixture를 직접 실행한다(stale 결과 파일을 읽으면 이전 build에 대해 통과한다 —
실제로 첫 시도에서 이 함정에 걸렸다).

Negative test: PE array subtotal에서 `layer_setup_energy`를 빼니 E11이 정확히 500 pJ 차이로
실패하고, eyeriss_energy의 E7은 여전히 통과했다 — RE8이 지적한 공백 그대로다.

## RE4 — MAC precision key에 accumulator precision 포함

Key를 `mac_energy_<input>_<weight>_<accumulator>`로 확장했다. 2-part key는 **partial
calibration**으로 계속 인정한다(bare scalar보다 낫다) 그러나 **absolute power 자격은 주지
않는다** — 모든 accumulator format이 같은 값을 공유하기 때문이다.

INT8×INT8 A/B (`gemmini_macacc_fp32` / `gemmini_macacc_fp16`): MAC count와 operand traffic이
bit-identical이고 compute energy만 146,801 vs 89,129으로 갈린다. 2-part key에서는 이 차이가
표현되지 않았다.

`gemmini_macacc_operand_only`(accumulator fp16, operand key만 선언)는 그 cost를 쓰면서
`(operands only; accumulator fp16 UNCALIBRATED)`로 보고하고 power는 unsupported다.

검증 — [validation/precision/check.py](../../validation/precision/check.py) PR6/PR7,
power checker의 PW7. Negative test로 2-part 조회를 되돌리니 PR6/PR7이 실패했다.

## RE6 — mesh multicast link traversal 계약 확정

**계약을 정의했다**: `noc_energy`는 **array 내부 router-to-router link 1개의 1회 traversal**
단가다. GLB↔array attach link는 GLB의 `transfer_energy`가, 각 PE의 ejection은 PE local buffer의
write가 이미 과금한다.

그 계약 아래서 `N`과 `N-1` 중 답은 **N-1**이다 — N개 router에 대한 spanning tree의 edge 수.
receiver를 세는 모델이 더하는 N번째 link가 정확히 attach link이고, 그러면 같은 wire를 두 번
과금한다. 16×16에서 256 → **255**.

- 닫힌 형태 대신 **edge set을 열거**하는 `spatial_multicast_edge_count()`를 두고,
  `spatial_noc_cost()`가 그것과 일치해야 한다. Unit test가 1×1 / 1×N / N×1 / N×M을 대조한다
  (1×1은 내부 link가 **0**이다 — 데이터가 전부 attach link로 들어온다).
- NoP 계약은 **다르고 그것이 맞다**: NoP의 source는 multi-chip staging buffer이고 chip 0으로
  들어가는 link는 다른 축이 과금하지 않으므로 **포함한다**. `nop_delivery_cost()`가 이미
  명시적으로 그 ingress link를 더하고 있었다.
- 두 계약을 **report에 출력**한다(`NoC link contract`, `NoP link contract`).

모든 shipped config가 `noc_energy = 0`이므로 baseline은 움직이지 않는다.

검증 — [validation/noc/check.py](../../validation/noc/check.py) NC1(255) / **NC6**(report가 계약
문장을 담는지). Negative test로 N을 되돌리니 NC1이 `131072 != 130560`으로 실패했다.

## RE7 — label과 specification unit 표기 일관성

- Network rollup이 `Network dynamic/static/total energy`를 출력한다. 그전에는 `network.txt`가
  layer와 같은 label을 써서 network total이 한 layer의 것으로 읽혔다.
- Specification 출력의 하드코딩 `pJ` 15곳을 `energy_units().label()`로 바꿨다
  (`pe.cc` 6, `global_buffer.cc` 4, `adder_tree.cc` 2, `systolic_array.cc`/`spatial_arch.cc`/
  `multi_chip.cc` 각 1). normalized fixture가 자기 입력은 pJ로, 출력은 normalized로 설명하고
  있었다.
- Checker가 두 scope를 구분한다: `network.txt`에 `Layer ...`가 남아 있으면 실패, layer file에
  `Network ...`가 있으면 실패.

검증 — energy checker **E12**. Negative test로 label과 `pJ` 한 곳을 되돌리니 양쪽이 잡혔다.

## Regression

Config 267 → 275. 외부 reference는 전부 불변:

| reference | 결과 |
| --- | --- |
| Gemmini Verilator RTL | MAPE 4.40% / max 7.86% (limit 8%/8%) |
| Eyeriss silicon latency | MAPE 4.26% / max 6.39% (limit 5%/8%) |
| SCALE-Sim cross-simulator gate | PASSED (fc 6.0%/8%, steady-conv 15.0%/20%, rho 0.929) |
| Eyeriss DRAM traffic | MAPE 22.95% / max 50.0% (기존 gate 범위 외, 변화 없음) |
| Eyeriss GLB band | 1/5 layer inside (PA9 missing source, 변화 없음) |

Gate 17종 전부 통과. Unit test는 sanitizer(ASan/UBSan) 포함 통과.

부수적으로 checker 두 개가 이번 변경에 걸려서 같이 고쳤다: `validation/dram/check.py`가 unit을
`pJ`로 하드코딩하고 있었고, `validation/leakage/check.py`가 `network.txt`에서 `Layer ...` label을
읽고 있었다. 둘 다 config 선언을 따르도록/두 scope를 모두 받도록 바꿨다.

## 남긴 것

**Phase 4** — measured/calibrated 외부 energy 및 traffic reference가 필요하다. 지금은
`energy_reference = synthetic-fixture-not-measured-silicon`처럼 fixture가 자기 출처를 명시하고,
RE2의 자격 판정이 미교정 config에 absolute 수치를 주지 않는다. 실제 silicon reference가 확보되면
그 문자열만 바꾸면 absolute 검증이 열린다.

**PE accumulator와 edge accumulator를 한 fixture에서 동시에 non-zero로** 만들 수는 없다 —
`edge_accumulation`이 소유권을 배타적으로 가르기 때문이다. AC4/AC9가 A/B로 양쪽 소유권을
검증하므로 E11의 `require_nonzero`에서는 edge 쪽만 요구한다.
