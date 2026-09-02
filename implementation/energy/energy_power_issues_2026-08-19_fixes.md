# 2026-08-19 energy·power 재평가 수정 기록 (Phase 1/2/3/5)

기준 문서: [../../assessment/ENERGY_POWER_ISSUES_2026-08-19.md](../../assessment/ENERGY_POWER_ISSUES_2026-08-19.md)

평가 문서의 권장 순서(Phase 1 → 5)를 따랐고, **Phase 4(외부 calibration)만 제외**했다 — 최소 한
accelerator의 measured component energy 또는 calibrated estimator reference 확보가 선행되어야
하므로 코드 수정으로 닫을 수 없다.

## 우선순위 조정

권장 순서에 동의했다. 한 가지만 옮겼다: **E8의 finite/non-negative guard와 "0 vs uncalibrated"
구분을 Phase 3 → Phase 1로 당겼다.** Phase 1 항목 3(unit provenance)과 같은 observability
문제이고, 어차피 같은 시점에 config 파싱을 건드리기 때문이다.

각 issue를 해결한 직후 검증을 붙였고, 모든 항목에 **negative test로 검출력을 확인**했다.

## Phase 1 — Energy 결과 계약

### E1a. 수집되지만 출력되지 않던 energy

세 항목이 component와 `stats_t`에 누적되고 network에도 합산되면서 **결과 파일에 전혀 나오지
않았다**: `access_energy_pe_array`, `access_energy_multi_chip`, `transfer_energy_dram`.

특히 `gemmini_dram_detail.cfg`의 `row_miss_energy = 20`은 내부 값만 바꾸고 사용자에게 아무것도
보여주지 않았다. 이제 `Interconnection energy (link + row activation)`으로 출력된다.

검증:
- **DR7** ([validation/dram/check.py](../../validation/dram/check.py)) 수기 항등식
  `streams × rows × row_miss_energy`. input 4×4×20=320, output 1×16×20=320.
  energy는 bank 병렬성의 이득을 받지 않으므로 8-bank/1-bank 양쪽에서 동일해야 하고, knob off면 0.
- **E6** ([validation/energy/check.py](../../validation/energy/check.py)) 14개 필수 energy 라인의 존재.
- 검출력: knob을 40으로 바꾸면 `640 vs hand 320`으로 실패.

### E1b. Component subtotal과 layer/network total

먼저 **이중 계상이 없는 경계를 확인하고 명문화**했다. 모든 전송이 세 개의 서로 다른 물리 자원을
과금하며 각각 그 자원의 소유 component에 청구된다.

```text
source buffer의 read port  -> SOURCE component의 access energy
건너는 link/fabric         -> 그 link를 소유한 component
destination buffer의 write -> DESTINATION component의 access energy
```

따라서 subtotal 합은 근사가 아니라 정확하다. leakage는 event count가 아니라 layer wall-clock의
함수이므로 별도 축으로 유지한다.

검증: **E7** — 요약의 subtotal을 **같은 리포트의 상세 per-datatype 라인에서 다시 합산**해 대조
(요약 자기 숫자를 쓰지 않으므로 둘이 반드시 일치해야 한다). **E8** — network total = layer 합.

### E7. Energy 단위 provenance

`eyeriss_energy.cfg`는 정규화 비용(MAC=1, RF=1x, GIN=2x, GLB=6x, DRAM=200x)인데 출력은 전부 `pJ`였다.
`energy_unit = pJ|normalized`와 `energy_reference`를 추가하고 **49개 출력 지점 전부**가 선언된
단위를 따르게 했다.

```text
Energy unit : normalized (relative to the declared reference; ABSOLUTE totals are NOT meaningful)
              -- provenance: Chen-ISCA2016-TableIV-normalized-to-one-MAC
```

checker 자신도 "이 fixture는 정규화이므로 절대 pJ를 검증하지 않는다"를 명시한다.

구현: [utils/energy_units.h](../../utils/energy_units.h) / [.cc](../../utils/energy_units.cc).
검증: **E9**. 검출력: 선언을 주석 처리하면 3건 실패.

### E12. Energy scope

`Energy scope: N of M layers`를 total 바로 위에 출력한다. 미지원 layer의 energy는 **추정도 0 채움도
아니고 그냥 없다**. 선언한 layer 수가 실제 rollup된 파일 수와 일치하는지도 검사해 장식용 문장이
되지 않게 했다(**E10**).

### E8. Config guard

22개 energy 키 중 negative guard가 있던 것은 3개뿐이었다. `validate_energy_settings()`가 config의
**모든** energy 설정을 검사한다 — finite, non-negative, 숫자 형식, 값 누락. 컴포넌트마다 자기 키를
기억할 필요가 없다.

또한 "silently free component"를 보이게 했다:

```text
 * PE array   0.00   0.00   0.00   <- no energy cost declared
Uncosted components   : 2 of 6 contribute no energy
```

검증: unittest negative fixture 4종(negative / non-numeric / negative compute / non-finite) +
**positive fixture 1종**(0은 정당한 modeled cost이므로 허용되어야 한다). runtime 경로에서도
거부됨을 확인했다.

## Phase 2 — 현재 산식의 직접 결함

### E2. Mesh multicast physical-link energy

16×16 array에서 multicast multiplier가 **평균 Manhattan hop 15**였다. 256개 endpoint로 fanout하는
delivery는 256개 link를 쓴다 — **약 17배 과소계상**. 반면 output writeback은 per-source unicast이므로
평균 hop이 맞다. 두 방향을 분리했다(평가 문서 완료조건 5).

```text
MULTICAST (operand distribution) : link traversals = h*w   (PE마다 자기 incoming link)
UNICAST   (output write-back)    : link traversals = 평균 Manhattan 거리
```

검증: [validation/noc/check.py](../../validation/noc/check.py) NC1~NC5.
`noc_energy = 2`인 16×16 / 8×8 fixture, **transaction 수가 동일함을 확인**해 차이가 shape에만
귀속되게 했다.

| shape | multicast in/wt | unicast out |
|---|---|---|
| 16×16 | 256 × 2 × **256** = 131,072 | 16384 × 2 × **15** = 491,520 |
| 8×8 | 256 × 2 × **64** = 32,768 | 16384 × 2 × **7** = 229,376 |

검출력: pre-E2 형태에서 7,680 vs 131,072.

### E3. Precision-dependent MAC energy

MAC energy가 format과 무관한 scalar 하나였다 — INT4 곱셈이 FP16 곱셈과 같은 값이었고 아무도
이상으로 판단하지 않았다. `mac_energy_<input>_<weight>` table을 추가했고, **미선언 조합은 임의
scalar를 재사용하되 UNCALIBRATED로 명시**한다(거부는 237개 config을 전부 깨므로 평가 문서가 허용한
대안을 택했다).

```text
MAC energy basis : computation_energy (UNCALIBRATED for int8 x int8: the same scalar is used
                   for every operand precision)
```

검증: [validation/precision/check.py](../../validation/precision/check.py) PR1~PR5.
MAC 수 동일(262,144)에서 compute는 선언 단가를, memory는 **datatype별** 폭(1:2:4)을 따르고
output은 `output_format` 고정이므로 불변 — precision이 전역이 아니라 tensor별로 추적됨을 보인다.
compute 비율 3.00 vs traffic 비율 1.67로 **두 축의 독립성**도 고정했다.

### E4. Accumulator format 연결

`accumulator_format`이 파싱·출력되고 소비되는 곳이 없었다. spill이 **output** datatype으로 크기가
정해져 fp32와 fp16 accumulator가 동일한 traffic/energy를 냈다.

두 가지를 고쳤다.
1. spill은 accumulator 정밀도, 최종 cast/pack은 output 정밀도 — **서로 다른 event로 분리**.
2. link **transaction** 단위로는 정밀도를 표현할 수 없다(MAC tile이 몇 element라
   `ceil(elements×width/link_bits)`가 fp32·fp16 모두 1로 반올림된다 — 이것이 애초에 accumulator
   format이 무력했던 이유다). 이들은 link 횡단이 아니라 **값에 대한 정밀도 변환**이므로
   **byte 기준**으로 과금한다.

또한 reduction tree energy를 `computation_energy`에서 분리해 자체 축으로 출력한다.

검증: [validation/accumulator/check.py](../../validation/accumulator/check.py) AC1~AC5 —
fp32 spill = 4×cast, fp16 = 2×cast, cast는 불변, energy 항등식.
검출력: pre-E4 형태에서 ratio 1.0.

### E5. Fold/setup/drain dynamic energy

`weight_fold_fill_cycle` / `layer_setup_cycle`은 latency를 늘리면서 같은 활동의 dynamic energy는
0이었다(leakage 제외). latency가 정밀해질수록 energy 쪽 비대칭이 커지는 구조였다.

`weight_fold_fill_energy` / `layer_setup_energy`를 추가하고 **latency와 동일한 선택 규칙**을
적용했다 — calibrated per-element 비용이 우선, analytical drain은 tile 경계당, 둘을 동시에
과금하지 않는다. latency를 공급한 쪽이 event count도 공급하므로 이중 과금이 구조적으로 불가능하다.

검증: [validation/fold_energy/check.py](../../validation/fold_energy/check.py) FS1~FS5.
16 events × 3.0 = 48, 단가 2배 → 정확히 2배, **compute schedule 불변**, layer total이 정확히
548만 이동(containment).

## Phase 3 — Non-zero regression

### E9 (reduction). `mac_reduction_energy`

모든 shipped config가 `mac_width = 1`이라 adder tree에 작업이 없었다 — 이 knob은 **어떤 배수로
틀려도, 잘못된 lane 수에 연결되어도 아무 검사도 알아채지 못했다**. 4-lane MAC fixture를 만들었다.

```text
262,144 scalar MACs / 4 lanes = 65,536 issues × 3 additions × 0.4 = 78,643.2
```

검증: [validation/reduction/check.py](../../validation/reduction/check.py) RD1~RD4.
검출력: N-1 대신 N additions로 바꾸면 104,857.6로 실패.

### E12. Multi-chip NoP energy

**2-chip mesh로는 판별이 불가능함을 먼저 확인했다** — multicast tree와 unicast route가 같은 2개
link를 덮어 energy가 동일하다. **1×4 배치**로 fixture를 다시 만들었다(unicast 1+1+2+3=7 traversal,
multicast tree 4).

리포트만으로 닫히는 항등식:

```text
energy = reported txns × nop_energy × total_traversals / bottleneck_copies
```

핵심은 **energy 비율 1.75(=7/4)와 transaction 비율 4.00(=4/1)이 다르다**는 점이다 — latency와
energy가 실제로 다른 양을 쓴다는 증거이며, 하나를 둘에 쓰는 모델은 두 datatype을 동시에
만족시킬 수 없다.

검증: [validation/nop_energy/check.py](../../validation/nop_energy/check.py) NP1~NP4.

### E11. Leakage production path

기존에는 `unit × elapsed` helper 단위 테스트 하나뿐이었고, production scaling·final-latency
rescale·network rollup은 전부 미검증이었다. 모든 config이 `static_energy = 0`이라 경로 전체가
죽은 코드일 수도 있었다.

DRAM만 10배 늦춘 A/B(= 같은 energy, 더 긴 window):

| | critical path | dynamic E | static E |
|---|---:|---:|---:|
| baseline | 59,010 | 1,260,831.96 | 460,278.00 |
| slow DRAM | 384,942 | **1,260,831.96 (불변)** | 3,002,547.60 |

static 비율과 critical-path 비율이 **6.5233으로 정확히 일치** — leakage가 최종 wall-clock에
과금됨을 고정한다.

검증: [validation/leakage/check.py](../../validation/leakage/check.py) LK1~LK5.
검출력: rescale을 제거하면 5.547 vs 6.523.

## Phase 5 — Power 모델

### 확정 계약

작지만 완결된 계약을 택했다.

```text
time  = critical-path cycles / authoritative frequency
power = energy / time
EDP   = energy × time        ED2P = energy × time²
```

**세 가지 전제를 가정하지 않고 출력한다.**

1. **ONE CLOCK.** timeline은 하나의 공유 cycle 축이므로(L5), cycle을 초로 바꿀 수 있는 것은 모든
   modeled component가 같은 clock일 때뿐이다. mixed domain은 임의 domain으로 나눈 숫자를 내놓지 않고
   **unsupported**로 보고한다.
2. **ABSOLUTE ENERGY.** 정규화 fixture에는 watt가 없다(E7). `energy_unit = pJ`일 때만 power를
   출력하고, 아니면 이유를 적는다.
3. **AVERAGE ONLY.** peak power는 개별 event의 concurrency가 필요하고 이 모델은 stage 수준
   overlap만 해결한다. 최대값의 합으로 근사하지 않고 **명시적 unsupported**로 둔다.

`static_energy`는 config에서 pJ/**cycle**이며 이미 layer 최종 wall-clock에 적분되어 있으므로
(E11의 rescale), 같은 시간창으로 나누면 이중 계상 없이 leakage power가 된다.

미포함 범위도 숫자 옆에 적는다: DRAM background/refresh, I/O termination, clock network,
controller/DMA, PHY — **core datapath only**.

```text
============= Power summary =============
Clock domain          : all modeled components share one clock
Authoritative clock   :          200.0 MHz
Not included          : DRAM background/refresh and I/O termination, clock network, controller/DMA, PHY; core datapath only
Peak power            : unsupported (needs per-event concurrency; this model resolves stage-level overlap only)
Elapsed time          :       0.295050 ms  (59010.0 cycles / 200.0 MHz)
Average dynamic power :          4.273 mW
Average leakage power :          1.560 mW
Average total power   :          5.833 mW
EDP                   :     507.813494 pJ*s
ED2P                  :       0.149830 pJ*s^2
```

### 검증

[validation/power/check.py](../../validation/power/check.py) PW1~PW6.

가장 강한 것은 **PW4**다. latency만 바꾼 A/B에서 두 축이 서로 다르게, 그러나 일관되게 움직여야 한다.

| | window | dynamic power | leakage power |
|---|---:|---:|---:|
| baseline | 0.2950 ms | 4.273 mW | 1.560 mW |
| slow DRAM | 1.9247 ms | **0.655 mW** (÷6.52) | **1.560 mW** (불변) |

dynamic energy는 불변이므로 window가 6.52배 길어지면 power는 정확히 그만큼 내려간다. leakage
energy는 나누는 그 window와 함께 커지므로 power는 불변이다. **한 번의 config 변경에서 두 결과를
동시에 얻는 것**이 시간 기준이 두 축에 일관되게 적용됨을 보인다.

검출력: 시간 기준을 compute schedule로 바꾸면 PW1과 PW4가 모두 실패한다.

## 검증 결과

| 검증 | 결과 |
|---|---|
| Unit/config validation | 통과, 267개 config |
| Gemmini RTL timing | MAPE 4.40%, max 7.86% (baseline 불변) |
| Eyeriss silicon latency | MAPE 4.26%, max 6.39% (baseline 불변) |
| SCALE-Sim cross-simulator gate | 통과 |
| Energy E1~E10 | 통과 |
| NoC energy NC1~NC5 | 통과 (신규) |
| Precision PR1~PR5 | 통과 (신규) |
| Accumulator AC1~AC5 | 통과 (신규) |
| Fold/setup FS1~FS5 | 통과 (신규) |
| Reduction RD1~RD4 | 통과 (신규) |
| NoP energy NP1~NP4 | 통과 (신규) |
| Leakage LK1~LK5 | 통과 (신규) |
| **Power PW1~PW6** | **통과 (신규)** |
| Traffic T1~T10 / DRAM DR1~DR7 / Bottleneck / Multi-chip / stats-timeline | 통과 |

외부 reference 3개 모두 불변이다.

## 남은 항목

### 외부 데이터가 필요한 것 (Phase 4)
- **E6** external traffic 정확도: batch filter reuse, psum GLB spill/reload 모델링.
  energy는 `event count × unit cost`이므로 traffic 오차가 memory energy로 거의 선형 전파된다.
  현재 상태는 축별 regression ceiling으로 고정돼 있다(DRAM 22.95%/50.0% RLC dense-equivalent,
  GLB band 1/5 layer 내부).
- 최소 한 accelerator의 **component별 measured energy 또는 calibrated estimator reference**.
  확보 전에는 absolute pJ 정확도를 주장할 수 없고, 현재 출력이 그렇게 말한다.
- precision별 compute/memory energy calibration — 현재 per-precision 값은 문법과 회귀만 있고
  값 자체는 fixture용 합성값이다(문서화됨).

### 모델 확장이 필요한 것
- **peak power**: per-event concurrency 필요. per-tile/per-packet event timeline이 더 정밀해진
  뒤에 진행하는 것이 안전하다(현재 명시적 unsupported).
- **E11 background/system power**: DRAM background/refresh, I/O termination, clock network,
  controller/DMA, PHY. 현재 범위는 core datapath로 출력에 한정 명시했다.
- **E12 sparse energy**: value/index/pointer traffic, zero gating, encoder/decoder,
  metadata buffer. timing entry가 sparse를 명시적으로 거부하므로 dense 비용으로 잘못 계산하는
  오류는 없다.
- **E8의 required unit cost 선언**: 어떤 component가 어떤 energy 키를 요구하는지의 schema가
  필요하다. 현재는 `uncosted components` 표시와 finite/non-negative guard로 실질 신호를 준다.
- 여러 clock domain 지원: per-domain cycle 정규화가 선행되어야 한다(현재 mixed domain은 power
  unsupported).
