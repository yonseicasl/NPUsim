# 2026-08-20 Phase 4 — 외부 absolute energy calibration (CACTI 7)

기준: [../../assessment/TIMING_SIMULATION_ENERGY_REEVALUATION_2026-08-19.md](../../assessment/TIMING_SIMULATION_ENERGY_REEVALUATION_2026-08-19.md)
의 Phase 4. 지금까지 유일하게 남아 있던 항목이다.

## 왜 남아 있었나

validation/의 다른 모든 gate는 모델을 **자기 자신과** 비교한다 — 수기 항등식, A/B 분리, containment.
외부 reference 셋(SCALE-Sim, Gemmini RTL, Eyeriss silicon)은 **latency와 traffic**을 검증한다.
absolute **energy**를 repo 밖의 무언가와 비교한 것은 없었고, 그래서 RE2는 어떤 shipped config의
총합도 absolute picojoule로 부르지 않고 wattage를 주지 않았다.

## 무엇을 실제로 얻을 수 있었나

기억에 의존한 근사값을 golden reference로 쓰지 않는다는 원칙으로 먼저 환경을 조사했다.

| 후보 | 결과 |
| --- | --- |
| **CACTI 7** | ✅ GitHub(HewlettPackard/cacti, commit `1ffd8dfb`)에서 clone·빌드 성공. 실제·재현 가능 |
| **Accelergy** | ✅ 설치·동작(commit `6911d156`, v0.4). 다만 CACTI plug-in의 실체가 CACTI이므로 같은 숫자에 배관만 추가 |
| Accelergy table plug-in | ❌ **reference로 거부.** `intmac.csv`/`DRAM.csv`가 예시 두 줄(45nm 16-bit MAC 3 pJ, LPDDR4 64-bit read 100 pJ)이고 유도 근거가 없다. 템플릿이지 데이터셋이 아니다 |
| **Simba (MICRO'19)** | ❌ DOI가 **HTTP 403**, arXiv에 없음. 검증 가능한 수치 출처가 없다 |

> Simba를 기억으로 채우지 않았다. 출처를 추적할 수 없는 golden reference는 없는 것보다 나쁘다 —
> 검증처럼 보이기 때문이다. 무엇이 필요한지는
> [validation/phase4/PROVENANCE.md](../../validation/phase4/PROVENANCE.md)에 적었다.

작업 중 arXiv ID 하나를 잘못 짚었다(2003.05666은 COVID-19 논문이었다). 확인하고 넘어갔다.

## Reference

[validation/phase4/cacti_sram_reference.csv](../../validation/phase4/cacti_sram_reference.csv) —
21행 = 7개 geometry × 3개 노드(45/32/22 nm). geometry는 **config가 선언한 값**을 쓴다. 모델이
맞아떨어지게 만드는 숫자가 아니라 "이 config가 주장하는 크기의 SRAM이 실제로 얼마인가"에 답한다.

22 nm 기준:

| geometry | read pJ | write pJ | leak pJ/cycle |
| --- | --- | --- | --- |
| pe_local_buffer_input (512 B) | 0.1819 | 0.1540 | 0.9415 |
| pe_local_buffer_weight (64 B) | 0.0404 | 0.0713 | 0.2065 |
| glb_input_128KB | 12.0501 | 10.5294 | 194.03 |
| glb_output_64KB | 8.4307 | 9.3146 | 104.20 |
| multi_chip_256KB | 18.9589 | 16.0138 | 384.83 |

## 발견 — shipped 값은 자기 geometry에서 재현되지 않는다

**기술 노드를 가정하지 않아도** 성립하는 불일치다. gate가 10건을 찾는다:

| section | 위반 |
| --- | --- |
| `[separate] read_energy` | input/weight 모두 128 KB인데 3.06 vs 1.07 (**2.86배**) |
| `[separate] read_energy` | output이 더 작은데(64 vs 128 KB) 더 비싸다 (3.06 vs 1.07) |
| `[separate] write_energy` | 같은 크기에 3.0 vs 0.96 (**3.12배**) |
| `[systolic_array] lb_read_energy` | weight가 더 작은데(64 vs 512 B) 더 비싸다 (0.44 vs 0.18) |
| `[systolic_array] lb_read_energy` | weight/output 모두 64 B인데 0.44 vs 0.33 (1.33배) |
| `[systolic_array] lb_write_energy` | 같은 크기에 0.51 vs 0.33 (1.55배) |

동일 용량 배열이 3배 다른 접근 에너지를 가질 수 없고, 더 작은 배열이 더 비쌀 수 없다. 즉 이 숫자들은
배열 geometry가 아닌 다른 것을 담고 있다 — per-datatype activity factor이거나 다른 단위 규약일
것이다.

**그 숫자들을 고치지 않았다.** 이 repo에 없는 post-layout Gemmini 데이터에서 왔을 수 있다. gate가
하는 일은 불일치를 **고정해서 계속 보이게** 만드는 것이다 — 어느 쪽을 조용히 맞춰서 없앨 수 없도록.

한 가지 유혹은 피했다: `0.18`이 CACTI 22nm의 512 B read `0.1819`와 1% 안에서 맞는다. 하지만 45nm의
64 B read도 정확히 `0.1800`이다. 한 숫자의 일치는 노드 추론의 근거가 되지 못하므로, 노드 추정
대신 **노드와 무관한 관계 위반**을 주장했다.

## 산출물

**`gemmini_cacti22.cfg`** — on-chip SRAM energy 전부가 CACTI에서 생성된, **repo 최초로 추적 가능한
출처를 가진 config**. `energy_unit = pJ`와 tool·commit·노드·cell type·온도를 담은
`energy_reference`를 선언한다. 그래서 RE2의 자격 판정을 통과하고 **정당하게 절대 전력을 얻는 첫
config**다: `pJ (absolute; totals and power are meaningful)`, **402.353 mW**.

DRAM access energy와 MAC compute energy는 gemmini.cfg에서 그대로 물려받았고 **CACTI 유래가 아니다**
(CACTI가 둘 다 모델하지 않고, Accelergy table plug-in은 템플릿뿐). `energy_reference`가 그렇게
명시한다.

**`validation/phase4/check.py`** P4-1~P4-5:

| 검사 | 내용 |
| --- | --- |
| P4-1 | fixture의 SRAM 값이 CACTI reference 행과 **정확히** 일치 — calibration이 출처에서 표류할 수 없다 |
| P4-2 | fixture의 geometry 일관성(같은 크기 → 같은 비용, 작은 쪽이 더 비싸지 않음) |
| P4-3 | **shipped gemmini.cfg의 위반이 여전히 존재**한다. 사라지면 실패한다 — 위반이 바람직해서가 아니라, 사라졌다면 config가 이 기록 없이 바뀌었거나 checker가 보기를 멈췄다는 뜻이므로 |
| P4-4 | RE2 end-to-end — calibrated fixture는 wattage를 받고 shipped config는 거부되며 ESTIMATED로 표시된다. Phase 4가 닫으려던 고리 |
| P4-5 | reference 자체의 sanity(노드가 미세해지면 에너지 감소, 배열이 커지면 증가) |

gate는 체크인된 CSV만 읽으므로 CI에서 CACTI가 필요 없다.

## 작업 중 잡은 제 실수

P4-2/P4-3의 첫 판은 `[systolic_array]`의 용량을 `input_buffer`로 찾았는데 gemmini는 `input_size`를
쓴다 — PE local buffer arm이 **조용히 건너뛰어졌고** docstring은 그 위반을 주장하고 있었다. 두 key
집합을 모두 시도하되 어느 것도 풀리지 않으면 "UNCHECKED"로 실패하게 고쳤다. 그 뒤 위반이 4건에서
10건으로 늘었다.

## 남은 것

**NoP link energy**(RE6의 계약)에는 여전히 외부 reference가 없다. Simba가 적임이고 유료다.
`PROVENANCE.md`가 닫으려면 필요한 것을 적었다 — 링크 통과당 또는 비트당 NoP 에너지와 그 노드, 같은
노드의 PE/GLB 접근 에너지(규모 비교용), 각 수치의 chip 수·array geometry. 그것이 있으면
`validation/phase4/check.py`에 SRAM arm과 같은 형태의 NoP arm이 붙는다.

---

# Phase 4 DRAM 축 — DRAMsim3 (2026-08-20)

CACTI는 DRAM을 모델하지 않아 off-chip 항이 미교정으로 남아 있었다. **DRAMsim3가 이미
`ext/DRAMsim3`에 vendored** 되어 있고 JEDEC device 파일(datasheet IDD 전류와 VDD)로 명령당 에너지를
계산한다 — 즉 NPUsim config가 `dram_config`에 **이름을 적어둔 바로 그 device**에 대한 정당한
reference다.

## 단위 확정

DRAMsim3의 에너지 누산 단위를 두 방법으로 교차 확인했다:

- `read_energy / num_read_cmds` = 2894.4 이고 소스의 공식
  `VDD*(IDD4R-IDD3N)*burst_cycle*devices` = 1.35 × (185−51) × 4 × 4 = 2894.4 로 **정확히 일치**
- 자체 보고 `average_power` = `total_energy / num_cycles` = 1310.84 로 **일치**

따라서 단위는 **V·mA·cycles = mW·cycle**이고, tCK(ns)를 곱하면 pJ가 된다.

## Reference

[validation/phase4/dramsim3_dram_reference.csv](../../validation/phase4/dramsim3_dram_reference.csv) —
NPUsim config가 이름을 적는 6개 device 전부.

| device | tCK | burst | read pJ/burst | read pJ/byte | activate pJ |
| --- | --- | --- | --- | --- | --- |
| DDR3_8Gb_x16_1600 | 1.25 | 64 B | 3618.0 | **56.531** | 5325.8 |
| DDR3_8Gb_x16_1866 | 1.07 | 64 B | 3281.9 | 51.280 | 5286.9 |
| DDR4_8Gb_x16_2400 | 0.83 | 64 B | 3123.5 | 48.804 | 8242.9 |
| DDR4_8Gb_x16_3200 | 0.63 | 64 B | 2975.6 | 46.494 | 9991.3 |
| LPDDR4_8Gb_x16_2400 | 0.83 | 128 B | 6246.9 | 48.804 | 6956.1 |
| HBM2_4Gb_x128 | 1.0 | 64 B | 804.0 | **12.562** | 828.0 |

P4-6이 이 표를 device 파일의 파라미터에서 **독립 재유도**한다(DDR3/DDR4 4종). 즉 "도구가 출력한
값"을 그대로 믿는 것이 아니다 — 행 하나를 훼손하면 검출된다.

## 발견 — 선언값이 device와 무관하다

| config | 이름 적은 device | 선언 pJ/byte | DRAMsim3 pJ/byte | 배수 |
| --- | --- | --- | --- | --- |
| eyeriss | DDR3-1600 | 8.000 | 56.531 | 7.07× |
| eyerissv2 | DDR4-3200 | 8.000 | 46.494 | 5.81× |
| fsd | LPDDR4-2400 | 8.000 | 48.804 | 6.10× |
| tpu | DDR3-1866 | 8.000 | 51.280 | 6.41× |
| tpuv3 | **HBM2** | 8.000 | **12.562** | 1.57× |
| maeri | DDR3-1600 (bw 256) | 2.000 | 56.531 | 28.27× |

**6개 device 전부에 8.0 pJ/byte가 선언돼 있다.** 그 device들은 per-byte 에너지가 4.5배 차이 난다
(HBM2 12.56 vs DDR3-1600 56.53). 즉 config가 HBM2를 적어도 DDR3급 에너지를 과금받고, 그 역도 같다.

## 그런데 timing은 device에서 유도된다 — 이것이 발견을 날카롭게 만든다

`gemmini_dram_detail`의 `t_ras_cycle = 7`은 DDR3-1600의 **28 tCK × 1.25 ns / 5 ns = 7.00**과 정확히
일치한다(`t_rp_cycle = 3` ≈ 2.75도). 즉 **같은 fixture에서 DRAM timing은 datasheet를 cycle 단위로
따르는데 energy는 전혀 따르지 않는다.** "아무도 datasheet를 안 봤다"가 아니라 "timing만 거기서
왔다"가 정확한 진술이다. P4-9가 이 대조를 고정한다.

`dram_config` 자체는 어느 쪽도 해소하지 못한다: `#ifdef DRAMSIM3` 뒤에서만 읽히고 `npusim.sh`는
`DRAMSIM3=0`이 기본이므로, 검증되는 모든 빌드에서 device 이름은 **아무 효과가 없다**. device를
HBM2로 바꿔도 보고 필드가 0개 움직인다(측정 확인).

> 이 때문에 knob gate의 `dram_config` 면제 사유("external DRAM model selection")는 기본 빌드에서는
> 너무 관대한 표현이다. P4-9가 그 사실을 명시적으로 기록한다.

## 산출물

`gemmini_cacti22.cfg`의 DRAM 항을 이름 적힌 device에서 채웠다:

| key | 값 | 출처 |
| --- | --- | --- |
| `read/write_energy` | 452.25 pJ/txn | 56.531 pJ/byte × 8 B |
| `row_miss_energy` | 5325.75 pJ | JEDEC activate |
| `row_buffer_size` | 8 KB | 1024 columns × 16 b × 4 devices |
| `t_ras_cycle` / `t_rp_cycle` | 7 / 2.75 | 28·11 tCK → 200 MHz |

이제 이 fixture는 **on-chip SRAM(CACTI) + off-chip DRAM(DRAMsim3)** 이 모두 외부 유래이고,
`418.685 mW`를 보고하며 `DRAM timing model: analytical open-page, JEDEC tRC=tRAS+tRP`를 쓴다.
MAC compute energy만 미교정으로 남는다.

검출력:

| 주입 | 결과 |
| --- | --- |
| 교정 fixture의 DRAM energy를 64로 되돌림 | **P4-7 실패** (`64.0 vs 452.25`) |
| reference 행 하나 훼손(3618 → 3000) | **P4-6 실패** — 재유도 값과 불일치, 그리고 재유도 device 수가 4→3으로 줄어 두 번째 실패도 발생 |

## Phase 4 현황

| 축 | 상태 |
| --- | --- |
| on-chip SRAM | ✅ CACTI 7 |
| off-chip DRAM device | ✅ DRAMsim3 (JEDEC IDD) |
| absolute power 자격(RE2 고리) | ✅ 닫힘 |
| off-chip link (`[dram] transfer_energy`) | ❌ DRAMsim3는 device만 모델 |
| MAC compute energy | ❌ 합성 또는 노드 명시된 공개 수치 필요 |
| NoP link energy (RE6 계약) | ❌ Simba 유료 |
