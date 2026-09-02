# NPUsim 하드웨어 정합성 평가 — Eyeriss silicon, Gemmini RTL, NVDLA RTL

**평가일:** 2026-08-25  
**평가 metric:** cycle, energy  
**평가 대상:** 현재 NPUsim working tree와 저장소에 보존된 외부 golden/provenance

## 1. 결론 요약

현재 NPUsim은 **Eyeriss silicon과 Gemmini RTL에 대한 cycle 정합성은 평가된 workload 범위에서
높다**고 판단할 수 있다. 그러나 두 결과 모두 calibration에 사용한 workload와 최종 비교
workload가 완전히 분리된 독립 holdout 검증은 아니다.

**절대 hardware energy 정합성은 Eyeriss, Gemmini, NVDLA 모두 아직 입증되지 않았다.** Eyeriss는
공개된 normalized energy 비율에 대한 상대 분해 검증이 있고, Gemmini는 CACTI/DRAMsim3 기반의
일부 component cost reference가 있지만, 어느 경우도 동일 hardware 실행의 전체 energy golden과
NPUsim total energy를 직접 비교하지 않는다. NVDLA energy는 아직 미착수 상태다.

| 기준 | Cycle fidelity | Relative/internal energy | Absolute hardware energy |
| --- | --- | --- | --- |
| Eyeriss silicon | **높음 — calibrated, AlexNet conv 5개 한정** | **중간~높음** | **미입증** |
| Gemmini RTL | **높음 — calibrated, WS GEMM 6개 한정** | **내부 회계는 높음** | **미입증** |
| NVDLA RTL | **예비 n=1, 증거 부족** | **미검증** | **미검증** |

현재 가장 정확한 종합 주장은 다음과 같다.

> NPUsim은 평가된 Eyeriss silicon convolution과 Gemmini cycle-exact RTL workload에서 각각
> 4.26%와 4.40%의 cycle MAPE를 보인다. Energy 모델은 event accounting, normalized breakdown
> 및 일부 SRAM/DRAM component cost에 대해 검증되었지만, Eyeriss·Gemmini·NVDLA 어느 것에
> 대해서도 전체 절대 hardware energy fidelity는 아직 입증되지 않았다. NVDLA cycle 검증은
> CRC가 통과한 단일 micro-workload의 +14.2% 예비 결과에 머물러 있다.

## 2. 평가 방법과 코드 상태

### 2.1 현재 NPUsim 평가 상태

평가 마지막 재실행 시점의 repository HEAD는 다음과 같다.

```text
d3ef86e038bf872c0911f260fa6855c8a556c099
```

현재 working tree에는 사용자 작업으로 보이는 staged 변경이 존재했다.

- `components/adder_tree.cc`
- `components/pe.cc`
- `components/spatial_arch.cc`

`models/model`은 이 source 변경보다 나중에 생성된 binary였고 `./npusim.sh build npusim`은 추가
재빌드가 필요하지 않다고 판정했다. 공식 timing 및 energy validation을 이 상태에서 다시
실행했다. 이 문서는 commit 하나만의 결과가 아니라 **위 HEAD와 당시 working-tree 변경이 반영된
binary의 결과**를 평가한다.

### 2.2 재실행한 validation

```bash
./unittest/run_timing_validation.sh
python3 validation/energy/check.py
python3 validation/phase4/check.py
python3 validation/power/check.py
python3 validation/unpriced/check.py
```

결과:

- Gemmini RTL cycle gate: PASS
- Eyeriss silicon cycle gate: PASS
- SCALE-Sim cross-simulator gate: PASS
- Normalized energy accounting gate: PASS
- CACTI/DRAMsim3 Phase 4 component-reference gate: PASS
- Power qualification/refusal gate: PASS
- Unpriced-event gate: PASS

이 PASS 결과의 의미를 구분해야 한다. Cycle gate는 외부 RTL/silicon golden과 실제 수치를
비교한다. Energy 관련 gate는 대부분 component unit cost의 provenance, event×cost 산술,
containment, total identity 및 잘못된 absolute-power 주장 방지를 검사한다. **세 accelerator의
hardware total energy golden과 직접 비교하는 gate는 아니다.**

### 2.3 NVDLA 결과의 branch 범위

NVDLA fidelity 산출물은 현재 `xPU` working tree에 병합되어 있지 않고 별도 worktree에 있다.

```text
.claude/worktrees/nvdla-fidelity
HEAD b8c51fe7ffbd3a8072a242629554e5cd19eec820
```

현재 main working tree에는 `configs/accelerators/nvdla_small.cfg`와
`validation/nvdla/`가 없다. 따라서 NVDLA 결과는 현재 NPUsim HEAD에 대한 최신 회귀 결과가 아니라,
별도 worktree에서 확보된 **첫 예비 데이터 포인트**로만 평가한다.

## 3. Eyeriss silicon

## 3.1 Cycle 정합성

### 결과

동일 조건:

- Eyeriss 12×14 spatial array
- AlexNet batch 4
- 16-bit
- 200 MHz core
- silicon processing latency: Chen et al., JSSC 2017 Table V
- NPUsim 비교량: `Computation cycle + Fold fill cycle`

| Layer | NPUsim schedule | Silicon processing | Signed error |
| --- | ---: | ---: | ---: |
| conv1 | 15.47 ms | 16.5 ms | -6.2% |
| conv2 | 36.83 ms | 39.2 ms | -6.1% |
| conv3 | 21.28 ms | 21.8 ms | -2.4% |
| conv4 | 15.96 ms | 16.0 ms | -0.3% |
| conv5 | 10.64 ms | 10.0 ms | +6.4% |

```text
MAPE             = 4.26%
Maximum |error|  = 6.39%
n                = 5 convolution layers
```

현재 acceptance limit은 MAPE 5%, max 8%이며 모두 통과한다.

### 평가

**평가된 다섯 AlexNet convolution에서의 calibrated cycle fidelity는 높다.** 실제 silicon
processing latency가 기준이라는 점에서 증거의 질도 높다. 특히 layer 크기와 group 조건이 다른
conv1~conv5 전체가 6.4% 이내에 있다.

다만 다음 제약 때문에 범위를 넓혀 주장해서는 안 된다.

1. `weight_fold_fill_cycle=0.11`이 동일한 Eyeriss silicon 지점을 이용해 보정되었다.
2. calibration과 독립적인 holdout network가 없다.
3. 비교 대상은 convolution processing latency이며 pooling, activation 등 AlexNet의 나머지
   layer는 timing total에서 제외된다.
4. 다른 mapping, precision, batch, Eyeriss v2로 일반화되지 않는다.

### 허용 가능한 주장

> NPUsim은 평가된 Eyeriss AlexNet convolution 5개에서 silicon processing latency를 MAPE
> 4.26%, 최대 오차 6.39%로 재현한다.

### 허용되지 않는 주장

- 모든 Eyeriss workload에서 6.4% 이하의 cycle 오차
- AlexNet end-to-end latency 전체의 4.26% 정합성
- 독립 holdout에서 검증된 predictive accuracy

근거:

- [validation/phase3/README.md](../validation/phase3/README.md)
- [validation/phase3/golden_eyeriss_silicon.csv](../validation/phase3/golden_eyeriss_silicon.csv)
- [validation/check_timing.py](../validation/check_timing.py)

## 3.2 Energy 정합성

### 현재 검증 결과

`eyeriss_energy.cfg`는 ISCA 2016 Table IV의 normalized cost를 사용한다.

```text
energy_unit      = normalized
energy_reference = Chen-ISCA2016-TableIV-normalized-to-one-MAC
```

현재 코드로 재실행한 주요 결과:

| Layer | RF access/MAC | RF:rest | DRAM energy share |
| --- | ---: | ---: | ---: |
| conv1 | 5.21 | 3.40:1 | 11.3% |
| conv2 | 5.45 | 3.42:1 | 4.3% |
| conv3 | 5.23 | 3.17:1 | 6.3% |
| conv4 | 5.28 | 3.41:1 | 7.6% |
| conv5 | 5.28 | 3.40:1 | 8.6% |
| fc1 | 7.98 | 0.41:1 | 64.8% |
| fc2 | 7.97 | 0.41:1 | 64.8% |
| fc3 | 7.97 | 0.40:1 | 64.9% |

다음 검사가 모두 통과했다.

- ALU energy = MAC count × normalized MAC cost
- DRAM energy/access cycle 단가비
- GIN energy = serialized transaction × unit cost
- component subtotal과 layer total identity
- network energy = 실행된 layer energy 합
- optional energy axes의 containment
- normalized unit/provenance 표시
- partial network energy scope 표시

### 평가

상대적인 component breakdown과 event-cost coupling은 양호하다. RF-dominated convolution과
DRAM-dominated FC라는 Eyeriss의 공개된 energy 경향도 재현한다. 따라서 다음 용도로는 사용할 수
있다.

- 동일 cost basis에서 mapping 간 상대 energy 비교
- component별 energy 비중 비교
- traffic 변화에 따른 energy 증감 분석
- cost knob sensitivity 분석

그러나 이것은 절대 pJ 또는 joule 검증이 아니다.

- layer별 silicon energy golden(J)이 없다.
- clock network와 controller energy가 직접 보정되지 않았다.
- zero gating과 입력 데이터 switching을 RTL/silicon activity와 비교하지 않았다.
- normalized 단가에 시간 단위를 결합해 wattage를 출력하는 것은 validation에서 명시적으로
  거부된다.

따라서 Eyeriss energy의 최종 판정은 다음과 같다.

| 구분 | 판정 |
| --- | --- |
| Internal accounting | 높음 |
| Relative component breakdown | 중간~높음 |
| Workload 간 상대 energy | 제한적으로 사용 가능 |
| Silicon absolute energy | 미입증 |
| Silicon power | 미입증 |

근거:

- [validation/energy/README.md](../validation/energy/README.md)
- [validation/energy/check.py](../validation/energy/check.py)
- [configs/accelerators/eyeriss_energy.cfg](../configs/accelerators/eyeriss_energy.cfg)

## 4. Gemmini RTL

## 4.1 Cycle 정합성

### 결과

동일 조건:

- Chipyard `GemminiRocketConfig`
- Verilator cycle-exact RTL
- 16×16 weight-stationary array
- INT8 operand / INT32 accumulation
- 비교량: `Computation cycle + Fold fill cycle`

| GEMM M×K×N | NPUsim cycle | RTL cycle | Signed error |
| --- | ---: | ---: | ---: |
| 64×64×64 | 3,518 | 3,810 | -7.7% |
| 128×128×128 | 11,358 | 10,530 | +7.9% |
| 256×256×256 | 71,390 | 69,653 | +2.5% |
| 16×512×512 | 32,990 | 32,506 | +1.5% |
| 512×512×64 | 69,598 | 74,413 | -6.5% |
| 512×512×512 | 540,894 | 543,240 | -0.4% |

```text
MAPE             = 4.40%
Maximum |error|  = 7.86%
n                = 6 GEMM workloads
```

현재 acceptance limit은 MAPE 8%, max 8%이며 통과한다.

### 평가

**평가된 여섯 WS GEMM에서의 calibrated cycle fidelity는 높다.** 작은 문제에서 setup/fold overhead가
지배하는 구간부터 512³ steady-state까지 포함하며 모든 점이 8% 이내다.

그러나 다음 두 parameter가 동일한 6개 RTL 지점에서 보정되었다.

- `weight_fold_fill_cycle=14`
- `layer_setup_cycle=2270`

따라서 이 결과는 높은 calibrated agreement이지만 독립 holdout accuracy는 아니다. 또한 다음은
검증 범위 밖이다.

- Gemmini convolution
- 다른 array dimension
- output-stationary dataflow
- sparse execution
- 다른 scratchpad/accumulator geometry
- DMA bandwidth/latency 변화에 대한 RTL sensitivity

### 허용 가능한 주장

> NPUsim은 평가된 6개 Gemmini 16×16 WS GEMM에서 cycle-exact Verilator RTL 대비 MAPE 4.40%,
> 최대 오차 7.86%를 보인다.

근거:

- [validation/phase2/README.md](../validation/phase2/README.md)
- [validation/phase2/golden_rtl_cycles.csv](../validation/phase2/golden_rtl_cycles.csv)
- [validation/phase2/compare.py](../validation/phase2/compare.py)

## 4.2 Energy 정합성

### Gemmini RTL energy golden의 부재

Gemmini Verilator validation은 cycle을 측정하지만 energy는 측정하지 않는다.

- RTL VCD/SAIF 기반 switching activity energy 없음
- 합성 netlist의 hierarchy별 power report 없음
- RTL 실행의 total dynamic/leakage energy golden 없음
- 동일 여섯 GEMM에 대한 component energy golden 없음

따라서 cycle-exact RTL이 존재한다는 사실을 energy fidelity의 근거로 확장할 수 없다.

### Shipped `gemmini.cfg`

기본 `gemmini.cfg`는 다음 이유로 absolute energy 자격을 얻지 못한다.

- `energy_unit`과 추적 가능한 energy provenance가 없음
- 동일한 128 KiB input/weight GLB가 3.06/1.07처럼 2.86배 다른 read cost를 가짐
- 더 작은 64 B weight local buffer가 512 B input buffer보다 더 비싼 cost를 가짐
- 동일 DRAM device에 대한 declared energy가 DRAMsim3 reference보다 크게 다름

현재 power validation은 이 config가 wattage를 출력하면 실패하도록 되어 있다. 실제 재실행에서도
wattage가 거부되고 total은 `ESTIMATED`로 표시되었다.

### CACTI/DRAMsim3 component reference

`gemmini_cacti22.cfg`에는 다음 외부 reference가 있다.

- On-chip SRAM: CACTI 7, 22 nm, declared geometry
- DRAM device: DRAMsim3 JEDEC IDD model

Phase 4 checker는 config의 SRAM/DRAM 단가가 reference CSV와 정확히 일치하고 geometry 관계가
타당함을 검증한다. 이 부분은 component-level external calibration으로 인정할 수 있다.

그러나 전체 Gemmini energy는 아니다. 현재 실행에서 다음 6개 active event가 unpriced로 검출되어
wattage가 거부된다.

1. layer setup
2. weight fold
3. accumulator reload/spill
4. final output cast
5. Format-IP payload
6. PE local-buffer leakage

MAC compute energy도 Gemmini RTL 합성값 또는 post-layout reference로 보정되지 않았다. Off-chip link와
controller/clock network도 완전한 reference가 없다.

### Synthetic power fixture의 의미

`gemmini_power_calibrated.cfg`는 모든 active event에 cost를 선언하므로 64³ fixture에서 9.372 mW를
출력한다. 그러나 provenance는 다음처럼 명시되어 있다.

```text
energy_reference = synthetic-fixture-not-measured-silicon
```

이 결과는 energy/time 산술과 power qualification 로직을 검증하기 위한 synthetic test다. Gemmini
RTL 또는 silicon의 9.372 mW를 의미하지 않는다.

### 최종 판정

| 구분 | 판정 |
| --- | --- |
| Event×cost accounting | 높음 |
| SRAM component cost | 외부 CACTI reference 존재 |
| DRAM device cost | 외부 DRAMsim3 reference 존재 |
| MAC/control/clock/fold/setup 전체 | 미교정 또는 불완전 |
| Gemmini RTL total energy | 미입증 |
| Gemmini silicon energy/power | 미입증 |

근거:

- [validation/phase4/PROVENANCE.md](../validation/phase4/PROVENANCE.md)
- [validation/phase4/check.py](../validation/phase4/check.py)
- [configs/accelerators/gemmini_cacti22.cfg](../configs/accelerators/gemmini_cacti22.cfg)
- [configs/accelerators/gemmini_power_calibrated.cfg](../configs/accelerators/gemmini_power_calibrated.cfg)
- [validation/power/check.py](../validation/power/check.py)
- [validation/unpriced/check.py](../validation/unpriced/check.py)

## 5. NVDLA RTL

## 5.1 Cycle 정합성

### 확보된 결과

별도 `worktree-nvdla-fidelity`에서 다음 조건의 NVDLA RTL 실행이 완료되었다.

- `nvdla/hw` master commit `1a65f1f5b48268accaa47c95f95c2601918be095`
- `nv_small`
- Verilator 5.022
- official trace `dc_1x1x8_1x1x8x1_int8_0`
- INT8 direct convolution
- NVIDIA embedded golden CRC PASS (`0x8f68a2ae`)

| Metric | RTL | NPUsim | Error |
| --- | ---: | ---: | ---: |
| End-to-end cycle | 169 | 193 | +14.2% |

비교 스크립트 재실행 결과:

```text
MAPE             = 14.2%
Maximum |error|  = 14.2%
n                = 1
```

이 값은 published architecture parameter로 만든 first-principles config에서 얻었으며 169 cycle에
맞춰 역산한 결과는 아니다. 따라서 first data point 자체는 의미가 있다.

### 왜 fidelity 주장으로는 부족한가

1. 표본이 n=1이다.
2. calibration/holdout 분리가 없다.
3. 현재 최신 NPUsim HEAD가 아니라 별도 과거-base worktree의 결과다.
4. `core_cycles`가 아니라 output DBB write를 completion proxy로 한 end-to-end cycle이다.
5. C=128, C=192처럼 multi-Atomic-C-group이 필요한 official trace는 all-zero output과 CRC mismatch가
   발생했다.
6. CRC mismatch의 원인이 해결되지 않아 큰 convolution의 cycle을 golden으로 신뢰할 수 없다.
7. 현재 main branch에는 NVDLA config, mapping, validation gate가 병합되어 있지 않다.

따라서 현재 NVDLA cycle 판정은 다음과 같다.

| 항목 | 판정 |
| --- | --- |
| RTL build/trace bring-up | 성공 |
| Single micro-layer functional golden | 성공 |
| First-principles cycle 비교 | +14.2%, n=1 |
| Calibration suite | 없음 |
| Independent holdout | 없음 |
| Current-head regression | 없음 |
| NVDLA cycle fidelity 주장 | 아직 불가 |

근거:

- [NVDLA first result](../.claude/worktrees/nvdla-fidelity/assessment/NVDLA_RTL_FIDELITY_FIRST_RESULT_2026-08-20.md)
- [NVDLA validation README](../.claude/worktrees/nvdla-fidelity/validation/nvdla/README.md)
- [NVDLA cycle golden](../.claude/worktrees/nvdla-fidelity/validation/nvdla/golden_rtl_cycle.csv)
- [NVDLA provenance](../.claude/worktrees/nvdla-fidelity/validation/nvdla/PROVENANCE.md)

## 5.2 Energy 정합성

NVDLA energy validation은 아직 시작되지 않았다.

- `nvdla_small.cfg`의 모든 energy cost가 0
- `golden_rtl_energy.csv` 없음
- VCD/FST/SAIF activity file 없음
- 합성 netlist power report 없음
- NVDLA SRAM macro power model 없음
- logic/clock/CBUF/CACC/DRAM component breakdown 없음
- relative energy comparison script 없음

따라서 현재 판정은 다음과 같다.

| 구분 | 판정 |
| --- | --- |
| Internal NVDLA energy accounting | 미검증 |
| RTL relative switching energy | 미검증 |
| Synthesized implementation energy | 미검증 |
| NVDLA absolute energy/power | 미검증 |

Plain Verilator cycle 실행은 energy golden을 제공하지 않는다. 절대 energy를 평가하려면 최소한
동일 measurement window의 activity, 합성 library/netlist power, CBUF/CACC SRAM characterization,
DBBIF traffic 기반 DRAM energy가 필요하다.

근거:

- [assessment/NVDLA_RTL_FIDELITY_PLAN_2026-08-20.md](NVDLA_RTL_FIDELITY_PLAN_2026-08-20.md)

## 6. 세 기준의 직접 비교

### 6.1 Cycle evidence 수준

| 기준 | External reference | 표본 | 현재 오차 | Calibration 독립성 | 평가 |
| --- | --- | ---: | ---: | --- | --- |
| Eyeriss | measured silicon | 5 conv | MAPE 4.26%, max 6.39% | 동일 지점 보정 | 높음, 범위 제한 |
| Gemmini | cycle-exact RTL | 6 GEMM | MAPE 4.40%, max 7.86% | 동일 지점 보정 | 높음, 범위 제한 |
| NVDLA | cycle-exact RTL | 1 micro-conv | +14.2% | 분리 불가 | 예비 단계 |

Eyeriss는 silicon reference라는 장점이 있고 Gemmini는 정확한 RTL cycle counter라는 장점이 있다.
두 경우 모두 표본 범위에서는 강한 agreement를 보이지만, independent holdout이 없으므로 새로운
shape/architecture parameter에 대한 predictive fidelity는 아직 별도 검증이 필요하다.

### 6.2 Energy evidence 수준

| 기준 | 외부 energy reference | 비교 방식 | Absolute total 평가 |
| --- | --- | --- | --- |
| Eyeriss | normalized published costs와 qualitative breakdown | relative component ratio | 불가 |
| Gemmini | CACTI SRAM + DRAMsim3 device 일부 | component unit-cost provenance | 불가 |
| NVDLA | 없음 | 없음 | 불가 |

세 architecture 중 hardware의 동일 layer 실행 total joule과 NPUsim total joule을 직접 비교한 사례는
현재 **0개**다.

## 7. 주장 가능 범위

### 현재 주장해도 되는 것

1. Eyeriss AlexNet convolution 5개에서 silicon processing latency MAPE 4.26%, max 6.39%.
2. Gemmini 16×16 WS GEMM 6개에서 cycle-exact RTL MAPE 4.40%, max 7.86%.
3. Eyeriss normalized energy breakdown에서 RF-dominated convolution, DRAM-dominated FC 경향 재현.
4. NPUsim energy의 component 합, network rollup, unit propagation, event-cost coupling 및 missing-cost
   detection은 강하게 검증됨.
5. CACTI/DRAMsim3에서 얻은 SRAM/DRAM cost를 NPUsim에 재현 가능하게 연결할 수 있음.
6. NVDLA nv_small RTL build와 CRC-passing micro-convolution 한 건이 확보되었고 first-principles
   NPUsim cycle 오차는 +14.2%였음.

### 현재 주장하면 안 되는 것

1. NPUsim의 모든 Eyeriss/Gemmini workload cycle 오차가 8% 이하라는 주장.
2. Eyeriss 또는 Gemmini의 독립 holdout predictive fidelity가 입증됐다는 주장.
3. Eyeriss, Gemmini, NVDLA의 절대 energy 또는 power fidelity가 높다는 주장.
4. `gemmini_power_calibrated`의 9.372 mW가 실제 Gemmini power라는 주장.
5. `gemmini_cacti22`가 complete Gemmini energy model이라는 주장.
6. NVDLA cycle fidelity가 검증 완료됐다는 주장.
7. NVDLA의 multi-Atomic-C convolution 결과를 golden으로 사용하는 것.

## 8. 권장 후속 순서

### P0 — NVDLA cycle 검증 정상화

1. NVDLA worktree를 현재 HEAD에 rebase/merge
2. multi-Atomic-C-group CRC mismatch 해결
3. current-head에서 official trace 반복 재현
4. C-tail, K-tail, stride, padding, CBUF-fit workload 수집
5. calibration/holdout 분리
6. 최소 8개 이상의 holdout에서 cycle gate 평가

### P1 — 독립 cycle holdout 추가

1. Gemmini: calibration에 쓰지 않은 GEMM 및 convolution RTL workload
2. Eyeriss: 공개 silicon point가 추가로 없으면 동일 architecture의 mapper-independent analytical
   invariant 또는 별도 RTL/FPGA reference를 보조 근거로 사용
3. calibration result와 holdout result를 보고서에서 분리

### P2 — NVDLA RTL energy pipeline

1. layer measurement window의 FST/VCD/SAIF 수집
2. synthesized netlist와 cell-library power 분석
3. CBUF/CACC SRAM macro 또는 CACTI characterization
4. DBBIF transaction을 DRAMsim3에 replay
5. logic/clock/SRAM/DRAM hierarchy breakdown 생성
6. calibration workload와 독립 holdout workload 분리

### P3 — Gemmini absolute energy

1. Gemmini RTL activity 수집
2. 동일 technology/PVT에서 MAC, control, clock network 합성
3. 현재 unpriced 6개 event의 독립 cost 확보
4. CACTI/DRAMsim3와 합성 reference의 technology boundary 통일
5. 동일 GEMM의 RTL energy와 NPUsim total energy 직접 비교

### P4 — Eyeriss absolute energy

Eyeriss는 공개된 layer별 실측 joule 또는 동일 operating point의 component power가 확보되지 않으면
absolute validation을 완료하기 어렵다. 그 전까지는 `normalized` 상태를 유지하고 relative energy
결과를 absolute pJ/W로 변환하지 않는 것이 올바르다.

## 9. 최종 판정

### Cycle

- **Eyeriss:** 높은 calibrated fidelity, 5개 AlexNet convolution 한정
- **Gemmini:** 높은 calibrated fidelity, 6개 16×16 WS GEMM 한정
- **NVDLA:** bring-up 성공과 n=1 예비 결과만 존재, fidelity 판단 보류

### Energy

- **Eyeriss:** relative breakdown은 유의미하지만 absolute silicon energy는 미입증
- **Gemmini:** internal accounting과 일부 SRAM/DRAM component reference는 강하지만 RTL total
  energy는 미입증
- **NVDLA:** energy metric 수집 자체가 미착수

따라서 현재 NPUsim의 가장 강한 외부 정합성 근거는 **cycle 축의 Eyeriss silicon과 Gemmini RTL**이다.
Energy 축은 계산 구조와 일부 component model은 충분히 견고해졌지만, hardware-level absolute
fidelity를 주장할 단계에는 아직 도달하지 않았다.
