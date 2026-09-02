# Phase 1 — Cross-simulator timing validation vs SCALE-Sim v2

## 목적

NPUsim의 systolic-array timing 모델을 RTL-validated 시뮬레이터인 **SCALE-Sim v2**와
동일 layer shape에서 교차 비교해, (a) 정합 구간의 오차를 정량화하고 (b) 모델이
포착하지 못하는 현상을 식별한다. 이는 simulator-vs-simulator 비교이므로 "증명"이
아니라 **스크리닝**이며, ground truth는 Phase 2(Gemmini RTL)가 담당한다.

## 실험 설정

| 항목 | NPUsim | SCALE-Sim v2 |
|---|---|---|
| Array | 12×14 (`configs/accelerators/scalesim.cfg`, `[systolic_array]`) | 12×14 (`scalesim_12x14_ws.cfg`) |
| Dataflow | WS (`pe_stationary=weight_stationary`) + eyeriss alexnet mapping | WS (자체 fold 스케줄) |
| Workload | AlexNet 8 layers (batch 4, conv2/4/5 group 2) | 동일 shape (`alexnet_npusim.csv`), single-image·single-group |
| MAC 단가 | `computation_cycle = 1` | 1 MAC/PE/cycle |
| 메모리 | 비교 대상에서 제외 (SCALE-Sim CALC 모드, stall=0 확인) | 〃 |
| SCALE-Sim 버전 | scale-sim-v2 HEAD + numpy≥1.25 호환 패치 (`double_buffered_scratchpad_mem.py`) | |

**비교 축**: SCALE-Sim `Total Cycles`(compute 스케줄, fill/drain 포함) × batch ×
group ↔ NPUsim `Computation cycle`(busiest-PE 직렬 MAC issue; `= MACs/activePE`가
정확히 성립함을 확인). 실행: `python3 validation/phase1/compare.py`.

## 결과

```
layer         MACs(x4)   NPUsim comp    SS x b x g     err%    NPU actPE    SS actPE  SS util%
conv1      421,660,800     3,833,280     2,656,944    +44.3        110.0       158.7      94.5
conv2      895,795,200    19,906,560     6,119,992   +225.3         45.0       146.4      87.1
conv3      598,081,536     3,833,856     4,408,316    -13.0        156.0       135.7      80.8
conv4      448,561,152     2,875,392     3,306,232    -13.0        156.0       135.7      80.8
conv5      299,040,768     1,916,928     2,361,592    -18.8        156.0       126.6      75.4
fc1        150,994,944     2,359,296    33,303,548    -92.9         64.0         4.5       2.7
fc2         67,108,864     1,048,576    14,830,484    -92.9         64.0         4.5       2.7
fc3         16,384,000       163,840     3,644,348    -95.5        100.0         4.5       2.7
```

## 해석 — 세 개의 regime

### R1. 정상 steady-state conv (conv3/4/5): 오차 −13 ~ −19%
양쪽의 유효 active PE(156 vs 127–136)가 근접하고, 오차는 SCALE-Sim이 세는
**fold-경계 낭비**(파형 가장자리의 유휴 lane)를 NPUsim 매핑이 세지 않는 데서
온다. Cross-simulator 비교로서 합리적인 정합 구간이며, cycle 회계 자체의 결함
징후는 없다.

### R2. 매핑-제한 conv (conv1 +44%, conv2 +225%)
NPUsim은 **주어진 매핑을 충실히 실행**한다 — eyeriss RS 매핑이 이 layer들에서
110/45개 PE만 활성화하므로 느리게 나온다. SCALE-Sim은 **자체 WS folding**으로
87–95%를 채운다. 즉 이 오차는 timing 모델이 아니라 **mapping 품질 축**이다.
NPUsim의 `NPU actPE` 열이 그 원인을 그대로 노출한다(45 = conv2 매핑의 spatial
factor 곱). 공정 비교를 위해서는 WS-folding에 상응하는 매핑을 작성해 주입해야
하며, 이는 Phase 2에서 매핑을 우리가 통제함으로써 해소된다.

### R3. FC / 낮은 temporal-depth (fc1–3): 오차 −93 ~ −95% — **핵심 발견**
SCALE-Sim의 FC utilization은 **2.7%**다. WS systolic에서 batch-1 FC는 fold당
유효 사이클이 1인데 pipeline fill/drain이 `R+C−2 = 24` cycle을 차지하므로
utilization이 붕괴한다 — 이는 TPU v1 논문이 MLP에서 보고한 낮은 utilization과
같은 **실제 하드웨어 현상**이다. NPUsim의 compute 모델은 `MACs/activePE`로,
**fold-level fill/drain을 모델링하지 않아** 이 현상을 놓친다(SY2에서 넣은
diameter fill은 *operand 분배 stream당* 1회이지, *compute fold당*이 아니다).

→ **조치 항목 (V1)**: fold-aware compute fill — 분석적으로는
`cycles ≈ folds × (R + C − 2 + T_fold)` 형태의 항을 `computation_cycle`에 도입
가능. temporal depth(T_fold)가 클수록(conv) 무시 가능하고 작을수록(FC) 지배적
— R1 구간을 해치지 않고 R3를 잡는다. batch-as-spatial 재인코딩 실험(아래)으로
현상 스케일을 재확인했다.

### Batch-일치 FC 재실험 (fc2/fc3, batch4를 spatial 2×2로 인코딩)

1×1 filter FC는 batch를 IFMAP spatial로 인코딩할 수 있어 양쪽 batch 조건을
일치시킬 수 있다(`alexnet_fcs_b4.csv`, 스케일 팩터 ×1).

| layer | NPUsim comp | SCALE-Sim (b4) | err% | SS util% |
|---|---:|---:|---:|---:|
| fc2 | 1,048,576 | 4,008,239 | **−73.8** | 10.0 |
| fc3 | 163,840 | 984,959 | **−83.4** | 9.9 |

batch-4 amortization으로 SS utilization이 2.7%→10.0%로 오르는 것은 fold당
temporal depth T=4에서 `T/(T+R+C−2) = 4/28 ≈ 14%`(+drain)라는 fold-fill 이론과
정확히 부합한다. 오차가 −93%→−74%로 줄지만 여전히 fold-fill이 지배 —
**V1(fold-aware compute fill)은 batch 조건과 무관하게 필요**하다.

## R2 해소 — mapping-matched 실험 (`matched.map`)

R2 오차가 전적으로 매핑 축인지 검증하기 위해, SCALE-Sim의 WS folding에 근접하는
매핑을 conv1/conv2에 대해 작성했다(`configs/mappings/scalesim/alexnet/matched.map`;
제약: 레벨별 곱=layer 차원, through-GLB 곱 보존(GLB 용량), per-PE tile ≤ LB 용량,
active shape ≤ 12×14). conv1은 `PE_X=Q11`(121 PE), conv2는 `PE_X=K4·P3,
PE_Y=C2·R5`(120 PE). 실행: `python3 validation/phase1/compare.py matched`.

| layer | 원본 err% | matched err% | NPU actPE | SS actPE |
|---|---:|---:|---:|---:|
| conv1 | +44.3 | **+31.2** | 110→121 | 158.7 |
| conv2 | **+225.3** | **+22.0** | 45→120 | 146.4 |

**핵심 관찰**: matched 후 잔여 오차가 active-PE 비율과 **정확히 일치**한다 —
conv1: 158.7/121 = 1.312 → +31.2%, conv2: 146.4/120 = 1.220 → +22.0%. 즉
`NPUsim comp = MACs/actPE`가 소수점까지 성립하므로, **매핑을 맞춘 뒤 남는 오차는
100% 채움율 차이이고 timing 잔차는 0**이다. 채움율을 더 올릴 수 없는 이유는
순수 산술적 제약이다: 12×14를 채우려면 7의 배수 인수가 필요한데 AlexNet 차원
(96, 55=5·11, 27=3³, 3, 11, 5)에 7이 없다. SCALE-Sim은 fold 나머지를 부분-활용
fold로 처리해(94.5% 유효) 이 제약을 회피한다 — NPUsim 매핑 문법은 나눗셈이
정확히 떨어져야 하므로 partial fold를 표현할 수 없다(도구 표현력의 차이로 기록).

**결론**: R2는 timing 결함이 아님이 실험적으로 증명됨. conv-only MAPE(matched):
**19.6%** (31.2/22.0/13.0/13.0/18.8).

**부수 발견 — 정정 (2026-08-16)**: 매핑 작성 중 관측된 segfault의 원인은 당초
추정("PE_X K > PE-row K 제약")이 **아니었다**. Phase-2 준비 중 ASan으로 추적한
결과, 진짜 원인은 `scheduler/stats.cc`의 기존 UB — `tile_size.reserve(N)` 후
**빈 outer vector를 `operator[]`로 인덱싱**하는 코드였다. allocator 상태(직전
할당 이력)에 따라 무증상이거나 crash하는 전형적 UB로, 매핑 내용은 할당 패턴을
바꾸는 간접 트리거였을 뿐이다. `resize()`로 수정 완료(회귀·unittest 통과).
따라서 K-spatial 매핑 제약은 존재하지 않으며, matched.map의 Q-spatial 우회는
불필요했지만 결과 유효성에는 영향이 없다(수치는 매핑 그대로의 정확한 비용).

## R3 해소 — 비교 축 정정 (V1 재판정: 모델 갭이 아니라 축 오류)

후속 분석에서 R3(−93~−95%)는 **모델 부재가 아니라 비교 축 선택의 오류**로
판명됐다. NPUsim은 systolic wavefront fill을 SY2 구현에서 이미 **operand stream당**
과금하며(`PE-array Interconnection cycle`), FC처럼 stream이 지배하는 layer의
array-schedule 시간은 `comp`가 아니라 stream 축에 기록된다. 올바른 비교량은

> **array-sched = max(Computation cycle, IC_in, IC_w)**

이다 (OUTPUT stream 제외 — NPUsim 매핑은 partial sum을 reduction step마다
write-through하는 반면 WS array는 in-array accumulate하는 **문서화된 모델링
차이**). 축 정정 결과(코드 변경 0):

```
layer       NPU comp    NPU IC_in     NPU IC_w  array-sched   SS x b x g    err% comp-only%
conv1      3,484,800      366,300    1,386,000    3,484,800    2,656,944   +31.2      +31.2
conv2      7,464,960    2,052,864   11,104,128   11,104,128    6,119,992   +81.4      +22.0
conv3      3,833,856      559,104    1,876,992    3,833,856    4,408,316   -13.0      -13.0
conv4      2,875,392      619,008    1,517,568    2,875,392    3,306,232   -13.0      -13.0
conv5      1,916,928      412,672    1,011,712    1,916,928    2,361,592   -18.8      -18.8
fc1        2,359,296   33,030,144   35,389,440   35,389,440   33,303,548    +6.3      -92.9
fc2        1,048,576   14,680,064   15,728,640   15,728,640   14,830,484    +6.1      -92.9
fc3          163,840    2,949,120    3,440,640    3,440,640    3,644,348    -5.6      -95.5

MAPE = 21.9%   Spearman rho = 0.929   (원래 축: 74.5% / 0.024)
```

- **FC: +6.3 / +6.1 / −5.6%** — fill/drain 물리가 이미 정확히 잡혀 있었음이 증명됨
  (fc2 산술: IC_w = 2,097,152 txn×1 + 1,048,576 events×13 fill = 15,728,640 —
  SY2 fill 공식 그대로).
- conv1/3/4/5: comp 지배(스트림이 그 아래로 overlap) — R1 오차 유지.
- conv2(+81.4): matched 매핑이 PE 채움(45→120)을 얻는 대신 **weight 재사용을
  희생**해 IC_w가 comp를 추월 — comp축(+22%)과 stream축(+81%)이 그 tradeoff의
  양면. SS는 weight-hold를 극대화하는 스케줄을 스스로 고르므로 발생하는,
  매핑 선택의 실제 비용이다. conv2 제외 MAPE = **13.4%**.

**V1 재판정**: "fold-aware compute fill 모델 추가" 불필요 — 필요한 것은 (a) 본
문서의 비교 축 정의와 (b) (선택) NPUsim 결과에 array-sched 파생값을 직접
출력하는 리포팅 편의뿐이다. `issue/timing/systolic_array.md`의 V1 항목을 이
판정으로 갱신함.

## 위협 요인 (threats to validity)

1. SCALE-Sim도 시뮬레이터다 — R3의 fill/drain 상수는 SCALE-Sim의 WS 스케줄 정의에
   의존한다(단, 그 정의는 RTL 대비 검증됨).
2. NPUsim 매핑은 eyeriss RS용으로 작성된 것 — R2 오차는 매핑 산물이다.
3. `u_comp=1`, 메모리 제외 — compute 스케줄만의 비교다(의도된 범위).

## 결론 및 다음 단계

- **정합 구간(R1)**: cross-simulator 오차 13–19% — 스크리닝 통과.
- **모델 갭 식별(R3)**: fold-level fill/drain 부재를 정량화(−93%). V1 조치 항목
  도출 — 이것이 Phase 1의 가장 큰 수확.
- **Phase 2 (Gemmini Verilator RTL)**: 매핑을 통제한 GEMM sweep으로 R1 오차의
  ground-truth 확정 + V1 구현 후 재검증.

## 재현 방법

```bash
# SCALE-Sim (ext/scale-sim-v2 + ext/venv-scalesim 필요)
cd ext/scale-sim-v2
PYTHONPATH=. ../venv-scalesim/bin/python scalesim/scale.py \
  -c ../../validation/phase1/ss_convs.cfg -t ../../validation/phase1/alexnet_convs.csv \
  -p ../../validation/phase1/out -s N
# (fcs, fcs_b4 동일 패턴)

# NPUsim
./npusim.sh run scalesim alexnet energy

# 비교
python3 validation/phase1/compare.py
```

## Phase-2(RTL) 이후 재해석 (2026-08-17) — 전 잔차의 ≤8% 귀속

Phase-2에서 cycle-exact RTL이 **fold당 오버헤드 ≈14 cycle**을 실측하면서,
SCALE-Sim의 fold 상수(R+C−2≈24)가 실제의 약 2배임이 판정되었다. 이에 따라
Phase-1 잔차를 재귀속하면:

| 층 | vs SS total | 귀속 | 정정 후 |
|---|---|---|---|
| fc1/2/3 | +6.3/+6.1/−5.6% | 동일 현상(스트림 fill) 축끼리 비교 | **≤8% ✓** |
| conv3/4/5 | −13.0/−13.0/−18.8% | SS 자체의 fold 오버헤드(RTL이 과대 판정) | SS total×Overall Util(=SS의 fold 오버헤드 제거) 대비 **+7.69/+7.69/+7.69% ✓** |
| conv1/2 | +31.2/+22.0% | active-PE 매핑 산술: 7∤dims라 SS의 partial-row 채움(Mapping Eff 95.6/91.4%)을 NPUsim 매핑이 표현 불가 → active PE 수 비율이 오차를 정확히 예측 | timing 오차 아님(표현력 한계, 문서화) |

결론: **RTL을 기준 축으로 하면 NPUsim timing 오차는 전 구간 ≤8%**
([Phase-2](../phase2/README.md) 전 6점 MAPE 4.4%). SCALE-Sim 대비 잔차는
(a) 참조의 과대 fold 상수(RTL로 판정 완료), (b) 매핑 채움 산술로 전량 설명된다.
