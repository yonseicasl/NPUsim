# 2026-08-19 latency 재평가 Phase 4 구현 기록 (L4 제외)

기준 문서: [../../assessment/TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md](../../assessment/TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md)
선행 기록: [Phase 1](latency_issues_2026-08-19_phase1_fixes.md),
[Phase 2](latency_issues_2026-08-19_phase2_fixes.md),
[Phase 3](latency_issues_2026-08-19_phase3_fixes.md)

Phase 4의 항목 12·13·14를 구현했다. **항목 11(L4, memory-inclusive `Critical-path latency`의
외부 reference 확보)은 측정 데이터가 없어 제외했다** — 코드 수정으로 닫을 수 있는 항목이 아니다.

## 용어 주의 — 두 개의 "phase" 번호

이 repo에는 서로 무관한 두 개의 phase 번호가 있다.

| 번호 체계 | 의미 |
|---|---|
| `validation/phase1` / `phase2` / `phase3` | **검증 대상**: SCALE-Sim v2 / Gemmini RTL / Eyeriss silicon |
| 본 문서의 Phase 1~4 | 평가 문서 「권장 수정 순서」의 **수정 우선순위 단계** |

대응하지 않는다. 예를 들어 본 문서의 Phase 3는 DRAM/systolic 상세 모델 작업이며
`validation/phase3`(Eyeriss silicon)와 무관하다.

## 항목 12 — acceptance gate 추가

### 12-A. SCALE-Sim v2 cross-simulator gate

`validation/phase1`은 첫 milestone부터 SCALE-Sim v2와 systolic array schedule을 비교해 왔지만
**사람이 읽는 리포트였을 뿐 자동으로 실행되지 않았다.** 우리가 가진 세 외부 reference 중 하나가
회귀 없이 방치돼 있었다는 뜻이다. gate로 만들었다.

SCALE-Sim의 `COMPUTE_REPORT.csv`가 `validation/phase1/out`에 커밋돼 있어 **SCALE-Sim 설치 없이**
동작한다(재생성 방법은 README에 있음).

비교 축은 문서화된 것 그대로다.

```text
array-sched = max(Computation cycle, PE-array input IC, PE-array weight IC)
```

OUTPUT stream은 제외한다 — NPUsim 매핑은 reduction step마다 partial sum을 write-through하고 WS
array는 in-array accumulate하는, 문서화된 모델링 차이다.

**regime별로 한계를 다르게 둔 이유.** 단일 tolerance는 전부 통과시키거나 timing 결함이 아닌 것을
실패시킨다. README의 regime 분석(Gemmini RTL로 확인됨)에 따라:

| regime | layer | 한계 | 근거 |
|---|---|---|---|
| FC | fc1~fc3 | **8%** | 양쪽 모두 fold fill/drain이 지배 — 같은 축의 진짜 timing 비교 |
| steady-state conv | conv3~5 | **20%** | SCALE-Sim 자체 fold 상수(R+C−2≈24)가 RTL 실측(~14)의 약 2배라 NPUsim이 체계적으로 낮게 읽힘 |
| mapping-limited conv | conv1~2 | **baseline 고정** | 순수 active-PE 산술(12×14는 7의 인수가 필요하나 AlexNet 차원에 7이 없음). 정확도로 gate하면 문서화된 표현력 한계를 gate하는 것이 됨 |

추가로 모든 layer를 `validation/phase1/npusim_baseline.csv`에 고정하고, Spearman rho ≥ 0.90을
요구한다.

현재 결과 (Phase 1~3 변경 후에도 문서화된 값과 **정확히 일치**):

```text
              fc: MAPE   6.0%  max   6.3%  limit 8%
     steady-conv: MAPE  15.0%  max  18.8%  limit 20%
 mapping-limited: MAPE  56.3%  max  81.4%  baseline-pinned
      all layers: MAPE  21.9%  Spearman rho 0.929 (limit 0.9)
```

구현: [validation/phase1/gate.py](../../validation/phase1/gate.py),
[validation/phase1/npusim_baseline.csv](../../validation/phase1/npusim_baseline.csv),
`unittest/run_timing_validation.sh`에 연결.

### 12-B. Bottleneck classification gate

외부 memory-inclusive latency reference가 없으므로(L4) **절대 critical path는 gate할 수 없다.**
reference 없이 검증 가능한 것을 gate한다: timeline이 bottleneck을 올바른 위치에 두는지, 그리고
**그것이 변경된 sub-model 때문인지**.

`gemmini.cfg`에서 각각 한 가지만 다른 세 config:

| config | 변경 | bottleneck | compute schedule |
|---|---|---|---|
| `gemmini` | 없음 | PE array (fabric-bound) | 3,518 |
| `gemmini_memory_bound` | DRAM 단가 ×10 | **DRAM** | **3,518 (불변)** |
| `gemmini_compute_bound` | MAC cycle 1→60 | **PE (compute+LB)** | 63,934 |

핵심은 **귀속**이다. memory-bound 변형의 compute schedule이 bit-identical하므로 재분류가 memory
model에서 왔음이 증명되고, compute-bound 변형의 모든 memory axis가 bit-identical하므로 그 반대도
증명된다. 세 config이 서로 다른 3개 regime을 덮으므로 두 label을 한 현상에 붙인 것이 아니다.

[validation/bottleneck/check.py](../../validation/bottleneck/check.py) BN1~BN6:

- BN1: 보고된 `Bottleneck level`이 실제로 busy가 가장 큰 stage인가 (자기 숫자와 어긋나는 분류는
  없는 것보다 나쁘다)
- BN2: critical path ≥ bottleneck busy, ≥ compute schedule
- BN3/BN4: 위의 귀속 논증
- BN5: 무변경 config은 둘 중 어느 쪽도 아님 (fabric-bound)
- BN6: multi-chip fixture의 compute-schedule identity

## 항목 13 — external traffic 범위 명시

평가 문서는 "오차 축소 **또는** 명시적 비지원으로 제한"을 요구했다. **후자를 택했다** — 오차를
줄이려면 batch filter reuse와 psum GLB spill/reload 모델링이 필요하고 이는 timing 수정이 아니라
functional/mapping 모델 확장이다.

### 문제였던 것

기존에는 `Eyeriss traffic: max=133.94%` 한 줄이었다. **어느 축이 왜 실패하는지 숨겼고**,
`--check-traffic`은 항상 실패하므로 아무 정보도 주지 못했다.

### 축별 분해 — 두 축이 서로 다른 원인이다

```text
 * GLB  MAPE= 30.50%  max= -81.62% (conv3)  ceiling 82%  OUT OF GATE SCOPE
        per layer: conv1 +6.9%, conv2 -11.1%, conv3 -81.6%, conv4 -26.5%, conv5 -26.4%
 * DRAM MAPE= 60.46%  max=+133.94% (conv5)  ceiling 134%  OUT OF GATE SCOPE
        per layer: conv1 +30.6%, conv2 +47.4%, conv3 +16.2%, conv4 +74.1%, conv5 +133.9%
```

- **133.94%는 DRAM 축(conv5)** 이다. DRAM 오차는 체계적 과대계상이며 chip이 활용한 batch reuse가
  클수록 커진다(conv3 +16% → conv5 +134%) — batch filter reuse 미모델링과 정확히 부합한다.
- **GLB 축은 conv3만 −81.6%** 로 다른 conv(−11~−27%)와 **다른 regime**이다. 균일한 모델 갭이 아니라
  그 layer의 매핑을 가리킨다 — 이후 reuse 모델링 작업자를 위한 구체적 단서로 기록했다.

### gate가 의미를 갖게 만든 방법

축별로 **regression ceiling**을 오늘 측정값에 고정했다. `--check-traffic`은

1. milestone 15% 안에 든 축은 15%로 gate하고,
2. 밖에 있는 축은 **더 나빠지지 않는지**를 gate한다.

따라서 flag는 항상 의미가 있고(무성 악화 불가), 15%는 영구적인 빨간불이 아니라 **도달 목표**로 남는다.
현재는 두 축 모두 milestone 밖이므로 gate 출력이 그 사실을 명시한다.

```text
--check-traffic gates: no axis is inside the milestone limit yet; both are held to their
                       regression ceilings only
```

검출력 확인: DRAM ceiling을 134→100으로 낮추면
`REGRESSION: 133.94% exceeds this axis's 100% ceiling`으로 exit 1. 원복 후 exit 0.

구현: [validation/check_timing.py](../../validation/check_timing.py) `TRAFFIC_AXES`.

## 항목 14 (L11) — 지원 layer 범위를 숫자 옆에 명시

### 문제였던 것

지원하지 않는 layer(pooling/activation/normalization)는 timing에서 제외되고 경고가 나왔지만,
그 경고는 `network.txt`의 **맨 끝(333행)** 에 있었다. 위쪽 `Layer timeline`의 latency를 읽는 사람은
그것이 partial인지 알 수 없었다 — "network end-to-end latency"로 오독될 수 있는 구조다.

### 수정

scope 문장을 **그것이 한정하는 숫자 바로 앞**으로 옮겼다.

```text
============ Layer timeline =============
Timing scope          : 8 of 13 layers  (PARTIAL: 5 unsupported layer(s) excluded; this is NOT an end-to-end network latency)
Compute-schedule latency : 23999585.3 cycles (validated metric)
```

제외된 layer의 시간은 **추정하지도 0으로 채우지도 않고 그냥 없다**는 점도 주석에 명문화했다.

### gate

`check_network_rollup()`이 다음을 검사한다.

- scope 문장이 존재하는가
- 선언한 layer 수가 실제로 rollup된 `layer_*.txt` 수와 일치하는가 (장식용 문장 방지)
- 제외가 있으면 `PARTIAL`로 표시되는가
- 전부 지원되면 `complete`로 표시되는가

기존 trailing 경고도 유지했다(하위 호환).

## 검증 결과

| 검증 | 결과 |
|---|---|
| Unit/config validation | 통과, 233개 config |
| **SCALE-Sim cross-simulator gate** | **통과 (신규)** — fc 6.0%/8%, steady-conv 15.0%/20%, rho 0.929/0.90 |
| Gemmini RTL timing | 통과, MAPE 4.40%, max 7.86% (baseline 불변) |
| Eyeriss silicon latency | 통과, MAPE 4.26%, max 6.39% (baseline 불변) |
| **External traffic 축별 scope** | **통과 (신규)** — 두 축 모두 regression ceiling 내 |
| **Network timing scope (L11)** | **통과 (신규)** |
| Traffic T1~T10 | 8개 workload 통과 |
| Multi-chip MC1~MC9 | 통과 |
| Buffering LB-A~LB-E | 통과 |
| Bypass BP1~BP6 | 통과 |
| Asymmetric stats/timeline (L10) | 통과 |
| DRAM DR1~DR6 | 통과 |
| Format FM1~FM6 | 통과 |
| **Bottleneck BN1~BN6** | **통과 (신규)** |
| Energy E1~E5 | 통과 |

## 남은 항목

### 코드로 닫을 수 없는 것
- **L4 / 항목 11**: memory-inclusive `Critical-path latency`의 외부 reference.
  측정 latency 데이터가 선행되어야 한다. 그때까지 `Critical-path latency`는 informational이고
  `Compute-schedule latency`만 외부 검증 대상이라는 분리를 유지한다(출력에 명시돼 있음).

### 모델 확장이 필요한 것
- **external traffic 오차 축소**: batch filter reuse, psum GLB spill/reload 모델링.
  현재는 축별 원인과 regression ceiling으로 범위를 고정해 두었다. GLB 축 conv3의 −81.6%가 다른
  conv와 다른 regime이라는 관찰이 착수 지점이다.
- **지원 layer 확장**: pooling/activation/normalization의 accelerator timing.
  구현되면 `Timing scope`가 `complete`로 바뀌고 그때 비로소 "network end-to-end latency"라고
  부를 수 있다.
- Phase 2~3에서 이미 기록한 항목: NoP router arbitration/virtual channel, per-chip mapping,
  DRAM address 기반 bank/row(선택지 2 = DRAMSim3 gate), systolic per-PE register 단위 event 모델.
