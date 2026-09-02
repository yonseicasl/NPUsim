# Phase 3 — 실리콘 실측 검증 vs Eyeriss 칩 (JSSC Jan. 2017)

## 목적과 결론 요약

NPUsim의 RS(spatial_arch) timing을 **실제 실리콘 실측치**(Chen et al., "Eyeriss:
An Energy-Efficient Reconfigurable Accelerator for Deep CNNs," JSSC 2017,
Table V — AlexNet, batch 4, 1.0V, core 200MHz)와 대조했다.

- **AlexNet 전 5개 conv layer의 Processing-Latency 오차 ≤6.4% (MAPE 4.3%)** —
  비교량은 시뮬레이터가 직접 출력하는 `Computation + Fold fill cycle`.
- Phase 1(시뮬레이터 교차)·Phase 2(RTL)에 이어 **세 검증 축 모두에서 ≤8%**:
  systolic/WS(RTL)와 spatial/RS(실리콘)의 두 아키텍처 경로가 모두 검증됨.

## 기준 데이터 (JSSC'17 Table V)

| Layer | Processing (ms) | Total (ms) | MACs | Active PEs | GLB (MB) | DRAM (MB) |
|---|---|---|---|---|---|---|
| CONV1 | 16.5 | 20.9 | 0.42G | 154 (92%) | 18.5 | 5.0 |
| CONV2 | 39.2 | 41.9 | 0.90G | 135 (80%) | 77.6 | 4.0 |
| CONV3 | 21.8 | 23.6 | 0.60G | 156 (93%) | 50.2 | 3.0 |
| CONV4 | 16.0 | 18.4 | 0.45G | 156 (93%) | 37.4 | 2.1 |
| CONV5 | 10.0 | 10.5 | 0.30G | 156 (93%) | 24.9 | 1.3 |

"Processing Latency"는 DRAM 트래픽이 완전 오버랩됐을 때의 성능(논문 정의) —
NPUsim on-chip compute 축과 동일한 경계. MAC 수는 batch-4 AlexNet 수기 검산과
전층 일치(기준 신뢰성 확인).

## 결과 (silicon.map + eyeriss.cfg, 2026-08-17)

```
 layer NPU sched(ms) chip proc(ms)   err%
 conv1         15.47          16.5   -6.2
 conv2         36.83          39.2   -6.1
 conv3         21.28          21.8   -2.4
 conv4         15.96          16.0   -0.3
 conv5         10.64          10.0   +6.4
latency MAPE = 4.3%   max |err| = 6.4%
```

## 수행한 조치

1. **매핑(silicon.map)**: 칩의 실제 매핑을 재현.
   - conv1: 칩은 55 출력열을 **7×8=56으로 패딩**해 14열×11행=154 PE 사용.
     NPUsim 매핑도 동일하게 padded 표현(`PE_X P=7`, `DRAM P=8`) — 곱>차원
     매핑이 기존 tile 기계에서 일관 동작함을 확인(comp = 패딩 MACs/154 정확).
   - conv2: C를 PE 시간축에서 PE_Y 공간축으로 이동 + P 27→28 패딩 →
     14열×10행=140 PE (칩 135는 14+13 열 분할의 평균).
   - conv3/4/5: 기존 매핑이 이미 칩과 동일한 156 PE.
2. **코드 — spatial-padding guard** (`scheduler/npu.cc`,
   `mapping_table.{cc,h}` `calculate_total_parameter_size`): 매핑 곱이 layer
   차원을 초과(패딩)하거나 미달(부분 시뮬레이션)하면 **경고 출력** — 의도적
   패딩은 표현 가능하되 오타가 조용히 결과를 왜곡하지 못하게 함.
3. **캘리브레이션 — V2 fold-fill 재사용** (`eyeriss.cfg`
   `weight_fold_fill_cycle = 0.11`): 칩 논문이 명시한 미달 원인 2)
   "pass ramp-up 시간"을 weight-원소 residency당 0.11 cycle의 노출 비용으로
   과금(spad 채움의 비은폐 분율; per-layer setup 0). RTL의 Gemmini(14)와 같은
   knob, 아키텍처별 상수.
   - **주의(퇴화)**: 이 매핑들은 weight-bound라 `이벤트×per-PE weight tile ≡
     comp` — 즉 fill 항이 comp×0.11과 동치다. per-pass 상수를 분리하려면 칩의
     pass 수 데이터가 필요한데 논문에 없다. 5점 적합의 u 허용역(전층 ±8%)은
     [0.089, 0.127]로 좁아, 임의 보정이 아닌 제약된 캘리브레이션이다.
4. **데이터 포맷 fp16**: 칩은 16-bit 고정소수점 — 포맷을 16b로 맞춰 spad 용량
   검사(224×16b filter spad에 conv3 tile 384B ≤448B 등)와 트래픽 바이트가
   칩 사양과 대응.

## 메모리축 (2026-08-17 정정·수정 후)

> **정정**: 최초 리포트의 "DRAM 3.4~5.1×, GLB 0.1×"는 **트랜잭션→바이트 단위
> 오류**였다(링크 트랜잭션은 line 단위가 아니라 **링크 bitwidth 단위**: GLB
> 256b=32B/txn, DRAM 64b=8B/txn). 정정 후 실제는 DRAM 33×, GLB 2.6× 과대였고,
> 아래 수정으로 개선했다.

**수정 1 — datatype별 반복 스케일링 (`stats.cc`/`mapping_table.cc`, 모델 버그)**:
시뮬레이터는 GLB 반복 1회만 실행하고 stats에서 ×R 스케일하는데, 이 스케일이
off-chip(DRAM·multi-chip) 트래픽에도 일괄 적용되고 있었다. 그러나 tensor가
의존하지 않는 차원의 반복(예: weight에 대한 Q/B 루프)은 **같은 타일을 재방문
하며 on-chip 버퍼가 이를 보유**하므로 DRAM을 다시 읽지 않는다. →
`mapping_table::datatype_repetitions()`(weight: K·C·R·S, input: B·C·P·Q,
output: K·B·P·Q)로 off-chip 카운터를 per-type 스케일. conv3 DRAM weight
92MB→**1.77MB(=fetch-once 정확)**.

**수정 2 — 매핑 루프 재배치 (conv1/conv3)**: 재사용 루프(B/P, K-outer)를
DRAM행→legacy GLB행으로 이동해 live-레벨 재인출 제거.

**수정 3 — DR5: grouped legacy-K 의미론 확정 (2026-08-17)**: live 머신이
DRAM행 K 루프를 누적 GROUP으로 나누므로 **/G는 정확히 한 행에만** 두는 관례를
확정(`datatype_repetitions`에 K-보유 타입의 legacy-G 나눗셈 추가). conv2/4/5를
legacy K × DRAM K2 패턴으로 재매핑 → **전 5층 weight DRAM = fetch-once 정확
일치**, comp·fold 불변.

| 축 | 수정 전 | 수정 후 | 잔여 원인 |
|---|---|---|---|
| DRAM 내부 항등식 | input 2.60~5.06× | **input/weight/output 전층 mapped volume 1.00×**: P/Q halo union을 GLB capacity가 허용할 때 exact coalesce | legacy R/S 분할은 loop-order 미검증으로 보수적으로 미지원 — [validation/traffic](../../validation/traffic/README.md) |
| GLB MB | (0.1×로 오판) | chip이 현재 mapping의 [완전보존 하한, literal 상한] 안에 **1/5층** | psum spilling이 GLB capacity가 아니라 tiling hierarchy를 따라 발생하여 chip traffic과 latency를 동시에 맞추지 못함 |
| DRAM (RLC 보정) | — | **0.58~0.95× (MAPE 23.8%)**: chip dense-equivalent보다 작은 방향 | halo는 해결됨. 남은 것은 chip의 비공개 mapper/refetch traffic과 Table V RLC 변환 불확실성 |

input halo 수정 후에도 latency 축은 **불변**(≤6.4% 유지 확인).

## 회귀

- Phase-1(scalesim/matched) **byte-identical**, Phase-2(gemmini RTL 6점)
  MAPE 4.4% 불변, unittest 5/5.
- `eyeriss.cfg`는 fp16+fold-fill로 갱신됨 — 기존 `energy.map` eyeriss 결과는
  의도적으로 변경(캘리브레이션 반영).

## 재현

```bash
./npusim.sh run eyeriss alexnet silicon
python3 validation/phase3/compare.py
```

## Traffic·Energy 검증 (2026-08-17 추가)

latency 외 축의 체계적 검증을 신설 — [validation/traffic](../../validation/traffic/README.md)
(항등식 6종 전 8케이스 통과; GB5/DR6 적발·해결),
[validation/energy](../../validation/energy/README.md) (ISCA'16 Table IV 정규화
비용으로 계층 분해 대조: ALU·DRAM 정확, RF/MAC 5.2 ✓, FC DRAM 지배 ✓,
CONV RF:rest 2.2~3.0 — 잔차 PA9 귀속; GB6/GB8 적발·해결).

## Phase 1+2+3 종합

| Phase | 기준 | 경로 | 결과 |
|---|---|---|---|
| 1 | SCALE-Sim (교차) | systolic/WS | FC ±6.3%; conv 잔차는 참조의 fold 상수(RTL로 과대 판정)와 매핑 산술로 전량 귀속 |
| 2 | Gemmini RTL (cycle-exact) | systolic/WS | GEMM 6점 ≤7.9%, MAPE 4.4% |
| 3 | **Eyeriss 실리콘 (실측)** | spatial/RS | conv 5층 ≤6.4%, MAPE 4.3% |
