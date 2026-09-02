# Energy-model validation (2026-08-17)

## 방법

에너지 = (검증된 traffic/접근 카운트) × (단위 비용)이므로, 카운트 검증
([traffic](../traffic/README.md)) 위에 **공개된 단위-비용 기준**을 얹어 계층별
에너지 분해를 대조한다. `configs/accelerators/eyeriss_energy.cfg`는 단위 비용을
ISCA'16 Table IV의 정규화 접근 에너지(1 MAC 기준: RF 1×, array/GIN 2×, GLB 6×,
DRAM 200× per 16b word)로 설정 — timing에는 무영향(latency 불변 확인).

```bash
./npusim.sh run eyeriss_energy alexnet silicon
python3 validation/energy/check.py
```

## 결과

```
 layer  RF/MAC   ALU%    RF%   GIN%   GLB%  DRAM%  RF:rest
 conv1    5.18   11.0   57.3    3.7   11.3   16.7    2.20:1
 conv2    5.27   12.8   67.7    2.8    8.5    8.2    2.80:1
 conv3    5.15   13.5   69.3    2.3    7.0    7.9    3.04:1
 conv4    5.18   12.7   65.7    2.7    8.5   10.3    2.75:1
 conv5    5.18   12.4   64.1    2.7    8.4   12.5    2.74:1
   fc1    6.00    1.3    7.9    5.9   18.0   66.9    0.31:1
```

| 검사 | 기준 | 판정 |
|---|---|---|
| E1 | ALU 에너지 == MAC 수 × 1 | **정확** (conv3: 598,081,536) |
| E2 | RF 접근/MAC ≈ 4 (+spad refill) | **5.15~5.27** ✓ (per-MAC operand 이동이 이미 모델링됨) |
| E3 | CONV RF:rest(DRAM 제외) ≈ 4:1 (JSSC 칩 실측 주장) | **2.20~3.04** — 잔차는 PA9(GLB read 재스트림 2~2.6×)로 전량 귀속: 제거 시 ≈4.2 ✓ |
| E4 | FC는 DRAM 지배 (ISCA Fig.10) | **67~71% DRAM** ✓ |
| E5a | DRAM access 에너지/사이클 == 단가비(logical access 기준) | **정확** (전 8층) |
| E5b | GLB GIN 에너지 == serialized 링크 txn × 단가 | **정확** (전 8층) |

## 검증이 적발·해결한 이슈

1. **mc→GLB fill의 스케일 비일관 (해결)**: DRAM은 datatype별 반복 스케일(DR4)로
   fetch-once를 credit하면서, 그 공급을 받는 GLB **write(fill) 측은 uniform ×R**
   — 출처 없는 데이터가 GLB에 52× 쓰이는 자기모순. GLB에 `fill_access_*`
   분리 계정 신설(`multi_chip.cc`가 fill로 과금, stats가 per-type 스케일 후
   access 합산) → conv3 GLB weight 184MB→94MB, E3 1.78~2.57 → 2.20~3.04.
2. **GLB Interconnection(GIN) 에너지 미출력 (해결)**: `transfer_energy_global_buffer`
   수집만 되고 출력 누락 → GLB 섹션에 "Interconnection energy" 추가.
3. **(외부 assessment 적발, 해결) GLB fill cycle의 timeline 누락**:
   `modeled_elapsed_cycles()`가 fill 축을 무시해 `busy ≥ max(access)` 불변식
   위반(전 8층) → busy 축을 `access + fill` 합으로 보수화(스케일 전 합산이라
   per-type ≤ uniform 관계로 스케일 후에도 불변식 보장). 전 8층 성립 확인.
4. **(외부 assessment 적발, 해결) E5 정의 오류**: 기존 E5는 access 에너지를
   *serialized 링크* 트래픽과 비교해 packet padding이 있는 FC에서 −4~−10%
   편차 — logical access 항등식(E5a)과 링크 항등식(E5b)으로 분리.
5. **checker 판정 부재 (해결)**: E1~E5b에 assertion·비정상 종료 코드 추가
   (`ALL ENERGY CHECKS PASSED` / exit 1).

## 한계·후속

- E3 잔차 = PA9(P4): batch-인터리브 필터 공유 + psum 왕복 — GLB read 측 2~2.6×.
- 절대(pJ) 검증은 미수행: 칩 전력은 clock network 33~45%·zero-gating 포함으로
  단위-비용 캘리브레이션 없이 비교 불가 — 상대 분해 검증이 올바른 축.
- eyeriss.cfg의 GLB `bypass=0:1:0`은 descriptor 경로에서 미적용(용량 검사만) —
  소형 정합 이슈로 기록.
