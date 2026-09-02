# Traffic-counter validation (2026-08-17)

## 방법

`check.py`가 검증된 3개 워크로드군(eyeriss/alexnet silicon 5 conv +
gemmini GEMM 3점)에 대해 **해석적 항등식**을 시뮬레이터 카운터와 전수 대조한다.
단위는 config에서 유도(링크 트랜잭션 = 링크 bitwidth: eyeriss GLB 32B/DRAM 8B,
gemmini GLB 16B/DRAM 8B; 데이터 fp16 2B / int8 1B).

| 검사 | 항등식 | 결과 |
|---|---|---|
| T1 | output DRAM == 출력 볼륨 (단일-쓰기 관례) | **전 8케이스 1.000** |
| T2 | weight DRAM ≥ 볼륨 (비율=재인출 계수) | **전 케이스 1.00 = fetch-once** |
| T3 | input DRAM == 전체 mapped input union(halo 포함) | **전 8케이스 1.000** |
| T4 | GLB↔array ≥ DRAM (타입별) | 전 케이스 ok |
| T5 | GLB weight 바이트 == 이벤트×array tile (PA1 부기) | **전 케이스 1.000** |
| T6 | comp == MACs/activePE | **전 케이스 1.000** |
| T11 | 결과의 halo 계약·working set·전후 transaction == T3 | **conv 5개 applied, GEMM 3개 not needed** |

```bash
python3 validation/traffic/check.py   # ALL CHECKS PASSED
```

## 검증이 적발·해결한 이슈

1. **GB (해결)**: T5=0.500 전 케이스 → `global_buffer.cc` data_transfer의
   **WEIGHT legacy host-address timing 경로가 `#else`로 활성**이어서 descriptor
   계정과 이중 계상(카운트·GLB access cycle/energy). INPUT은 이미 `#if 0`
   처리되어 있었음 — WEIGHT도 동일 처리 → T5 = 1.000.
2. **DR6 (해결)**: T1 사전값 0.5~0.97 → 레이어 최종 output tile flush 미계상.
   `multi_chip::flush_output_writeback()`(live 루프 종료 시 1회, repetition
   스케일 전) 추가 → 전 케이스 1.000.
3. **E20-3 input halo (2026-08-20 해결)**: legacy GLB의 P/Q tile을 각각 독립
   fetch하던 input traffic을 전체 convolution window의 정확한 union으로 병합했다.
   단순한 반복계수의 소수 스케일이 아니라 padded P/Q·stride·filter extent를 직접
   계산하므로 conv1/2의 2차원 타일링과 마지막 transaction tail도 정확하다.
   구현은 다음 두 조건을 동시에 만족할 때만 적용된다.

   - GLB가 sliding ring working set을 실제로 수용할 것. shared GLB는 resident
     weight/output tile과 double-buffer copy까지 함께 예약한다.
   - legacy R/S factor가 1일 것. R/S 분할의 loop-order 의미론은 P/Q sliding과
     다르므로 검증되지 않은 경우에는 보수적으로 적용하지 않는다.

   logical tile request는 유지하되 중복 payload만 coalesce한다. DRAM source access는
   exact 정수 access 수×unit cost로 cycle/energy를 재구성하고, multi-chip·GLB fill은
   동일한 coalescing factor로 함께 줄인다. 결과 파일에는
   replicated→unique elements, working-set bytes, DRAM 전후 transaction을 남긴다.

## 외부 실측과 남은 차이

- **칩 Table V DRAM은 RLC 압축 후 수치**: ifmap zeros 72~79% 층에서 칩 값이
  dense 볼륨보다 작은 이유(예: conv5 칩 1.3MB ≈ dense in×(1−0.776)+wt+out).
  NPUsim은 dense 계상 — 비교 시 정의차로 참작해야 함.
- halo 적용 후 NPUsim input은 전층 정확한 dense union이다. RLC dense-equivalent
  chip과의 잔여 오차는 이제 NPUsim이 **더 작은 방향**이며, 이는 halo가 아니라
  공개되지 않은 silicon mapper의 추가 refetch/scheduling traffic 또는 Table V의
  dense-equivalent 변환 불확실성이다. 외부 DRAM 축은 따라서 계속 informational이다.
