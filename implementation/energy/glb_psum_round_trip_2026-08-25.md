# 2026-08-25 GLB psum round trip 수정 (E20-3b)

E20-3의 남은 절반. `validation/phase3/compare.py`가 4/5 layer에서 "chip GLB traffic이 model의
literal-streaming 상한보다 위"라고 FAIL을 내던 축을 닫았다.

## 진단 — 가설이 아니라 계측으로

기존 주석은 원인을 "psum spilling이 GLB capacity가 아니라 tiling hierarchy를 따른다"로 적고 있었다.
capacity 쪽을 먼저 의심할 근거는 약했다: conv3의 GLB utilization은 **1.8%**다. capacity가
binding이면 traffic은 늘어야 하는데 model은 chip보다 **적게** 옮기고 있었다.

그래서 counter를 직접 갈랐다. `[shared] write_energy`의 output 성분을 3.00 → 3000.00으로 흔들고
energy delta에서 access 수를 역산했다(read 3.06은 고정이므로 두 미지수가 분리된다):

```
baseline  R*3.06 + W*3.00 = 101,300,060.16
probe     R*3.06 + W*3000.00 = 98,125,929,308.16
        -> W = 32,707,584 accesses,  R = 1,038,336 accesses
```

tile 832 element, line 8 bit -> event당 1664 access이므로

| conv3 (layer_4) | event 수 |
| --- | ---: |
| psum store (array -> GLB) | **19,656** |
| psum load (GLB -> array) | **312** |
| 최종 off-chip store를 위한 GLB read | 312 |

63배 비대칭이다. **써낸 psum을 읽어오지 않으면 누적이 성립하지 않는다.** 이건 capacity 문제가
아니라 왕복의 한쪽만 세고 있던 문제다.

## 원인

`global_buffer_t::data_transfer()`의 load 쪽이 store 쪽과 **다른 판단**을 쓰고 있었다.

* store: `pe_array_t`의 `equal_output_tile` — `tiles_suggest_retention && psum_retention_valid`.
  세 갈래로 multi-chip/DRAM output tile과 weight tile까지 본다. E20-3이 여기를 고쳤다.
* load: `skip_transfer[OUTPUT]`를 **tile 크기 비교만으로 latch**. 게다가 layer 안에서 다시 false로
  돌아오지 않는다. 한 번 걸리면 이후 모든 reload가 사라진다.

E20-3이 store 쪽만 고쳐 두 판단이 서로 모순된 상태로 남아 있었다.

## 수정

두 다리를 **하나의 판단**으로 묶었다. 별도 술어를 유도하지 않고 store 쪽 판정을 그대로 비춘다:

```cpp
// components/global_buffer.cc
skip_transfer[data_type_t::OUTPUT] = pe_array->equal_output_tile;
```

처음에는 latch 조건에 `psum_retention_valid`만 AND했는데, 그건 tile 비교를 남겨두는 것이라
512x512x512 GEMM에서 store/load가 여전히 어긋났다(아래 OPEN). tile 비교가 `tiles_suggest_retention`
보다 약하다는 걸 T12가 잡아내서 지금 형태로 바꿨다.

## 결과

| layer | GLB upper (전) | GLB upper (후) | chip | band |
| --- | ---: | ---: | ---: | --- |
| conv1 | 19.8 | **22.1** | 18.5 | ok -> ok |
| conv2 | 69.0 | **103.8** | 77.6 | FAIL -> **ok** |
| conv3 | 40.9 | **73.1** | 50.2 | FAIL -> **ok** |
| conv4 | 27.5 | **43.3** | 37.4 | FAIL -> **ok** |
| conv5 | 18.3 | **28.9** | 24.9 | FAIL -> **ok** |

conv3의 73.1 MB는 수정 전에 손으로 예측한 값과 정확히 일치한다(store 19,656 + load 19,656 =
39,312 event x 52 link x 32 B + input 7.67 MB).

**compute-schedule latency는 전혀 움직이지 않았다** — Eyeriss MAPE 4.26% / max 6.39%,
Gemmini RTL 4.40% / 7.86%, SCALE-Sim gate PASSED. psum traffic은 검증된 지표를 건드리지 않는다.
DRAM 축도 불변(MAPE 23.80%).

## 회귀 고정 3중

1. **`GLB_BAND_CEILING_LAYERS` 4 -> 0.** CI(`--check-traffic`)가 이 상한을 강제하므로 어느 layer든
   band 밖으로 나가면 실패한다.
2. **`Psum round trip : N loads / N stores`** 를 GLB 리포트에 노출. 숫자가 아니라 메커니즘이 보인다.
3. **traffic gate T12** — `loads <= stores <= loads + (output tile 수)`. 여유분은 output tile당 마지막
   store 하나(그건 reload 대신 off-chip으로 나간다). 예전 상태(19,656:312)를 압도적으로 잡는다.

`validation/phase3/npusim_baseline.csv`도 재생성했다. 이 파일은 halo 수정과 E20-3 시점부터 이미
stale이었다(conv3 glb 9.225216 = E20-3 이전 값, conv1 dram 6.532416 = halo 이전 값).

## 같이 바로잡은 것

`validation/phase3/traffic_reference.py`의 docstring이 아직 **반증된 주장**을 담고 있었다 — "upper
bound 위의 chip 값은 model이 아예 갖지 못한 traffic SOURCE의 증거". compare.py와 check_timing.py는
2026-08-20에 고쳤는데 정작 둘이 single-source라고 선언하는 모듈이 남아 있었다. 두 번의 반증
(tiling probe, 그리고 이번의 half-counted round trip)을 적어 넣었다.

## OPEN — T12가 새로 찾은 것

**gemm512x512x512: 0 loads / 96 stores.** reduction 차원이 array *위*가 아니라 **GLB 레벨 자체**에
있을 때(`GLB = ..., C 32, ...`) load 경로에 진입 자체를 하지 않는다. `skip_transfer`를 강제로 false로
두는 probe를 해봐도 0이므로 원인은 `output_read_global_buffer`의 offset 처리이고, 모든 방문이
first touch로 취급된다. E20-3b가 고친 범위(reduction이 array 위 + GLB 레벨이 output 인자를 가짐)와
다른 경우다. Eyeriss CONV 5개와 gemm64/gemm256은 모두 해당 없음.

T12에 **관측값을 못박아** known-open으로 등록했다(`PSUM_OPEN = {"gemm512x512x512": (0, 96)}`).
어느 쪽 숫자든 변하면 gate가 실패하므로 조용히 흘러가지 않는다.

## 이 세션에서 별개로 확인한 것

`validation/memory_concurrency`의 MC2/MC5가 실패한다(`expected DRAM->Multi-chip depth 128, got 1`).
**이 수정과 무관하다** — 패치를 되돌리고 재빌드해도 동일하게 실패하는 것을 확인했다. 별도 사안.

---

# 후속 (같은 날) — OPEN 종결과 MC2/MC5

## MC2/MC5 — knob이 아예 없었다

게이트와 config 3개는 있는데 **모델에 `max_outstanding_requests`가 존재하지 않았다.** 알 수 없는 키는
파서가 조용히 무시하므로 dead가 아니라 아예 없는 상태였다. `[multi_chip]`에 knob을 추가하고
(`0` = 미설정, 기존 double_buffer 유도 depth 유지) DRAM<->multi-chip boundary depth를 직접 덮도록
`stats_t`에 실어 보냈다. `pipeline_timeline_cycles()`는 이미 임의 depth N을 지원하고 있었다.

| config | depth | DRAM stall | critical |
| --- | ---: | ---: | ---: |
| gemmini_memory_bound | 1 | 1,152.0 | 384,942.0 |
| ..._outstanding1 | 1 | 1,152.0 | 384,942.0 |
| ..._outstanding | **128** | **0.0** | **383,790.0** |

MC1~MC5 전부 통과. critical 감소 1,152가 stall 감소와 정확히 일치한다(MC5의 산술 검사).

## gemm512 OPEN — 규칙이 틀렸고, 그 다음 셋이 더 나왔다

T12가 잡은 `0 loads / 96 stores`를 파고들자 결함이 연쇄로 나왔다.

**1. retention 규칙이 GLB row를 못 봤다.** `reduction_tiled_above_array()`가
`mapping_table[GLOBAL_BUFFER]`를 읽는데 그 행은 1로 채워져 있고 실제 GLB 인자는
`legacy_global_buffer_mapping`에 있다. **검증된 GEMM은 전부 reduction을 GLB에 매핑**하므로 전부
"array 위에 reduction 없음"으로 읽혔다.

**2. 규칙 자체가 틀렸다.** reduction이 array 위에 있다는 것만으로는 psum이 나가지 않는다. 나가게
만드는 것은 **그 reduction 안쪽에 output loop가 있어서** array가 다른 출력 타일을 돌다 돌아오는
경우다. gemm512의 GLB는 C=32뿐이고 그 아래 output 인자가 없으므로 array는 그냥 계속 누적하면 된다.
`psum_must_leave_array()`로 이름과 정의를 바꿨다 — 안쪽 레벨부터 바깥으로 훑으며 "이 reduction 안에
output loop가 있는가"를 한 번에 판정한다(같은 레벨 공존은 순서를 모델링하지 않으므로 보수적으로 취급).

**3. tile 크기 3분기가 남은 오차였다.** stationary type별 세 갈래와, 그 중 하나의
`weight tile이 다르면`이라는 의미 불명 대리 조건. 이제 두 조건뿐이다 — GLB output tile == array
output tile, 그리고 `psum_retention_valid`. 둘 다 loop nest에 대한 질문이다.

**4. residency latch가 풀리지 않았다.** retention이 걸리면 array가 **다른 출력 타일로 옮겨가도**
래치가 유지돼, 첫 완성 타일만 쓰고 이후를 전부 삼켰다. GLB가 새 출력 offset을 시작할 때 풀도록 했다
(`dram_t`가 multi-chip에 하는 것과 같은 방식).

**5. 마지막 타일이 나가지 않았다.** in-loop write-back은 타일이 **교체될 때** 발생하므로 layer가 끝날
때 array가 들고 있는 타일은 영원히 안 나간다 — DR6이 한 단계 위에서 닫은 것과 같은 구멍.
`pe_array_t::flush_psum_writeback()`을 추가하고 회계 본문을
`account_psum_store_to_global_buffer()`로 공유해 두 경로가 어긋날 수 없게 했다.

**6. store 경로가 bypass를 무시했다.** GLB `bypass=1:1:1`인데도 GLB access를 과금했다. load 경로는
P1-B 이후 지키고 있던 규약을 store 쪽만 안 지키고 있었고, in-loop store가 그 fixture에서 안 뜨는
바람에 숨어 있었다. flush가 뜨자 `energy_schema` ES4가 잡았다.

**7. 스케일링이 reduction 깊이를 곱하고 있었다.** psum이 array에 남으면 이 경계를 건너는 것은
**완성된 타일**이고 서로 다른 출력 타일 수만큼만 반복된다. uniform GLB 인자로 곱하면 reduction
깊이만큼 부풀려진다(gemm512에서 8배). RE1이 off-chip cast에 이미 적용한 논거 그대로
`m_datatype_repetitions[OUTPUT]`을 쓴다. **전송의 양끝**(GLB 쪽 열과 PE-array 쪽 read)에 같은
인자를 걸어야 하며, PE-array의 OUTPUT 열은 누적 write와 store read가 섞여 있어 store 몫을
`psum_store_access_*`로 따로 집계했다. 한쪽만 고쳤을 때 `array_buffer` AB2가 4배:1배 불일치를 잡았다.

### 결과

| case | 전 | 후 | 물리적 기대 |
| --- | --- | --- | --- |
| gemm64 / gemm256 | 0 / 0 | 0 loads / **1** store | 출력 타일 1개 |
| gemm512 | 0 / 96 | 0 loads / **4** stores | 출력 타일 4개 |
| conv1 | 10560 / 10560 | 10560 / **15840** | +5280 = 출력 타일 수 |
| conv2 | 79488 / 81216 | 79488 / **82944** | +3456 |
| conv3 | 19656 / 19656 | 19656 / **19968** | +312 |
| conv4 | 19344 / 19656 | 19344 / **19968** | +624 |
| conv5 | 12896 / 13104 | 12896 / **13312** | +208 |

**8개 case 전부 `stores - loads == 출력 타일 수`를 정확히 만족한다.** 각 타일은 교체될 때마다
쓰이고 돌아올 때마다 읽히며, 완성될 때 한 번 더 쓰인다. T12를 부등식에서 **등식**으로 조였고
known-open 항목은 삭제했다.

### 게이트 두 곳도 결함을 고정하고 있었다

* **T4가 output을 검사하지 않았다.** docstring은 "per type"인데 튜플에 `Output data`가 없었다 —
  하필 write-back이 빠지면 on-chip 양이 off-chip보다 **적어질 수 있는** 유일한 datatype이. 추가했다.
* **`array_buffer`의 `READ_ACCESSES["Output data"] = 0`.** array가 출력을 temporal buffer에 쓰기만
  하고 한 번도 읽지 않는데 4096 B가 전부 DRAM에 도착하는 상태를 고정하고 있었다. 물리적으로 불가능하다.
  이제 타일 하나 읽어내는 4096 access이며, "partial sum을 중간에 읽어오지 않는다"는 원래 의도는
  **정확히 한 타일**이라는 더 강한 형태로 남겼다.

### 회귀

23개 gate 전부 통과(memory_concurrency 포함), unit test 통과,
`check_timing --check-baseline --check-traffic` exit 0. Gemmini RTL 4.40%/7.86%,
Eyeriss latency 4.26%/6.39%, SCALE-Sim PASSED — **셋 다 불변**. GLB band는 5/5 inside를 유지하며
상한만 이동했다(conv3 73.1 -> 73.6). baseline 재생성 완료.
