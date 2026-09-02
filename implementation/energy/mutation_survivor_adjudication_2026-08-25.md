# 2026-08-25 mutation survivor 9건 판정

`validation/mutation_partial.json`에 13 mutant(4 killed / 9 survived)까지만 돌고 중단된 캠페인의
**survivor 9건을 개별 판정**했다. 결과 원본: `validation/mutation_survivors_2026-08-25.json`.

## 방법

기억이나 추론이 아니라 **현재 스위트에 다시 돌렸다**. 캠페인 이후 소스가 많이 바뀌었기 때문에
line 번호가 아니라 **소스 스니펫으로 anchor**를 잡았고(2건은 중복이라 line 지정),
각 mutant마다 build → `unittest` → 23 gate → `check_timing --check-baseline --check-traffic`을
순서대로 돌려 **처음 실패한 검사를 killer로 기록**했다.

결과: **1 KILLED, 8 SURVIVED.**

이어서 살아남은 8건이 *왜* 살아남았는지를 갈랐다. 핵심 질문은 "게이트가 안 보는 것인가, 코드가 안
도는 것인가"이고, 이는 **리포트 diff**로 바로 갈린다 — 변이 후 리포트가 비트 단위로 동일하면 그 줄은
실행되지 않은 것이다.

## 판정

| mutant | 판정 | 근거 |
| --- | --- | --- |
| `interconnect/tallest_row` max→min | 등가 (죽은 코드) | `(void)tallest_row;` — 대입 후 한 번도 읽지 않는다 |
| `interconnect/pipeline_max` max→min | **unittest가 죽임** | `validation_test.cc`가 10/2/5/18/3으로 고정 |
| `interconnect/freq_guard` <=→< | 등가 (UB 실현값 일치) | 가드를 빼면 `8.0*bw/0.0 = inf`의 UB 캐스트가 마침 0 |
| `interconnect/pow2_check` -1→-2 | **증명 가능한 등가** | `w>=8`에서 두 술어가 항상 일치(0~2^20 전수 확인) |
| `pe/transfer_cycle` +=→-= | **진짜 구멍 → 수정** | 아래 참조 |
| `pe/lb_weight_read_energy` | 도달 불가 | 131회 실행에서 hit 0 |
| `pe/output_transfer_energy` | 도달 불가 | 131회 실행에서 hit 0 |
| `pe/array_write_back_cycle` | 죽은 코드 (주석) | `/* ... */` 블록 안 — harness가 주석을 추적하지 않았다 |
| `pe/lb_output_read_cycle` | 도달 불가 | 131회 실행에서 hit 0 |

**"도달 불가"의 근거**: 세 지점에 hit counter를 심고 `configs/mappings` 전체를 훑어
**131개 accelerator × network × mapping 조합**을 실행했다. 세 지점 모두 총 hit **0**이다.
`data_transfer_to_mac`/`flush_data`의 구형 address-walking 경로이며, analytical 경로
(`account_descriptor_dense_mac_transfer`)로 대체된 뒤 남은 잔재다. 컴파일은 되지만
(`#ifdef FUNCTIONAL`의 `#else` 쪽) 어떤 shipped config로도 실행되지 않는다.

## 유일한 진짜 구멍 — `transfer_cycle[type] += ... → -=`

이 변이만 리포트를 **실제로 바꾼다**. 부호가 뒤집히면 PE별 MAC↔LB transfer cycle이 음수가 되는데,
`stats_t`가 **0에서 시작하는 `std::max`**로 집계하므로 리포트에는 깔끔한 `0.0`이 찍힌다.
eyeriss와 gemmini 양쪽에서 그렇다.

**아무도 눈치채지 못한 이유가 네 겹이다.**

* `KN4`(dead-knob sweep)는 통과한다 — knob이 여전히 다른 숫자를 **움직이기** 때문이다. liveness에는
  부호 개념이 없다.
* `KN7`(latency monotonicity)은 **critical path만** 본다. 이 축은 critical path를 정하지 않는다.
* `T7`은 축들을 `max()`로 합쳐 stage busy를 만든다. PE stage는 다른 축이 지배하므로 link 축이 0으로
  주저앉아도 busy는 그대로다.
* 그리고 **어떤 게이트도 `Cycle (MAC-Local buffer)` 블록을 읽지 않았다.**

즉 이 축은 **출력은 되지만 검증되지 않는** 상태였다.

### 수정 — `validation/pe_transfer`

| check | 내용 |
| --- | --- |
| **PT1** | serialized MAC↔LB transaction이 있고 `transfer_cycle_pe > 0`이면 transfer cycle은 **양수**여야 한다. 변이가 위반하는 지점이다 — 일이 일어났고 가격도 선언됐는데 비용이 0이다. |
| **PT2** | `transfer_cycle_pe`를 2배로 하면 transfer cycle이 정확히 2배가 되고 transaction 수는 안 움직인다 |
| **PT3** | PE stage의 **link busy 축 == 세 datatype transfer cycle의 합** (gemmini·eyeriss 두 config에서 정확히 성립) |
| **PT4** | `transfer_energy_pe`는 **68개 shipped config 전부 0**이라 그 곱셈기가 한 번도 실행되지 않는다 — `noc_energy`(17배 틀린 채 발각 안 됨)와 같은 형태. 가격을 매겨 축을 정확히 검사하고 cycle은 안 움직임을 확인한다 |

**검증**: 변이를 다시 적용하니 PT1이 세 datatype 모두에서 실패한다 —
`65536 MAC<->LB transactions at transfer_cycle_pe=1.0 report 0.0 transfer cycles`.

## 조치

* `utils/interconnect_timing.cc`의 죽은 `tallest_row` 변수 **삭제** (변이 지점 자체가 사라진다).
* `validation/pe_transfer` **신규 게이트** 추가 (총 24개).
* 판정 근거를 `validation/mutation_survivors_2026-08-25.json`에 mutant별 사유와 함께 기록.

## 캠페인 harness에 남는 개선점

이번 판정이 harness 자체의 결함 두 가지를 드러냈다.

1. **unit test를 kill set에 넣지 않았다.** 9건 중 1건은 unit test가 6초 만에 죽인다.
2. **`/* */` 블록 주석을 추적하지 않는다.** `#ifdef` region은 추적하도록 고쳤지만 주석은 아니어서,
   컴파일러가 보지도 않는 텍스트를 변이했다.

두 가지를 고치면 60 mutant 완주 시 survivor 수가 의미 있게 줄어든다. 다만 이번 9건에서 보듯
**남는 survivor의 다수는 "게이트 구멍"이 아니라 "도달 불가 코드"**일 가능성이 높고, 그건 테스트를
추가할 대상이 아니라 **삭제를 검토할 대상**이다.

---

# 정정 및 죽은 코드 삭제 (같은 날, 나중)

## 정정 — "도달 불가 3건"은 2건이 틀렸다

위 표의 "도달 불가" 3건을 삭제하려고 구조를 확인하다 **제 규정이 틀렸다는 것**을 발견했다.
`pe_t::data_transfer_to_mac`과 `pe_t::flush_data`는 둘 다 이 모양이다:

```cpp
#ifndef FUNCTIONAL
    ... analytical 회계 ...
    return;            // pe.cc 859행 / 2068행
#endif
    ... 수백 줄 ...     // FUNCTIONAL 빌드가 실행하는 구현
```

즉 early return 뒤의 코드는 **`FUNCTIONAL=1` 빌드의 구현체**다(`npusim.sh:60,81`이 지원하는 옵션).
131회 sweep은 기본 빌드만 돌렸으므로 그 "hit 0"은 early return의 당연한 귀결이지 죽은 코드의
증거가 아니었다. `pe/output_transfer_energy`와 `pe/lb_output_read_cycle` 두 건은 **삭제하면 안 되는
살아있는 코드**다. 지웠다면 functional simulation 경로를 통째로 없앨 뻔했다.

(같은 맥락에서, 처음에 `-DFUNCTIONAL` 컴파일이 통과한다고 본 것도 틀렸다 — `$?`가 파이프라인의
값이라 0으로 보였을 뿐, 실제로는 nebula include 경로가 빠져 전후 모두 실패하고 있었다.)

## 그렇다면 무엇이 정말 죽었나

정확한 규칙은 이것이다. **early return **뒤**에 있는 `#ifdef FUNCTIONAL` 그룹의 `#else` 분기**는
두 빌드 어디에서도 실행되지 않는다:

* FUNCTIONAL 빌드 — `#else`는 애초에 컴파일되지 않는다.
* analytical 빌드 — 컴파일되지만 `return;` 뒤라 도달하지 않는다.

이 규칙으로 전 소스를 훑어 **6개 영역 633줄**을 찾았다.

| 파일 | 함수 | early return | 죽은 `#else` | 줄 |
| --- | --- | ---: | ---: | ---: |
| `components/pe.cc` | `data_transfer_to_mac` | 859 | 1229–1330 | 102 |
| `components/pe.cc` | `data_transfer_to_mac` | 859 | 1728–1830 | 103 |
| `components/adder_tree.cc` | `data_transfer` | 262 | 697–803 | 107 |
| `components/adder_tree.cc` | `data_transfer` | 262 | 1251–1356 | 106 |
| `components/spatial_arch.cc` | `data_transfer` | 241 | 680–788 | 109 |
| `components/spatial_arch.cc` | `data_transfer` | 241 | 1233–1338 | 106 |

여기에 `pe.cc`의 `/* ... */` 주석 블록 13줄(`pe/array_write_back_cycle` 변이 지점)을 더해
**총 646줄을 삭제**했다. survivor 중 `pe/lb_weight_read_energy`(1775행)가 두 번째 영역 안에 있었으므로
그 변이 지점도 함께 사라졌다.

## 검증

* **analytical 리포트가 비트 단위로 동일**하다 — eyeriss conv3, gemm512, gemmini_all_knobs 세 리포트를
  삭제 전후로 `diff`해 차이 0.
* **FUNCTIONAL 빌드가 여전히 컴파일된다** — 세 파일을 올바른 include 경로로
  `g++ -fsyntax-only -DFUNCTIONAL`에 통과시켰다.
* 24개 gate 전부 통과, unit test 통과, `check_timing --check-baseline --check-traffic` exit 0,
  SCALE-Sim PASSED.

## 남는 것

`pe/output_transfer_energy`, `pe/lb_output_read_cycle` 두 지점은 FUNCTIONAL 빌드 코드로 **남긴다**.
이들을 mutation으로 검증하려면 캠페인이 `FUNCTIONAL=1` 빌드도 돌려야 하는데, 현재 회귀 suite에는
functional 경로를 검사하는 게이트가 하나도 없다. 이는 **테스트 구멍이 아니라 검증 범위의 공백**이며,
별도 항목으로 둔다.
