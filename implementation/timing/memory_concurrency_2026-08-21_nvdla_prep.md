# 2026-08-21 memory-interface concurrency (outstanding-request) 구현 기록

기준 문서: [../../assessment/NVDLA_RTL_FIDELITY_PLAN_2026-08-20.md](../../assessment/NVDLA_RTL_FIDELITY_PLAN_2026-08-20.md)

## 배경

NVDLA fidelity 계획서 8.1/8.2절은 memory 조건을 plusarg로 명시하고 sensitivity suite에서
쓸어야 할 축으로 다음을 요구한다.

- read latency: 32/64/128/180 cycles
- **max outstanding request: 128** (baseline), sensitivity 대상
- bandwidth utilization: 25%/50%/100%

이 중 read latency와 bandwidth는 기존 `[dram]` unit cost/row-buffer 설정으로 이미 표현 가능하다.
그러나 "max outstanding request"(DBBIF가 동시에 진행 중일 수 있는 요청 수)는 착수 전 코드
감사 결과 **구조적으로 표현 불가능**했다: `scheduler/stats.cc`의 per-tile pipeline timeline
(`utils/interconnect_timing.h`의 `pipeline_timeline_cycles()`)은 이미 임의의 N-deep staging
buffer를 일반적으로 지원하는데도, DRAM↔multi-chip 경계의 depth는 오직
`multi_chip_t::double_buffer` 하나의 bool에서 `? 2 : 1`로만 도출되고 있었다 — 1 또는 2 외의
값을 만들 방법이 없었다. `components/dram.cc:239-250`의 `describe_timing_limits()`도 이를
"request queueing... not modeled"로 명시하고 있었다.

## 구현

기존 메커니즘을 그대로 재사용하고, depth를 도출하는 지점 하나만 확장했다.

- `components/multi_chip.h`/`.cc` : `unsigned max_outstanding_requests` 필드 추가.
  `[multi_chip]` 섹션의 `max_outstanding_requests` 키에서 읽으며 기본값은 0("미설정").
- `scheduler/stats.h` : `timeline_boundary_max_requests[4]` 추가 (경계 0 = DRAM↔multi-chip만
  현재 사용).
- `scheduler/stats.cc` :
  - `update_stats()`에서 `timeline_boundary_max_requests[0] = m_multi_chip->max_outstanding_requests`.
  - `finalize_layer_timeline()`의 depth 계산을
    `timeline_boundary_overlap[b] ? 2 : 1` 에서
    `timeline_boundary_max_requests[b] > 0 ? timeline_boundary_max_requests[b] : (timeline_boundary_overlap[b] ? 2 : 1)`
    로 확장.

0(미설정)일 때는 기존 `double_buffer` 기반 1/2 depth를 그대로 재현하므로, 이 키를 선언하지 않는
기존 config는 전부 이전과 **bit-identical**하다 (아래 회귀 결과 참고). NVDLA config를 작성할 때
`[multi_chip] max_outstanding_requests = 128`을 선언하면 DBBIF의 outstanding-request 한도를
그대로 표현할 수 있다.

## 검증

### 회귀 (기존 config, 변경 없음을 확인)

수정 전/후 동일하게 재실행, 세 gate 모두 숫자까지 bit-identical:

- `validation/bottleneck/check.py` — gemmini 59,010 / gemmini_memory_bound 384,942 /
  gemmini_compute_bound 85,510 / gemmini_2chip 50,166 (critical-path, cycles) — ALL PASSED
- `validation/traffic/check.py` — T1-T11 전부 동일 — ALL PASSED
- `validation/buffering/check.py` — LB-A~E 전부 동일 (LB-E: 59010.0 -> 59889.5) — ALL PASSED

### 신규 기능 (`validation/memory_concurrency/check.py`, MC1-MC5)

`gemmini_memory_bound.cfg`에 한 줄만 다른 두 변형을 추가했다.

| config | max_outstanding_requests | DRAM->Multi-chip depth | DRAM stall | Critical path |
| --- | ---: | ---: | ---: | ---: |
| `gemmini_memory_bound` | (미설정) | 1 | 1,152.0 | 384,942.0 |
| `gemmini_memory_bound_outstanding1` | 1 (명시) | 1 | 1,152.0 | 384,942.0 |
| `gemmini_memory_bound_outstanding` | 128 | 128 | 0.0 | 383,790.0 |

- MC1 명시적 depth=1이 암묵적 기본값과 완전히 동일함을 확인 (override 경로가 별개 코드가
  아님을 증명).
- MC2 depth override가 DRAM↔multi-chip 경계에만 적용되고 다른 세 경계는 불변임을 확인.
- MC3 Computation cycle/Fold fill cycle이 세 실행 모두 bit-identical — 재분류가 concurrency
  지식 하나에서만 왔음을 증명 (기존 BN3/BN4와 동일한 isolation 논증).
- MC4 depth를 올렸을 때 stall/critical path가 절대 나빠지지 않음을 확인.
- MC5 stall 감소분(1,152.0 → 0.0)과 critical path 감소분(384,942.0 → 383,790.0)이 정확히
  일치함을 확인 — pipeline timeline 자체의 back-pressure 산술을 end-to-end로 검증.

ALL MEMORY-CONCURRENCY CHECKS PASSED.

## 범위와 한계

- 여전히 **tile 단위** depth cap이다. NVDLA RTL의 실제 per-request address 기반 큐잉/bank
  conflict는 여전히 모델링하지 않는다 (`dram_t::describe_timing_limits()`의 disclaimer는
  그대로 유효 — 고치지 않았다. 이 fix는 그 disclaimer가 말하는 "request queueing"의 더 거친
  근사 하나를 추가했을 뿐, 그 disclaimer 자체를 무효화하지 않는다).
- 현재 경계 0(DRAM↔multi-chip)에만 배선했다. 다른 경계(예: multi-chip↔GLB)에 같은 개념이
  필요해지면 `timeline_boundary_max_requests[1..3]`과 해당 컴포넌트에 동일한 패턴을 반복하면
  된다 — 구조는 이미 일반적이다.
- 이 fix는 계획서 8.2절의 sensitivity suite를 **표현 가능**하게 만들 뿐이다. 실제 NVDLA RTL
  대비 outstanding-request 축의 정합성(계획서 8.3절 gate)은 Phase 0/1 golden data 없이는
  아직 검증되지 않았다.
