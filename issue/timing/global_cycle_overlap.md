# Global cycle / compute-memory overlap: timing simulation 이슈

> **재판정 (2026-09-02): 종결.** 이 문서의 "후속 정밀화 범위"였던 per-tile event
> timeline까지 구현·검증 완료됨 — `pipeline_timeline_cycles()`(L5)가 stage별 per-tile
> 비용 벡터와 boundary staging depth(1=단일버퍼, 2=이중, N=queue; MC2
> `max_outstanding_requests` 오버라이드 포함)로 back-pressure stall attribution까지
> 계산하고, off-chip 재인출률(L5)·GLB bypass(L6)·SFU post-processing 6번째
> stage(plan_sfu Phase 5, queue_depth 계약)가 같은 timeline 위에 있다.
> `unittest/validation_test.cc`(hand-calc stall/기간 검증)와
> `stats_timeline_test.cc`(SFU streaming 4케이스)가 계약을 고정한다. 아래 본문은
> 2026-08-15 stage-granular 단계까지의 기록이다.

## 판정

현재 timing 모델에는 **전역 clock(공유 timeline)이 없다.** 각 컴포넌트(PE, PE array, GLB, multi-chip, DRAM)가 자신의 `cycle`/`energy` 카운터를 독립적으로 누적하고, `stats_t`가 이를 컴포넌트 내부는 `max`, 네트워크 레이어 간은 `+`로 집계한다. 계층 간 전송과 연산이 **하나의 시간축 위에서 겹치는(overlap) 관계가 표현되지 않으므로**, compute-bound/memory-bound 판정, double-buffering 이득, 파이프라인 latency hiding을 정량화할 수 없다.

- 심각도: ~~P1 (구조적)~~ → **stage-granular analytical timeline 구현됨 (2026-08-15)**. per-tile event timeline은 후속 정밀화 범위.
- 성격: 결함 수정이 아니라 timing 엔진의 실행 모델 확장

## 구현 결과 (2026-08-15)

`stats.cc::update_stats`가 5개 계층(DRAM/multi-chip/GLB/PE-array/PE)의 busy time(`modeled_elapsed_cycles`)을 공유 timeline에 배치한다. 경계 overlap 규칙은 access-cycle 회계와 동일한 **destination double-buffer 관례**: `DRAM|MC`는 `multi_chip.double_buffer`, `MC|GLB`는 GLB `double_buffer`, `GLB|PE-array|PE`는 동일 스트림이라 항상 overlap. **layer critical-path latency = overlap 구간은 max, 직렬 경계는 sum**의 segment-combine. 이 latency가 leakage window와 `mac_available_cycle`을 구동하며, 결과 파일에 "Layer timeline" 섹션(latency + 계층별 busy ratio + bottleneck)으로 출력된다. repetition 스케일·network 합산 포함.

**검증 (eyeriss alexnet layer_0):**
- 산식 정합: latency 674,744,400 = DRAM 253.1M + MC 168.7M + max(GLB, PE-array, PE) 253.0M (수작업 일치)
- double-buffer 반응: `double_buffer=1` → 674.7M → **253.1M** (= max로 붕괴)
- bound 판정: 기본 config는 DRAM(37.5%) bottleneck·PE 3.4%; `computation_cycle=100`이면 **PE(47.6%)로 전환**

## 현재 동작의 근거

- 실행 루프는 `scheduler/npu.cc:312`의 `while(!is_idle())`로, tile 이동 상태 기계를 반복 호출할 뿐 **증가하는 global cycle 변수가 없다.** 각 `data_transfer()`가 자기 컴포넌트의 `transfer_cycle[type]` 등에 비용을 더한다.
- 집계는 `scheduler/stats.cc`에서 이뤄진다. 컴포넌트 내부(PE들 간, GLB들 간)는 `std::max`(예: `stats.cc:431,472,517,539`), 레이어 간은 `+=`(예: `stats.cc:751~`)다. **PE↔GLB↔DRAM 계층을 하나의 critical-path latency로 합치는 지점이 없다.**
- 각 계층의 cycle은 결과 파일에 개별 항목으로만 출력된다(computation/access/transfer cycle이 컴포넌트별로 따로). "총 실행 latency" = max(compute, memory) 형태의 단일 수치가 산출되지 않는다.
- 그 결과 `pe_elapsed_cycles`(정적 에너지 창)도 PE 내부 축의 `max`만 사용하고 memory stall 시간을 포함하지 못한다([static_energy.md](static_energy.md), `stats.cc:400-421`). 이는 global timeline 부재의 직접적 파생 결함이다.

## 왜 문제인가

| 표현하지 못하는 것 | 결과 |
|---|---|
| compute와 memory 전송의 시간축 상 overlap | latency가 "각 컴포넌트 cycle의 합"이나 "특정 계층의 max"로 왜곡됨 |
| double-buffering / prefetch에 의한 latency hiding | 버퍼링 유무가 총 latency에 반영되지 않음 |
| bound 판정(compute-bound vs memory-bound) | 병목 계층을 식별할 근거가 없음 |
| back-pressure / stall 전파 | 생산자-소비자 속도차가 무시됨 |
| 정적 에너지 시간창 | leakage가 실제 wall-clock이 아닌 PE 부분 창에 곱해짐 |

## 수정 방향 (제안)

목표는 **컴포넌트가 소비/생산하는 시간을 공유 timeline에 배치**해, 총 latency를 "합"이 아니라 계층 파이프라인의 **critical path**로 얻는 것이다. 단계적으로:

1. **Timeline 표현 도입.** 각 계층 경계(DRAM→multi-chip→GLB→PE-array→PE)에 대해 producer `busy interval`과 consumer `busy interval`을 tile 단위로 기록한다. 최소 구현은 계층별 `available_cycle`(다음 tile을 받을 수 있는 시점)과 `ready_cycle`(데이터가 준비된 시점) 두 스칼라를 유지하는 것으로 충분하다.
2. **Overlap 규칙.** consumer의 시작 시점은 `max(consumer_available, producer_ready)`. double-buffer가 있으면 producer가 tile n+1을 채우는 동안 consumer가 tile n을 소비하므로 두 interval이 겹친다. 버퍼 깊이(`double_buffer`, `exist_temporal_buffer`)를 overlap 허용량으로 사용한다.
3. **Layer latency = critical path.** 레이어 종료 cycle을 파이프라인 정상상태의 병목 stage 처리율과 fill/drain으로 계산한다. 현재의 `max`/`sum` 집계를 이 timeline 종료 cycle로 대체한다.
4. **Bound 리포팅.** 계층별 busy 비율(`busy_cycle / layer_cycle`)을 출력해 병목 계층을 표시한다.
5. **정적 에너지 재연결.** [static_energy.md](static_energy.md)의 elapsed-cycle 창을 새 layer critical-path latency로 교체한다.

기존 analytical 비용(`datatype_transfer_timing()`, `pipelined_transfer_cycles()`, `spatial_noc_cost()`)은 stage당 latency/throughput 입력으로 재사용한다. 즉 **비용 계산식은 유지하고, 그 비용들을 배치하는 timeline만 추가**한다.

## 상호작용하는 기존 이슈

- [static_energy.md](static_energy.md) — leakage 시간창은 이 timeline이 생기면 정확해진다.
- [systolic_array.md](systolic_array.md) — wavefront fill/drain은 PE-array 내부 timeline을 요구하므로 같은 프레임에서 다뤄진다.
- [global_buffer.md](global_buffer.md), [multi_chip.md](multi_chip.md) — contention/back-pressure는 timeline 위에서만 정의된다.

## 완료 조건

- [~] DRAM↔multi-chip↔GLB↔PE-array↔PE 경계가 공유 timeline에 배치된다 — **stage-granular busy interval + 경계 overlap flag로 구현**; per-tile producer/consumer interval은 후속 정밀화. **(2026-08-15 종결 판단: per-tile event timeline은 `while(!is_idle())` tile 상태기계에 clock을 관통시키는 엔진 재설계 과제로, 별도 설계 리뷰 후 착수 — 현 stage-granular 모델이 latency·bound·overlap의 1차 계약을 제공하며 모든 비용식(descriptor·hop·row-buffer·bank)은 재사용 가능)**
- [x] double-buffer 유무가 총 layer latency를 변화시킨다(overlap 반영). (2026-08-15)
- [x] 레이어 결과에 단일 critical-path latency와 계층별 busy 비율(병목 표시)이 출력된다. (2026-08-15)
- [x] compute-bound와 memory-bound 합성 config가 서로 다른 병목 계층으로 보고된다. (2026-08-15)
- [x] 정적 에너지 시간창이 새 critical-path latency를 사용한다(`mac_available_cycle`도 동일). (2026-08-15)
