# Timing Simulation Energy 재평가 및 권장 수정 순서 (2026-08-20)

## 1. 평가 범위

- 대상: 2026-08-20 현재 timing simulation의 energy 및 average-power 경로
- 기준 문서:
  - [TIMING_SIMULATION_ENERGY_REEVALUATION_2026-08-19.md](./TIMING_SIMULATION_ENERGY_REEVALUATION_2026-08-19.md)
  - [ENERGY_POWER_ISSUES_2026-08-19.md](./ENERGY_POWER_ISSUES_2026-08-19.md)
- 중점 재검토 항목:
  - absolute pJ 및 average-power 자격 판정
  - component cost completeness와 실제 event의 연결
  - accumulator reload/spill/final-cast 경계
  - 외부 traffic 및 component energy calibration
  - knob, architecture, NoP energy 회귀 검출력

이번 평가는 코드를 수정하기 위한 문서가 아니라, 현재 구현에서 아직 부족한 부분과 다음 수정 순서를 고정하기 위한 재평가다.

## 2. 결론

현재 구현은 **내부 event-count × unit-cost 산술 정합성**과 **knob/architecture coverage** 측면에서 크게 개선됐다. 이전에 미실행 상태였던 Format-IP, DRAM flat row-miss, local-buffer leakage, MAERI adder-tree 비용 경로에는 fixture와 수기 항등식이 추가됐다. CACTI 7과 DRAMsim3를 사용한 외부 SRAM/DRAM reference도 새로 확보됐다.

그러나 아직 **absolute layer energy 또는 average power가 완전히 교정됐다고 판정할 수는 없다.** 가장 중요한 이유는 다음과 같다.

1. 알려진 미교정 비용이 남아 있어도 report가 total energy를 absolute pJ로, average power를 mW로 출력한다.
2. 실제 event가 발생했지만 관련 unit cost가 빠진 상태를 component schema가 UNDERCOUNT로 판정하지 못한다.
3. Eyeriss 외부 traffic 오차가 남아 있어 SRAM/DRAM 단가가 정확해도 총 dynamic energy의 event count가 실제 칩과 다르다.
4. retained output의 accumulator reload 경계와 2D mesh multicast energy에는 검증 공백이 남아 있다.

따라서 현재 상태는 다음처럼 분류하는 것이 적절하다.

| 범위 | 현재 판정 |
| --- | --- |
| 내부 energy 합계 및 선형성 | 높음 |
| precision 및 accumulator width 반영 | 높음 |
| knob/architecture 실행 coverage | 상당히 개선 |
| SRAM/DRAM unit-cost provenance | 부분 완료 |
| active-event cost completeness | 미완료 |
| absolute layer energy | 미완료 |
| average power의 절대적 의미 | 미완료 |
| peak/system power | 의도적으로 미지원 |

## 3. 이번 버전에서 해결된 부분

### 3.1 이전의 7개 미실행 knob

다음 경로에는 이제 전용 config와 회귀 검사가 존재한다.

- `row_miss_cycle`
- `format_payload_energy`
- `format_metadata_cycle`
- `format_metadata_energy`
- `lb_static_energy`
- `adder_cycle`
- `adder_energy`

관련 checker는 payload/metadata energy-time 분리, DRAM flat activation latency, PE+LB leakage 합, MAERI의 `N-1` reduction additions와 tree-depth latency를 수기 항등식으로 고정한다.

### 3.2 Architecture coverage

TPU, TPUv3, Simba, MAERI, Eyeriss v2, FSD, Gemmini에 synthetic GEMM mapping이 추가됐다. 독립 실행한 architecture checker에서 다음 항목이 통과했다.

- 모든 architecture의 MAC count = `64 × 64 × 64 = 262,144`
- layer total = dynamic + static = 6개 component subtotal 합
- `0 < compute schedule <= critical path`
- DRAM traffic이 distinct tensor volume 이상
- mapping chip factor와 chip utilization 일치
- MAERI adder energy = `4096 outputs × 63 additions × 0.6`

### 3.3 외부 unit-cost reference

- on-chip SRAM: CACTI 7, commit과 geometry 및 node 기록
- off-chip DRAM device: DRAMsim3와 JEDEC IDD parameter에서 access/activation energy 재유도
- reference CSV의 값이 config와 표류하지 않는지 자동 검사

이는 이전의 synthetic unit cost보다 한 단계 높은 provenance를 제공한다. 다만 전체 accelerator energy calibration과는 구분해야 한다.

## 4. 남은 이슈

### E20-1. 미교정 결과가 absolute pJ와 mW 자격을 얻음

**우선순위: P0**

현재 `energy_units_t::is_absolute()`는 다음 두 조건만 확인한다.

```text
energy_unit == pJ
energy_reference가 비어 있지 않음
```

Power gate는 여기에 현재 precision 형식의 `mac_energy_*` key가 존재하는지만 추가로 확인한다. 해당 값이 실제로 외부 calibration된 값인지, 다른 active component의 비용이 빠졌는지는 확인하지 않는다.

`gemmini_cacti22`가 이 문제를 직접 보여 준다. config와 provenance 문서는 MAC compute energy가 외부 calibration되지 않았다고 명시하지만, report는 다음을 동시에 출력한다.

```text
MAC energy basis      : mac_energy_int8_int8_fp32
Energy unit           : pJ (absolute; totals and power are meaningful)
Layer total energy    : 123532954.56 pJ
Average total power   : 418.685 mW
```

같은 report 안에는 다음 미교정 또는 미과금 event도 존재한다.

```text
Layer-setup energy    : 0.00 pJ ... [unit cost UNCALIBRATED]
Accumulator reload    : 1048560 bytes
Accumulator spill     : 1048576 bytes
Accumulator energy    : 0.00 pJ
Output cast           : 4096 bytes
Output cast energy    : 0.00 pJ
```

즉 현재의 `418.685 mW`는 산술적으로는 `energy/time`과 일치하지만, 분자에 알려진 미교정·미과금 항이 있으므로 absolute power로 사용할 수 없다.

#### 필요한 수정

1. config 전체의 free-text `energy_reference` 하나가 아니라 component 또는 cost-key별 calibration 상태를 관리한다.
2. `DECLARED`, `MODELED_ZERO`, `CALIBRATED`, `UNPRICED_ACTIVE`를 구분한다.
3. active event에 `UNPRICED_ACTIVE`가 하나라도 있으면 layer/network total에 `ESTIMATED`를 표시한다.
4. 다음 조건을 모두 만족할 때만 average power와 EDP/ED2P를 출력한다.

```text
single clock domain
AND every active component has an absolute calibrated cost
AND current operand/accumulator precision is calibrated
AND energy_cost_schema().undercounted() == 0
AND no active event is UNPRICED_ACTIVE
```

#### 완료 조건

- `gemmini_cacti22`는 현재 상태에서 wattage를 거부해야 한다.
- provenance에 MAC, accumulator, cast, control, link calibration이 추가된 완전한 fixture만 wattage를 받아야 한다.
- `pJ + reference + precision MAC key`가 있어도 DRAM 또는 NoP cost가 빠진 negative fixture는 wattage를 받지 않아야 한다.

### E20-2. Feature-dependent cost completeness 부재

**우선순위: P0**

현재 component requirement schema는 base cost 일부만 required로 두고 다음을 optional로 취급한다.

- `accumulator_spill_energy`
- `output_cast_energy`
- `layer_setup_energy`
- `weight_fold_fill_energy`
- `row_miss_energy`

또한 다음 key는 energy-key schema에는 추가됐지만 component completeness requirement에는 포함되지 않는다.

- `format_payload_energy`
- `format_metadata_energy`
- `lb_static_energy`
- `adder_energy`

그 결과 feature가 실행됐지만 단가가 없는 상태와 feature가 비활성인 상태를 구분할 수 없다. 예를 들어 `gemmini_cacti22`에서는 setup, accumulator spill/reload, final cast가 실제 발생하지만 `Uncosted components` 경고 없이 component가 costed 상태로 남는다.

#### 필요한 수정

Cost completeness를 config 선언만으로 결정하지 말고 실행된 event와 결합한다.

| feature/event | active 조건 | active일 때 필요한 cost |
| --- | --- | --- |
| layer setup | setup event count > 0 | `layer_setup_energy` |
| weight fold | fold event count > 0 | `weight_fold_fill_energy` |
| accumulator reload/spill | reload/spill bytes > 0 | `accumulator_spill_energy` 또는 분리된 read/write cost |
| final cast | cast bytes > 0 | `output_cast_energy` |
| row activation | activation count > 0 | `row_miss_energy` |
| Format-IP payload/metadata | 해당 transaction > 0 | payload/metadata energy cost |
| MAERI reduction | additions > 0 | `adder_energy` |
| LB leakage | physical LB 존재 및 leakage scope 포함 | `lb_static_energy` 또는 명시적인 modeled zero |

비용을 의도적으로 0으로 모델링하려면 key를 명시적으로 `0`으로 선언하도록 해야 한다. key 부재와 declared zero를 같은 상태로 취급하면 안 된다.

#### 완료 조건

- active event의 cost를 하나씩 제거하는 negative fixture가 모두 `UNPRICED_ACTIVE` 또는 `UNDERCOUNT`를 검출한다.
- inactive feature의 cost 부재는 undercount로 세지 않는다.
- explicit zero는 `MODELED_ZERO`로 남되 report가 해당 사실을 표시한다.

### E20-3. 외부 traffic 오차가 총 energy에 전달됨

**우선순위: P1**

외부 unit cost가 정확해도 event count가 실제 칩과 다르면 total dynamic energy는 정확하지 않다. 현재 Eyeriss 외부 비교는 다음 상태다.

| 항목 | 현재 값 |
| --- | ---: |
| DRAM traffic MAPE | 22.95% |
| DRAM 최대 오차 | 50.00% |
| GLB reference band 내부 | 5개 CONV 중 1개 |

원인은 기존 분석과 동일하다.

- batch filter reuse 미모델링으로 off-chip refetch 과대계상
- psum GLB spill/reload traffic source 부재
- report 내부의 traffic 항등식 통과는 모델 내부 일관성만 보장하며 silicon traffic 정확도를 보장하지 않음

#### 필요한 수정

1. batch/filter retention 범위를 mapping과 buffer capacity에 연결한다.
2. partial sum의 GLB spill/reload source를 명시적인 event로 추가한다.
3. 추가 event가 timing과 energy 양쪽에서 동일한 byte/transaction count를 사용하도록 한다.
4. Eyeriss DRAM/GLB 15% gate를 informational 상태에서 실제 pass/fail gate로 승격한다.

#### 완료 조건

- Eyeriss 지원 CONV layer의 DRAM 및 GLB traffic 오차가 합의된 한계 안에 든다.
- SRAM/DRAM energy 변화가 traffic count 변화와 정확히 선형으로 연결된다.
- 추가된 psum traffic이 accumulator local event와 중복 과금되지 않는다.

### E20-4. Retained-output accumulator reload 경계

**우선순위: P1**

현재 output reload 경로는 `skip_transfer[OUTPUT]` 확인 전에 `account_accumulator_reload()`를 호출한다.

```cpp
account_accumulator_reload(tile_size_mac[OUTPUT]);
if(!skip_transfer[OUTPUT]) {
    // 실제 LB -> MAC transfer
}
```

따라서 output이 MAC에 retained되어 실제 LB→MAC transfer가 생략되는 mapping에서도 reload byte/energy가 과금될 가능성이 있다. 현재 accumulator fixture는 일반 reload 경로의 precision 비율과 절대 cast count를 잘 검증하지만, retained-output mapping은 실행하지 않는다.

새 hierarchy checker도 architecture별 output request counter의 의미가 다르며 단순 conservation law를 아직 세울 수 없다고 `OPEN`으로 기록한다.

#### 필요한 수정

1. reload accounting을 실제 transfer branch 안으로 이동하거나, retained state와 physical reload event를 별도 상태로 표현한다.
2. spill도 실제 eviction boundary에서만 과금되는지 동일하게 점검한다.
3. output residency를 report의 counter와 직접 연결한다.

#### 완료 조건

다음 세 fixture가 필요하다.

- fresh accumulator: create > 0, reload = 0
- spilled accumulator: reload bytes = 실제 LB→MAC accumulator bytes
- retained accumulator: 실제 transfer = 0, reload bytes/energy = 0

그리고 다음 항등식을 고정해야 한다.

```text
reload bytes == physical accumulator read-back elements × accumulator storage bytes
spill bytes  == physical accumulator eviction elements × accumulator storage bytes
```

### E20-5. External calibration의 남은 축과 shipped config migration

**우선순위: P1**

현재 외부 calibration 상태는 다음과 같다.

| 축 | 상태 |
| --- | --- |
| on-chip SRAM | CACTI 7 reference 확보 |
| DRAM device access/activation | DRAMsim3/JEDEC reference 확보 |
| MAC compute | 미교정 |
| accumulator/reduction/cast logic | 미교정 |
| NoP link | 미교정 |
| off-chip PHY/link | 미교정 |
| setup/control | 미교정 |

또한 `gemmini_cacti22`는 별도 calibration fixture일 뿐이고, 기존 shipped architecture config들의 비용을 교체하지 않는다. Phase 4 checker는 기존 Gemmini config의 geometry-cost 불일치 10건이 계속 존재하는지 확인한다. DRAM도 여러 architecture가 device 종류와 무관하게 같은 8 pJ/byte 값을 사용한다.

`dram_config`는 기본 `DRAMSIM3=0` 빌드에서 runtime energy cost를 자동으로 바꾸지 않는다. 따라서 device 이름과 unit cost가 자동 연동된다고 해석하면 안 된다.

#### 필요한 수정

1. MAC energy는 명시된 node/voltage/frequency/precision의 synthesis 또는 검증 가능한 공개 수치로 교정한다.
2. NoP 및 off-chip link는 pJ/bit 또는 pJ/traversal reference를 확보한다.
3. setup/control/cast/reduction은 별도 component 또는 명시적인 미지원 축으로 분리한다.
4. shipped config를 다음 둘 중 하나로 정리한다.
   - 외부 reference에 맞춰 교정하고 component별 provenance를 선언한다.
   - 그렇지 않으면 `normalized` 또는 `uncalibrated`로 명시하고 wattage를 금지한다.
5. `dram_config`에서 비용을 자동 유도할지, config에 고정 값을 선언할지 계약을 하나로 정한다.

### E20-6. 2D mesh multicast energy coverage

**우선순위: P2**

현재 NoP multicast hand identity는 1×4 mesh를 사용한다. 이 fixture는 row 방향 tree만 실행하므로 2D mesh의 column-height 계산이 잘못돼도 검출하지 못한다.

진행 중 mutation test에서는 `column_height[x] = max(column_height[x], y)`의 `max`를 `min`으로 바꾼 변이가 기존 gate를 통과했다. 이는 최소한 2개 row를 사용하는 multicast fixture가 없다는 검출력 공백을 뜻한다.

#### 필요한 수정

- 2×2와 비정방형 2×3 mesh multicast fixture 추가
- spanning-tree distinct-edge count 수기 계산
- broadcast와 partitioned datatype의 energy/serialization 분리
- 동일 active-chip 수를 가진 1×N과 2D grid의 total traversal 차이 검증

#### 완료 조건

- 2D column-height `max→min` 변이가 NoP energy checker에서 실패한다.
- multicast total energy가 실제 distinct router-to-router edge 수와 일치한다.

## 5. 권장 수정 순서

### Phase 1 — Absolute 자격 판정부터 차단

가장 먼저 잘못된 pJ/mW 주장을 막아야 한다. 이후 component 값을 고치는 동안 불완전한 결과가 absolute로 노출되지 않도록 하는 안전장치다.

1. component/cost별 calibration 상태 추가
2. `undercounted() == 0`을 total/power 자격에 연결
3. `UNPRICED_ACTIVE` 상태 추가
4. incomplete fixture에서 total을 `ESTIMATED`로 표시
5. `gemmini_cacti22`의 wattage를 임시로 거부
6. negative power fixture 추가

### Phase 2 — Active-event cost completeness

1. setup/fold/accumulator/cast/row/format/adder event count와 cost declaration 연결
2. absent, explicit zero, calibrated non-zero 구분
3. component annotation 및 network rollup까지 상태 전파
4. 각 optional feature의 cost 제거 negative test 추가

Phase 1과 Phase 2가 끝나야 `absolute`, `calibrated`, `power meaningful`이라는 report 문구를 신뢰할 수 있다.

### Phase 3 — Accumulator/output 경계

1. retained-output mapping 추가
2. reload accounting을 physical transfer boundary에 연결
3. spill을 physical eviction boundary에 연결
4. output request counter 의미를 architecture 공통 report 계약으로 정리
5. create/reload/spill/cast 절대 count를 함께 검증

### Phase 4 — External traffic 정합성

1. psum GLB spill/reload source 구현
2. batch filter reuse/retention 구현
3. traffic와 energy가 같은 event source를 사용하도록 통합
4. Eyeriss DRAM/GLB external gate 활성화

이 단계는 unit-cost calibration과 별개다. 단가와 event count가 모두 맞아야 absolute dynamic energy가 의미를 갖는다.

### Phase 5 — 남은 component calibration 및 config migration

1. MAC compute
2. accumulator/reduction/cast/control
3. NoP link
4. off-chip PHY/link
5. shipped architecture config migration
6. 완전히 교정된 fixture에서만 average power 재활성화

### Phase 6 — Coverage 강화

1. 2D mesh multicast fixture
2. frequency=0 및 derived-bitwidth validation
3. source/destination/link transaction이 모두 다른 overflow fixture
4. mutation survivor 개별 판정 및 kill test 추가
5. 전체 checker를 직렬 또는 격리된 build/result directory에서 실행

## 6. 최종 완료 체크리스트

다음 조건을 모두 만족하기 전에는 timing simulation energy 완료를 선언하지 않는다.

- [ ] active event에 빠진 energy cost가 하나도 없다.
- [ ] absent cost와 explicit zero가 report에서 구분된다.
- [ ] total energy의 `absolute` 표시는 모든 active component calibration 완료를 요구한다.
- [ ] average power는 incomplete/undercounted config에서 출력되지 않는다.
- [ ] `gemmini_cacti22`의 MAC, accumulator, cast, setup, NoP/link 비용 provenance가 완성된다.
- [ ] retained output이 reload/spill energy를 내지 않는 fixture가 통과한다.
- [ ] Eyeriss DRAM/GLB traffic이 외부 목표 범위 안에 든다.
- [ ] 2D mesh multicast distinct-edge identity가 통과한다.
- [ ] shipped config의 energy unit과 provenance가 실제 비용 출처와 일치한다.
- [ ] peak/system power 미지원 범위가 report에 계속 명시된다.

## 7. 검증 상태에 대한 주의

이번 재평가 시점에 같은 workspace에서 별도의 mutation campaign이 source file을 변이하고 `npusim.sh build npusim`을 반복 실행하고 있었다. 그 과정에서 일괄 checker 실행의 binary와 result가 서로 다른 mutant에 의해 오염되는 현상이 관찰됐다. 예를 들어 energy cost 변경이 timing field를 움직이거나 simulator binary가 순간적으로 실행되지 않는 비정상 결과가 발생했다.

따라서 본 문서의 판정은 다음 근거를 사용했다.

- 안정적으로 완료된 309개 config 정적 validation
- 독립 실행으로 통과한 accumulator 및 architecture checker
- 개별적으로 완료된 energy/schema/format/DRAM/Phase-4/power/precision/reduction/traffic checker 결과
- source와 generated report의 직접 대조
- mutation campaign의 stable backup source와 이미 보고된 survivor 문맥

Mutation campaign이 종료되고 원본 source/build가 복원된 뒤에는 전체 checker를 한 번 더 깨끗한 상태에서 실행해야 한다. 이 재실행은 위 코드·모델 이슈를 해소하는 것은 아니지만, 현재 회귀 suite의 최종 baseline을 확정하는 데 필요하다.

## 8. 최종 판단

현재 버전은 **normalized/relative breakdown 연구와 내부 energy accounting 검증**에는 충분히 유용하다. CACTI/DRAMsim3 reference 도입으로 일부 absolute component cost도 추적 가능해졌다.

그러나 현재 report가 출력하는 absolute layer pJ와 average mW는 아직 완전한 accelerator 또는 chip power로 해석하면 안 된다. 우선 `absolute` 자격 판정을 active-event calibration에 연결하고, 그 다음 accumulator/output 경계와 external traffic을 해결하는 순서가 가장 안전하다.
