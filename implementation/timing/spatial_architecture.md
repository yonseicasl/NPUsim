# Spatial architecture timing implementation

## 상태: 완료 (명시된 timing abstraction 범위)

BUS/store-and-forward/crossbar는 serialized-unicast, MESH는 ingress `(0,0)` 기반 Manhattan routed-unicast로 구현했다. multicast는 destination별 replication이며 contention은 destination transaction 직렬화로 정의한다.

## 구현된 변경

### Topology, route 및 transaction 정책

- `noc = bus`, `store_and_forward`, `crossbar`, `mesh`를 지원한다.
- MESH latency는 active array의 최대 Manhattan hop, energy는 active destination의 평균 Manhattan hop을 사용한다.
- BUS/store-and-forward/crossbar는 1-hop serialized fabric으로 계산한다.
- dense transfer cycle은 `num_access_pe` destination transaction을 직렬화한다. fanout/multicast는 destination별 replication energy로 계산한다.

### Tile state 및 overlap

- tile update마다 `read_tile_granular_pe_input/weight/output`를 reset한다.
- input/weight/output pipeline overlap에서 0 access를 별도 처리한다.
- active PE shape를 physical height/width와 비교해 검증한다.

## 검증

`validation_test.cc`의 spatial interconnect micro-test는 다음을 확인한다.

- BUS/store-and-forward/crossbar/MESH topology 지원
- 2×3 MESH의 latency multiplier 3 및 energy multiplier 1.5
- 4×8 physical array의 유효/무효 active shape
- 0, 1, 4 transaction pipeline cycle

```bash
./unittest/run_validation.sh
```

`spatial_arch.cc`는 `-Wall -Wextra -Werror`로 독립 컴파일한다.

## 모델 경계

router queue, virtual channel, adaptive route 및 cycle-accurate back-pressure는 이 abstraction의 대상이 아니다. 이를 요구하는 hardware 비교에는 별도 topology model이 필요하다.
