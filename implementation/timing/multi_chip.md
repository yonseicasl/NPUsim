# Multi-chip / NoP timing implementation

## 상태: dense datatype 완료, routed NoP 부분 완료

Multi-chip은 현재 BUS 기반 serialized NoP model을 제공한다. INPUT/WEIGHT/OUTPUT의 GLB↔multi-chip↔DRAM dense 경계는 runtime datatype descriptor transaction을 사용하며, payload/metadata/serialized traffic을 별도로 집계한다.

## 구현된 변경

### Configuration 및 topology 계약

- `height`, `width`가 0이 아니고 곱셈 overflow가 없음을 확인한다.
- frequency, bandwidth, bitwidth의 유효 범위를 검사한다.
- active chip shape가 physical `height × width`를 넘지 않는지 확인한다.
- `nop = bus`만 허용하며 mesh 등 미구현 NoP는 초기화에서 명시적으로 거부한다.

### Transaction 및 timing 안전성

- INPUT/WEIGHT **load**의 chip, GLB, DRAM access 및 link transfer는 descriptor storage transaction으로 계산한다.
- `line_size`는 bits 단위의 8 이상 power-of-two로 검증하고, link `bitwidth`도 bits로 통일한다.
- shared/separate temporal buffer capacity는 runtime datatype storage byte와 overflow-safe 합으로 검사한다.
- NoP 및 DRAM 경계의 payload/MXFP metadata/serialized transaction을 결과에 보고한다.
- GLB↔chip NoP transfer energy는 cycle과 동일하게 `num_access_multi_chip`을 사용한다.
- input/weight/output overlap의 0 transaction을 별도 처리해 unsigned underflow를 막는다.

### Static energy

- multi-chip transfer callback이 inactive chip에만 static energy를 더하던 경로를 제거했다.
- GLB static energy는 stats에서 physical GLB 전체에 modeled layer elapsed cycle 기준으로 한 번 집계한다.

> ⚠️ **정정 (2026-08-15):** OUTPUT **write-back**의 NoP link 비용(`multi_chip->write_back_cycle`)은 주석 처리되어 0으로 집계된다. GLB static energy의 "modeled layer elapsed cycle"은 실제로는 **PE-array elapsed**이며 NoP/DRAM 시간을 포함하지 않는다. 또한 패키지 temporal buffer 자체의 static energy는 미집계다. 이슈: [issue/timing/multi_chip.md](../../issue/timing/multi_chip.md), [issue/timing/static_energy.md](../../issue/timing/static_energy.md).

## 검증

```bash
./unittest/run_validation.sh
```

`multi_chip.cc`를 `-Wall -Wextra -Werror`로 독립 컴파일한다. interconnect micro-test는 BUS 지원, 미구현 NoP 거부, active shape와 0/1/N pipeline 계약을 검증한다.

## 남은 완료 조건

- route/hop/router/link latency·energy, contention, multicast/fanout, injection/ejection 및 back-pressure를 갖는 routed NoP model
