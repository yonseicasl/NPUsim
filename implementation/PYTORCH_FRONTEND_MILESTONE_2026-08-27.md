# PyTorch Frontend 1차 구현 결과

작성일: 2026-08-27

## 결론

`plan/plan_pytorch_frontend_migration.md`의 Phase 1~3 중 첫 실행 가능한 vertical slice를 구현했다. 현재 NPUsim은 PyTorch capture artifact를 strict lowering하여 생성한 `npusim.exec.v1` JSON을 C++에서 직접 읽고, 기존 scheduler와 timing/traffic/energy model로 실행할 수 있다.

현재 경로는 다음과 같다.

```text
PyTorch torch.export
  -> npusim.graph.v1
  -> Python strict lowering
  -> npusim.exec.v1
  -> C++ workload_graph_t
  -> transitional Nebula adapter
  -> NPUsim cycle/traffic/energy core
```

PyTorch는 export 시점에만 필요하다. 이미 생성된 capture/executable artifact의 검증·compile·simulation에는 PyTorch가 필요하지 않다.

## 구현 범위

### Capture 및 lowering

- positional/keyword literal argument 보존
- tensor stride, storage offset, contiguous/strided layout 보존
- parameter/buffer qualified name 보존
- `aten.linear.default` 지원
- `aten.conv2d.default` 지원
- non-transposed `aten.convolution.default` 지원
- Linear/Conv2d의 single-consumer ReLU·Leaky ReLU fusion
- last-axis Softmax 지원
- symbolic executable shape와 unsupported operation의 fail-fast
- Capture IR 및 Executable IR content hash 생성·Python 검증

### C++ 실행 경로

- `scheduler/workload_graph.*`의 `npusim.exec.v1` loader 추가
- schema, tensor dtype/byte, topology, operation 상태, geometry/shape 교차 검증
- executable operation ID와 mapping `op_id` 결합
- mixed/duplicate/missing/extra/wrong-kind mapping의 fail-fast
- `models/model run-ir ACCELERATOR EXECUTABLE MAPPING [RESULT_NAME]` 추가
- timing-only 실행에서는 Nebula dataset loading을 수행하지 않음
- functional build에서 tensor artifact 없이 IR을 실행하면 fail-fast
- nonlinear activation/Softmax를 SFU가 없는 accelerator에서 실행하면 fail-fast
- layer/network 결과에 frontend schema와 graph/executable provenance 기록

## 검증 결과

다음 명령이 통과한다.

```bash
unittest/run_pytorch_frontend_validation.sh
./npusim.sh build npusim
models/model run-ir \
  configs/accelerators/gemmini_sfu.cfg \
  frontend/fixtures/linear_relu.exec.json \
  frontend/fixtures/linear_relu.map \
  pytorch_linear_relu
```

`linear_relu` fixture의 실제 end-to-end 결과는 다음과 같다.

- executable operation: `npusim.linear`, fused activation `relu`
- geometry: batch 64, input features 64, output features 64
- compute-schedule latency: 3,518 cycles
- critical-path latency: 59,394 cycles
- SFU busy cycle: 256 cycles
- network total energy: 1,259,882.71 uncalibrated units
- energy scope: 1/1 layer complete

마지막 energy 값은 accelerator config가 `energy_unit`과 calibration provenance를 선언하지 않으므로 절대 pJ가 아니다. 이번 테스트의 목적은 PyTorch frontend가 기존 timing/traffic/energy accounting 경로에 완전히 연결되는지를 확인하는 것이다.

## 현재 제한과 다음 작업

1. C++ core는 아직 validated executable IR을 Nebula layer object로 임시 변환한다. 완전한 Nebula 제거는 operation descriptor가 DRAM/layer functional storage를 직접 소유한 뒤 진행해야 한다.
2. 현재 adapter는 linear tensor chain만 지원하며 Softmax는 마지막 operation이어야 한다.
3. asymmetric stride, dilation, transposed convolution, nonzero output padding은 adapter 표현 한계 때문에 거부한다.
4. dynamic shape는 capture에는 남길 수 있지만 실행 전에 concrete binding이 필요하다.
5. weight/input 값을 담는 tensor artifact가 없어 IR 경로는 timing-only이다.
6. C++ loader는 declared executable hash를 provenance로 기록하지만 canonical JSON hash 재계산은 Python `validate-executable` 단계가 담당한다.
7. 현재 환경에는 PyTorch package가 없어 실제 `torch.export` 호출은 실행하지 못했다. PyTorch 없이 동작하는 fixture compile, Python validator/lowering, C++ loader, full simulator 실행은 검증했다.

권장 다음 순서는 실제 PyTorch 환경의 Linear/Conv export regression, Nebula와의 frontend differential parity gate, shape-binding artifact, tensor-value artifact, general DAG/alias operation 지원, 마지막으로 Nebula adapter 제거이다.
