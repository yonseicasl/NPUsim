# Nebula에서 PyTorch Frontend로의 전환 계획

## 목적

현재 NPUSim의 network·layer·functional 입력 경로를 담당하는 Nebula 의존성을 제거하고, PyTorch model을 재현 가능한 framework-independent artifact로 변환하여 NPUSim의 scheduler와 timing/traffic/energy simulator가 직접 소비할 수 있도록 한다.

이 전환의 목표는 PyTorch 또는 libtorch를 NPUSim C++ process에 직접 embed하는 것이 아니다. PyTorch는 model capture와 reference 실행을 담당하고, NPUSim core는 PyTorch가 설치되지 않은 환경에서도 versioned IR artifact만으로 실행될 수 있어야 한다.

권장 전체 경로는 다음과 같다.

```text
PyTorch model
    │ torch.export
    ▼
Capture Graph IR
    │ NPUSim lowering
    ▼
Executable IR + shape binding + mapping binding
    │
    ▼
NPUSim C++ timing/traffic/energy core
```

## 현재 상태

현재 repository에는 PyTorch frontend의 첫 단계가 이미 존재한다.

- [`frontend/pytorch/export.py`](../frontend/pytorch/export.py): `torch.export` 기반 graph capture
- [`frontend/pytorch/graph_ir.py`](../frontend/pytorch/graph_ir.py): `npusim.graph.v1` schema와 validator
- [`frontend/pytorch/cli.py`](../frontend/pytorch/cli.py): graph export와 validation CLI
- [`unittest/python/test_pytorch_frontend.py`](../unittest/python/test_pytorch_frontend.py): graph IR 기본 regression
- [`frontend/fixtures/minimal_graph.json`](../frontend/fixtures/minimal_graph.json): PyTorch 없이 검증 가능한 최소 fixture

현재 frontend는 tensor shape, dtype, logical byte, dependency, operation name과 provenance를 JSON artifact로 저장한다. 그러나 이 artifact는 아직 C++ simulator가 소비하지 않는다.

### 남아 있는 Nebula 결합

현재 C++ core는 다음 지점에서 Nebula에 직접 의존한다.

- [`scheduler/npu.cc`](../scheduler/npu.cc)는 `nebula::convolutional_t`를 생성하고 Nebula network config를 읽는다.
- [`scheduler/npu.h`](../scheduler/npu.h)는 `nebula::network_t*`와 `nebula::layer_t*`를 core state로 보유한다.
- [`components/dram.h`](../components/dram.h)는 `connect_layer(nebula::layer_t*)`로 tensor pointer를 전달받는다.
- functional mode는 `network->load_data()`와 `layer->forward()`에 의존한다.
- accelerator 실행은 `network->layers[index]`와 mapping section index를 직접 결합한다.
- [`models/main.cc`](../models/main.cc)는 network 이름으로 Nebula config 경로를 구성한다.
- [`npusim.sh`](../npusim.sh)는 Nebula header/library를 include/link하고 Nebula build/download를 관리한다.
- [`models/setup.sh`](../models/setup.sh)는 Nebula weight/dataset downloader를 사용한다.

Nebula 의존은 주로 `npu_t`, `dram_t`, executable entry point와 build script에 집중되어 있다. PE, PE array, global buffer, multi-chip, stats의 핵심 timing model은 대부분 mapping/scheduler descriptor를 사용하므로 framework-neutral operation descriptor를 도입하면 재사용할 수 있다.

### 현재 PyTorch IR의 한계

현재 `npusim.graph.v1` exporter는 graph topology를 보존하지만 simulator 실행에 필요한 다음 정보가 부족하다.

- convolution stride, padding, dilation, groups
- softmax axis, Leaky ReLU slope 등 scalar/list literal argument
- tensor stride, storage offset, memory format과 alias 관계
- view/reshape가 no-cost alias인지 materialization인지에 대한 정보
- structured/multi-output operation
- parameter qualified name과 값 또는 value artifact 위치
- dynamic shape symbol의 concrete simulation binding
- operation support/lowering 상태
- C++ loader와 executable operation geometry
- mapping과 graph operation의 stable ID binding

현재 모든 tensor layout이 `contiguous`로 기록되며, weight 값은 저장하지 않는다. 따라서 지금 상태의 graph IR을 C++에서 바로 읽더라도 올바른 Conv/GEMM geometry와 functional tensor를 재구성할 수 없다.

## 전환 원칙

1. **PyTorch는 export-time dependency로 한정한다.** C++ core는 CPython, libtorch, Python ABI에 의존하지 않는다.
2. **Capture와 executable 의미를 분리한다.** PyTorch graph를 보존하는 IR과 hardware simulator가 실행하는 canonical IR을 별도로 둔다.
3. **Framework operation을 core에 노출하지 않는다.** C++ core는 `aten.*` 문자열이 아니라 `npusim.conv2d`, `npusim.gemm` 같은 canonical operation만 이해한다.
4. **Operation ID로 mapping을 결합한다.** network/mapping section의 ordinal index 결합을 제거한다.
5. **Unsupported operation은 fail-fast한다.** 지원되지 않는 node를 묵시적으로 0-cycle skip하지 않는다.
6. **Timing-only와 functional mode를 분리한다.** Timing-only 실행에는 weight 값이 필수가 아니며, functional 검증은 별도 tensor artifact와 PyTorch oracle을 사용한다.
7. **Frontend 교체는 hardware 결과를 바꾸지 않는다.** 동일 operation geometry와 mapping이면 cycle, traffic, energy가 기존 Nebula 경로와 동일해야 한다.
8. **Nebula는 parity가 확보된 뒤 제거한다.** 초기 migration 동안 legacy frontend로 유지하고 differential validation에 사용한다.

## 목표 아키텍처

### 1. Capture Graph IR

현재 `npusim.graph.v1`을 PyTorch capture artifact로 유지하거나 호환 가능한 다음 schema version으로 확장한다.

Capture IR의 책임은 PyTorch export graph를 손실 없이 직렬화하는 것이다.

최소 필드는 다음과 같다.

- graph/node/tensor stable ID
- raw source operation과 source location
- tensor dependency와 control dependency
- tensor logical shape, symbolic dimension, dtype
- physical stride, storage offset, layout/memory format
- tensor kind: input, parameter, buffer, constant, activation
- tensor alias/view/in-place relationship
- parameter qualified name, logical byte, checksum
- tensor operand와 literal argument의 분리
- positional/keyword literal argument
- PyTorch/Python/frontend version
- model structure hash와 export option
- unsupported structured result 정보

Capture IR은 timing model 지원 여부를 결정하지 않는다. 이는 다음 lowering 단계의 책임이다.

### 2. Executable IR

새 schema `npusim.exec.v1` 또는 동등한 executable artifact를 정의한다.

Executable IR은 C++ core가 직접 읽는 framework-neutral 입력이다.

```json
{
  "schema_version": "npusim.exec.v1",
  "operations": [
    {
      "id": "features_0_conv",
      "kind": "npusim.conv2d",
      "inputs": ["input", "conv_weight", "conv_bias"],
      "outputs": ["conv_output"],
      "geometry": {
        "batch": 1,
        "input_channels": 3,
        "input_height": 224,
        "input_width": 224,
        "output_channels": 64,
        "filter_height": 7,
        "filter_width": 7,
        "stride_height": 2,
        "stride_width": 2,
        "padding_height": 3,
        "padding_width": 3,
        "dilation_height": 1,
        "dilation_width": 1,
        "groups": 1
      }
    }
  ]
}
```

초기 canonical operation 집합은 현재 NPUSim timing 범위에 맞춰 작게 시작한다.

| Canonical operation | 초기 처리 |
|---|---|
| `npusim.conv2d` | 기존 convolution scheduler로 lowering |
| `npusim.linear` / `npusim.gemm` | 기존 connected/GEMM scheduler로 lowering |
| `npusim.relu` | SFU operation 또는 허용된 fused activation |
| `npusim.leaky_relu` | SFU operation 또는 허용된 fused activation |
| `npusim.softmax` | standalone SFU multi-pass operation |
| `npusim.reshape` / `npusim.view` | alias이면 no-cost, materialization이면 traffic operation |
| `npusim.transpose` | layout-only 또는 materialization을 명시 |
| `npusim.unsupported` | simulation 시작 전 명시적 rejection |

Executable IR의 operation은 반드시 다음 상태 중 하나로 분류한다.

- `modeled`: cycle/traffic/energy 모델 존재
- `no_cost_alias`: storage alias로 hardware event 없음
- `external_boundary`: sampling 등 accelerator simulation 경계
- `rejected`: 현재 지원되지 않음

묵시적 skip 상태는 허용하지 않는다.

### 3. Shape Binding

PyTorch export는 symbolic dimension을 보존할 수 있지만, NPUSim의 한 번의 timing 실행에는 concrete geometry가 필요하다.

별도 shape binding artifact를 둔다.

```json
{
  "batch": 4,
  "sequence_length": 128,
  "symbols": {
    "s0": 4,
    "s1": 128
  }
}
```

Unresolved symbol이 하나라도 있으면 mapping과 timing을 시작하지 않는다. 여러 shape를 평가할 때는 동일 Capture IR에 여러 binding을 적용한다.

### 4. Optional Tensor Artifact

Timing-only mode는 weight 값 없이 tensor metadata만 사용한다. Functional mode에서는 별도의 tensor artifact를 사용한다.

권장 tensor artifact 구성은 다음과 같다.

- tensor manifest: ID, qualified name, dtype, shape, offset, byte size, checksum
- binary tensor payload: input/parameter/reference output
- model/frontend/version/hash provenance

C++ core는 이 artifact 포맷만 읽고 PyTorch object를 직접 받지 않는다.

## 단계별 실행 계획

### Phase 0 — 범위 고정과 differential baseline

1. 현재 Nebula 경로의 Conv/FC/SFU 결과를 golden baseline으로 저장한다.
2. PyTorch model로 동일한 geometry를 만드는 small fixture를 준비한다.
3. frontend와 무관해야 하는 비교 항목을 고정한다.
4. 기존 network layer index와 mapping index 결합 사례를 목록화한다.
5. 현재 지원 layer와 timing scope를 명시한다.

필수 baseline workload는 다음과 같다.

- 작은 Conv2d
- stride/padding이 있는 Conv2d
- grouped/depthwise Conv2d
- Linear/GEMM
- Conv+ReLU
- Conv+Leaky ReLU
- Linear+Softmax
- mapping padding이 있는 Conv

비교 항목:

- B/K/P/Q/C/R/S와 groups/stride/padding
- tile size와 repetition
- MAC 및 reduction operation 수
- input/weight/output traffic
- layer critical-path cycle
- component별 dynamic/static energy
- SFU element/invocation 수

완료 조건: Nebula baseline과 PyTorch fixture의 semantic geometry가 문서화되고, 이후 migration 결과를 자동 비교할 golden이 존재한다.

### Phase 1 — Capture Graph IR 강화

1. [`frontend/pytorch/export.py`](../frontend/pytorch/export.py)에서 tensor operand와 literal argument를 분리한다.
2. Conv/GEMM/activation/softmax의 positional·keyword attribute를 보존한다.
3. tensor stride, storage offset, memory format을 추가한다.
4. alias/view/in-place semantics를 추가한다.
5. parameter qualified name과 structural/value checksum을 구분한다.
6. tuple/list/dict result를 explicit tensor output으로 표현한다.
7. symbolic shape와 concrete example size를 혼동하지 않도록 binding contract를 추가한다.
8. schema version migration과 backward compatibility validator를 추가한다.

예시 node:

```json
{
  "id": "conv2d",
  "op": "aten.convolution.default",
  "tensor_inputs": ["x", "weight", "bias"],
  "attributes": {
    "stride": [1, 1],
    "padding": [1, 1],
    "dilation": [1, 1],
    "transposed": false,
    "output_padding": [0, 0],
    "groups": 1
  }
}
```

완료 조건: PyTorch reference graph와 Capture IR의 dependency, tensor shape/dtype/layout, Conv attribute, parameter byte가 golden manifest와 일치한다.

### Phase 2 — Python lowering과 executable compiler

다음 module을 추가한다.

- `frontend/pytorch/op_registry.py`
- `frontend/pytorch/lowering.py`
- `frontend/pytorch/executable_ir.py`
- `frontend/pytorch/compile.py`

Lowering의 책임은 다음과 같다.

1. PyTorch/ATen operation을 stable NPUSim canonical operation으로 변환한다.
2. Input/output shape를 독립적으로 재계산하고 exported metadata와 대조한다.
3. Conv/Linear를 기존 mapping parameter B/K/P/Q/C/R/S로 변환한다.
4. View/reshape/transpose의 alias 또는 materialization 여부를 판정한다.
5. activation fusion을 안전한 경우에만 수행한다.
6. 지원 operation coverage report를 만든다.
7. unsupported operation이 있으면 executable artifact 생성을 실패시킨다.
8. source node와 lowered operation의 provenance mapping을 저장한다.

Activation fusion은 다음 조건을 모두 만족할 때만 허용한다.

- Conv/Linear output의 consumer가 하나
- 다음 operation이 지원되는 activation
- 중간 output에 branch 또는 alias consumer가 없음
- activation 순서가 SFU placement/output cast contract와 일치
- fusion 전후 dtype/layout 의미가 동일

Fusion된 operation은 source node 목록과 적용된 rule/version을 기록한다. 조건을 만족하지 않으면 standalone SFU operation으로 유지한다.

완료 조건: Linear/ReLU와 Conv/ReLU fixture가 canonical executable IR로 변환되고, geometry와 fusion provenance가 golden과 일치한다.

### Phase 3 — Framework-neutral C++ workload model

다음 C++ module을 추가한다.

- `scheduler/workload_graph.h`
- `scheduler/workload_graph.cc`
- `scheduler/graph_ir_loader.h`
- `scheduler/graph_ir_loader.cc`

목표 core type은 다음과 같다.

```cpp
struct tensor_desc_t {
    std::string id;
    std::vector<size_t> shape;
    tensor_format_t format;
    tensor_kind_t kind;
    size_t logical_bytes;
};

struct operation_desc_t {
    std::string id;
    operation_kind_t kind;
    std::vector<std::string> inputs;
    std::vector<std::string> outputs;
    operation_geometry_t geometry;
    activation_op_t activation;
};

struct workload_graph_t {
    std::vector<tensor_desc_t> tensors;
    std::vector<operation_desc_t> operations;
};
```

Loader는 다음을 검증한다.

- schema/version compatibility
- graph hash와 artifact checksum
- tensor/node ID uniqueness
- topological order와 dangling dependency
- tensor shape/dtype/logical-byte identity
- operation별 required attribute
- symbolic dimension의 concrete binding
- supported operation과 mapping requirement

JSON parser는 개발 환경의 우연한 Python/miniforge header에 의존하지 않는다. Repository에 고정된 parser를 vendoring하거나 명시적인 system/build dependency로 관리하고 version을 기록한다.

완료 조건: PyTorch 없이 C++ unit test가 `minimal_graph` 또는 executable fixture를 읽고 Linear operation geometry를 재구성한다.

### Phase 4 — NPU와 DRAM의 Nebula type 제거

[`scheduler/npu.h`](../scheduler/npu.h)의 다음 state를 framework-neutral type으로 바꾼다.

```cpp
// 제거 대상
nebula::network_t *network;
nebula::layer_t *layer;

// 목표
std::unique_ptr<workload_graph_t> workload;
const operation_desc_t *operation;
```

[`scheduler/npu.cc`](../scheduler/npu.cc)의 실행 순서는 다음과 같이 변경한다.

```text
load executable IR
→ validate graph and concrete shapes
→ bind operation IDs to mappings
→ iterate topological operations
→ configure scheduler from operation geometry
→ execute accelerator components
→ execute SFU/post-processing operation
→ update per-operation/network stats
```

제거할 Nebula 호출:

- `network->init()`
- `network->load_data()`
- `network->num_layers`
- `network->layers[index]`
- `network->input_data`
- `layer->forward()`
- `nebula::activation_type_str`

[`components/dram.h`](../components/dram.h)의 `connect_layer(nebula::layer_t*)`는 descriptor 기반 API로 바꾼다.

```cpp
void bind_operation(
    const operation_desc_t &operation,
    const tensor_desc_t &input,
    const tensor_desc_t *weight,
    const tensor_desc_t &output);
```

Timing-only path에서는 실제 `input_data`, `weight`, `output_data` pointer를 사용하지 않는다. Tile shape, datatype descriptor, mapping과 action type만으로 access/transaction을 계산한다.

Functional path의 tensor pointer는 optional tensor artifact가 활성화된 경우에만 별도 tensor store에서 얻는다.

완료 조건: minimal Linear/GEMM executable IR이 Nebula network/layer object 없이 기존 scheduler와 component timing path를 실행한다.

### Phase 5 — Mapping을 stable operation ID에 결합

현재 network layer index와 mapping section index의 직접 결합을 제거한다.

Mapping section에 stable operation ID를 추가한다.

```ini
[layer]
op_id = encoder_layer_0_mlp_fc1
op_kind = linear
```

Binding 규칙:

- executable operation ID로 mapping을 검색한다.
- operation kind와 mapping layer kind가 일치해야 한다.
- B/K/P/Q/C/R/S coverage를 operation geometry와 검증한다.
- mapping이 필요한 modeled operation에 mapping이 없으면 fail-fast한다.
- alias/view no-cost operation에는 mapping이 필요하지 않다.
- SFU operation은 SFU profile 또는 별도 SFU mapping contract를 사용한다.
- ordinal index fallback은 migration 기간의 legacy mode에서만 허용한다.

향후 동일한 operation에 shape별 mapping을 적용할 수 있도록 `op_id + shape/profile selector` 구조를 고려한다.

완료 조건: graph에 activation/view node가 삽입되거나 fusion되더라도 Conv/Linear mapping binding이 바뀌지 않는다.

### Phase 6 — CLI와 실행 artifact 전환

현재 실행 형식은 accelerator/network/mapping 이름에 의존한다. 새 CLI는 파일 artifact를 명시적으로 받도록 한다.

```bash
python3 -m frontend.pytorch.cli export \
  --factory my_model:build \
  --output artifacts/model.graph.json

python3 -m frontend.pytorch.cli compile \
  --graph artifacts/model.graph.json \
  --shape-binding artifacts/shape.json \
  --output artifacts/model.exec.json

./models/model run \
  --accelerator configs/accelerators/gemmini.cfg \
  --graph artifacts/model.exec.json \
  --mapping configs/mappings/gemmini/model.map
```

결과에는 다음 provenance를 기록한다.

- Capture IR hash
- Executable IR hash
- PyTorch/frontend/lowering version
- model structure/parameter hash
- shape binding hash
- mapping/accelerator config hash
- operation coverage와 fusion summary

Legacy positional CLI는 migration 기간에 wrapper로 유지할 수 있지만 내부에서는 동일 executable workload interface를 사용한다.

완료 조건: PyTorch가 없는 환경에서 export된 executable IR과 mapping만으로 NPUSim timing 실행이 가능하다.

### Phase 7 — Functional validation 경로

C++ core에 libtorch를 link하지 않고 PyTorch를 외부 oracle로 사용한다.

1. PyTorch exporter/reference runner가 input, parameter, reference output tensor artifact를 생성한다.
2. NPUSim functional mode가 artifact의 tensor manifest와 binary payload를 읽는다.
3. NPUSim component output 또는 operation boundary output을 별도 artifact로 기록한다.
4. Python checker가 PyTorch reference와 NPUSim output을 비교한다.
5. integer는 bit-exact, floating point는 profile별 absolute/relative/ULP tolerance를 사용한다.

초기에는 final output과 Conv/Linear/SFU boundary를 비교한다. 모든 intermediate를 항상 저장하면 artifact가 커지므로 debug/validation option으로 제한한다.

Nebula `layer->forward()`를 제거하기 전까지 동일 output을 Nebula와 PyTorch가 함께 수정하지 않도록 functional owner를 명확히 한다.

완료 조건: 작은 Linear, Conv, activation fixture의 NPUSim functional output이 PyTorch reference와 지정 tolerance 내에서 일치한다.

### Phase 8 — DAG와 operator 범위 확장

현재 실행은 이전 layer output을 다음 layer input으로 전달하는 선형 chain을 가정한다. PyTorch graph를 일반적으로 지원하려면 tensor ID 기반 DAG execution이 필요하다.

추가 과제:

- multiple input/output
- branch와 residual connection
- tensor reference count와 liveness
- concat
- element-wise add/multiply
- pooling
- BatchNorm folding 또는 standalone normalization
- LayerNorm/RMSNorm
- GELU/SiLU
- batched matmul
- embedding
- attention과 KV cache

View/reshape/transpose는 arithmetic operation 수가 없더라도 layout materialization과 traffic이 있는지 구분한다. BatchNorm folding은 eval-mode parameter와 provenance가 있을 때만 수행한다.

완료 조건: 대표 CNN residual block이 graph dependency를 보존하며 실행되고, 모든 node가 modeled/no-cost/external/rejected 중 하나로 분류된다.

### Phase 9 — Nebula dependency 제거

PyTorch/executable IR 경로가 기존 baseline과 parity를 달성한 뒤 다음 의존을 제거한다.

- `scheduler/npu.h`의 Nebula network/layer header
- `components/dram.h`의 Nebula layer header
- `models/main.cc`의 Nebula include
- `npusim.sh`의 Nebula clone/build/update 단계
- Nebula include path와 `-lnebula`
- Nebula 때문에 연결된 불필요한 OpenBLAS/OpenCV library
- `models/setup.sh`의 Nebula dataset/weight downloader
- 최종적으로 `ext/nebula`

기존 CNN regression은 다음 중 하나로 보존한다.

1. Nebula network config를 executable IR로 변환하는 일회성 converter
2. 기존 validation workload를 고정된 executable IR fixture로 변환

Nebula directory 삭제는 별도 destructive migration으로 수행하고, 모든 baseline과 build matrix가 통과한 뒤 진행한다.

완료 조건: 깨끗한 환경에서 Nebula를 clone/build/link하지 않고 NPUSim library, executable, unit test와 대표 validation을 실행할 수 있다.

## 수정 대상 파일과 module

| 대상 | 주요 변경 |
|---|---|
| `frontend/pytorch/export.py` | literal argument, layout, alias, parameter provenance export |
| `frontend/pytorch/graph_ir.py` | 확장 schema와 version migration validation |
| `frontend/pytorch/cli.py` | compile/lower/reference artifact command |
| `frontend/pytorch/op_registry.py` | ATen→NPUSim canonical op registry |
| `frontend/pytorch/lowering.py` | geometry 추출, fusion, unsupported-op gate |
| `frontend/pytorch/executable_ir.py` | executable schema와 validator |
| `scheduler/workload_graph.{h,cc}` | framework-neutral tensor/operation graph |
| `scheduler/graph_ir_loader.{h,cc}` | C++ executable IR loader |
| `scheduler/npu.{h,cc}` | Nebula network/layer state 제거, operation iteration |
| `scheduler/scheduler.{h,cc}` | operation geometry에서 mapping parameter 구성 |
| `scheduler/mapping_table.{h,cc}` | stable `op_id` binding과 coverage validation |
| `components/dram.{h,cc}` | layer pointer를 tensor/operation descriptor로 대체 |
| `components/sfu.{h,cc}` | canonical activation/softmax operation 연결 |
| `scheduler/stats.{h,cc}` | graph/operation coverage와 provenance report |
| `models/main.cc` | graph artifact 기반 CLI |
| `npusim.sh` | migration build mode, 최종 Nebula link 제거 |
| `models/setup.sh` | PyTorch artifact/reference fixture workflow |
| `unittest/python/*` | export/lowering/schema test |
| `unittest/*` | C++ IR loader, operation binding, differential test |
| `validation/frontend/*` | Nebula↔PyTorch parity와 unsupported-op regression |

## 테스트 전략

### Frontend schema test

- static/dynamic shape
- dtype와 logical byte identity
- Conv stride/padding/dilation/groups
- Linear transpose semantics
- Softmax axis와 activation attribute
- parameter/buffer/constant 구분
- alias/view/materialization
- multi-output operation
- dangling tensor, duplicate ID, unresolved symbol rejection

### Lowering test

- ATen Conv/Linear to canonical geometry
- Conv+ReLU/Leaky fusion의 positive/negative case
- branch가 있는 activation의 fusion 금지
- no-cost view와 materialized transpose 구분
- unsupported operation compile failure
- source graph와 executable operation provenance identity

### C++ loader test

- minimal executable IR load
- schema/hash/checksum mismatch rejection
- operation required attribute validation
- tensor dependency/topological order validation
- shape binding 적용
- PyTorch가 설치되지 않은 환경에서 실행

### Differential timing test

동일 operation geometry와 mapping에 대해 Nebula와 PyTorch IR 경로를 비교한다.

- tile/repetition
- MAC/reduction/SFU event count
- component traffic
- component busy cycle와 layer latency
- dynamic/static energy
- network rollup

Frontend만 달라진 경우 수치는 exact match를 기본으로 한다. 차이가 발생하면 hardware model 변경이 아니라 geometry/lowering/mapping 차이로 분류하고 source operation까지 추적한다.

### Functional test

- Linear와 Conv output
- ReLU/Leaky/Softmax output
- dtype conversion과 saturation/rounding
- branch/residual output
- tensor artifact checksum
- PyTorch reference tolerance

### Negative test

- unsupported ATen operation
- missing mapping
- operation kind와 mapping kind 불일치
- unresolved dynamic dimension
- invalid tensor byte/layout
- unsafe activation fusion
- weight-required functional run에서 tensor artifact 누락
- graph node의 silent timing exclusion 방지

## 호환성 및 migration 정책

- 초기에는 `legacy-nebula`와 `pytorch-ir` frontend를 모두 지원한다.
- 두 frontend는 가능한 한 동일 `workload_graph_t`와 scheduler 경로로 수렴시킨다.
- 기존 Nebula config의 cycle/traffic/energy baseline은 PyTorch parity gate로 사용한다.
- `[sfu]` 및 datatype profile은 framework frontend와 무관한 accelerator config로 유지한다.
- 기존 mapping의 ordinal binding은 deprecated warning 후 stable `op_id` binding으로 이전한다.
- PyTorch version 차이로 ATen graph가 달라질 수 있으므로 lowering registry와 export version을 결과에 기록한다.
- Capture IR schema와 Executable IR schema는 독립적으로 versioning한다.
- Weight 값을 포함하지 않은 timing-only artifact는 functional fidelity 근거로 사용하지 않는다.

## 위험과 의사결정 게이트

1. **ATen decomposition 불안정성**: PyTorch version마다 operation graph가 달라질 수 있으므로 core가 ATen을 직접 해석하지 않고 Python lowering registry에서 정규화한다.
2. **Literal argument 손실**: stride/padding/groups가 Capture IR에 보존되기 전에는 Conv timing을 C++ IR 경로에 연결하지 않는다.
3. **Mapping index 불일치**: operation ID binding이 구현되기 전에는 일반 PyTorch graph를 기존 mapping section 순서와 직접 결합하지 않는다.
4. **DAG를 선형 layer로 축소**: branch/residual graph를 이전 node output chain으로 변환하지 않는다. 지원 전에는 명시적으로 거부한다.
5. **Layout traffic 누락**: reshape/transpose가 alias인지 materialization인지 판정할 수 없으면 no-cost로 처리하지 않는다.
6. **Functional owner 중복**: Nebula/PyTorch/NPUSim이 동일 activation이나 layer output을 함께 수정하지 않는다.
7. **Frontend와 hardware 변경 혼합**: migration PR에서는 hardware timing model 수정을 분리하고 differential baseline으로 frontend 차이만 검증한다.
8. **성급한 Nebula 제거**: representative parity와 clean build gate가 통과하기 전에 `ext/nebula`를 삭제하지 않는다.

## 권장 구현 순서

```text
Nebula differential baseline
→ Capture IR attribute/layout 보강
→ Linear/ReLU executable lowering
→ C++ executable IR loader
→ stable op_id mapping binding
→ Linear timing exact parity
→ Conv2d/fused-SFU timing exact parity
→ functional tensor artifact + PyTorch oracle
→ DAG/operator 확장
→ legacy config fixture 변환
→ Nebula include/link/build 제거
```

## 첫 번째 milestone

첫 구현 milestone은 다음 범위로 제한한다.

```text
PyTorch Linear + ReLU model
→ torch.export Capture IR
→ npusim.linear + fused/standalone ReLU Executable IR
→ C++ IR loader
→ 기존 GEMM mapping을 op_id로 binding
→ cycle/traffic/energy 출력
→ 동일 Nebula fixture와 exact comparison
```

이 milestone의 완료 조건:

- C++ core가 PyTorch 없이 executable IR을 읽는다.
- Linear의 M/N/K와 tensor byte가 PyTorch model과 일치한다.
- ReLU가 fusion rule에 따라 SFU event로 정확히 한 번 반영된다.
- Mapping은 network index가 아니라 operation ID로 연결된다.
- 기존 Nebula fixture와 cycle/traffic/energy가 exact match한다.
- unsupported operation은 simulation 시작 전에 실패한다.

첫 milestone이 통과한 뒤 Conv2d, grouped convolution, standalone softmax 순서로 확장한다.

## 최종 완료 조건

- PyTorch model이 versioned Capture IR과 Executable IR로 변환된다.
- 모든 simulation operation은 canonical NPUSim operation으로 표현된다.
- Operation attribute, tensor shape/dtype/layout/alias와 mapping이 손실 없이 연결된다.
- C++ core가 `nebula::network_t`와 `nebula::layer_t` 없이 동작한다.
- Timing-only 실행은 PyTorch와 weight 값 없이 재현 가능하다.
- Functional validation은 tensor artifact와 PyTorch oracle로 수행된다.
- Unsupported operation과 unresolved shape/mapping이 silent skip되지 않는다.
- 동일 workload/mapping의 Nebula와 PyTorch IR cycle/traffic/energy가 parity를 만족한다.
- 결과에 graph, lowering, shape, mapping, accelerator config provenance가 기록된다.
- 최종 build와 representative validation이 Nebula clone/include/link 없이 통과한다.

## 2026-08-27 구현 상태

첫 vertical slice는 구현 완료했다.

- `npusim.graph.v1`에서 `npusim.exec.v1`로의 strict lowering
- Linear, Conv2d, fused ReLU/Leaky ReLU, terminal Softmax canonicalization
- C++ executable IR loader와 tensor/geometry 교차 검증
- mapping `op_id` binding과 잘못된 binding의 fail-fast
- `model run-ir` timing/traffic/energy 실행 경로
- Linear+ReLU Gemmini+SFU end-to-end 실행 및 provenance 출력

Nebula differential exact-parity gate, dynamic shape binding, tensor-value artifact, 일반 DAG/alias 처리와 Nebula 제거는 아직 완료되지 않았다. 상세 결과와 재현 명령은 `implementation/PYTORCH_FRONTEND_MILESTONE_2026-08-27.md`에 기록했다.
