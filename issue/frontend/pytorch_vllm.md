# PyTorch frontend 및 vLLM LLM-serving 연동 이슈

## 판정

NPUsim의 `Pytorch` 빌드 플래그와 CPython 임베딩 코드는 **PyTorch frontend가 아니다**. 현재 실행은 항상 Nebula의 `convolutional_t`를 생성하고, 조건부 Python 경로는 `main.build_network()`의 반환 graph/tensor를 NPUsim network나 mapping으로 변환하지 않는다. timing 실행도 convolution과 fully-connected만 대상으로 하므로, PyTorch 모델과 vLLM serving workload는 현재 NPUsim의 입력으로 사용할 수 없다.

이 이슈의 목표는 vLLM을 NPUsim 내부에서 실행하거나 GPU server를 대체하는 것이 아니다. 목표는 **PyTorch model과 vLLM serving workload를 versioned graph IR 및 serving trace로 변환하여 NPUsim의 architecture/mapping/timing model에 재현 가능하게 전달하는 것**이다.

- 심각도: **P0 frontend 계약 부재**, **P0 LLM serving workload 미지원**, **P1 Transformer operator·memory semantics 부재**
- 현재 공식 범위: Nebula 기반 CNN configuration의 Conv/FC timing·energy 추정
- 목표 범위: PyTorch graph ingestion, decoder LLM의 prefill/decode, vLLM scheduling decision을 보존한 serving trace replay

## 확인된 문제

| 우선순위 | 문제 | 영향 |
|---|---|---|
| P0 | `Pytorch=1`은 고정된 Python 3.8 header/library를 연결하고 `main.build_network()`만 호출한다 | PyTorch module, graph dependency, parameter, dtype, dynamic shape가 NPUsim 입력으로 전달되지 않는다 |
| P0 | `npu_t::init()`은 항상 `nebula::convolutional_t`를 생성한다 | frontend와 무관하게 Nebula CNN config parser에 종속된다 |
| P0 | canonical graph/workload IR 및 versioned serialization이 없다 | framework adapter, lowering, mapping generator, simulator 사이의 입력 계약이 없고 수동 cfg/map 작성이 필요하다 |
| P0 | vLLM adapter 또는 serving trace interface가 없다 | request arrival, prompt/generation length, continuous batching, prefill/decode 혼합, KV block event가 timing 입력에서 사라진다 |
| P1 | timing path는 Conv/FC만 실행한다 | embedding, GEMM/batched matmul, normalization, RoPE, masked softmax, attention, activation, sampling 및 KV cache read/write 비용이 누락된다 |
| P1 | layer/mapping이 정적 section 순서로 결합된다 | decode token 증가, variable sequence length, batch compaction에 맞는 shape validation과 mapping 선택이 불가능하다 |
| P1 | KV cache layout, block size, dtype/quantization, residency, eviction/prefix reuse 모델이 없다 | LLM decode의 지배적인 capacity·bandwidth·tail latency를 추정할 수 없다 |
| P1 | serving-level metric과 queue model이 없다 | TTFT, TPOT/ITL, request latency percentile, token throughput/goodput, queue delay 및 admission 결과를 낼 수 없다 |
| P2 | framework/trace/model revision, lowering rule, mapping version provenance가 없다 | PyTorch/vLLM 또는 scheduler policy 변경 뒤 결과 재현·비교가 어렵다 |
| P2 | frontend functional oracle 및 serving-trace golden regression이 없다 | model 변환의 semantic 보존과 scheduling trace 보존을 검증할 수 없다 |

## 코드 근거

[`npu.cc`](../../scheduler/npu.cc#L167)는 network를 `nebula::convolutional_t`로 생성한다. `Pytorch` 분기는 [`npu.cc`](../../scheduler/npu.cc#L169)에서 Python module과 `build_network`를 호출하지만 반환 `pValue`를 읽거나 기존 network를 교체하지 않는다. [`npusim.sh`](../../npusim.sh#L103)은 `/usr/include/python3.8`와 `-lpython3.8`을 직접 지정한다.

실제 accelerator timing은 [`npu.cc`](../../scheduler/npu.cc#L317)에서 convolution/connected만 허용하고, 그 밖의 layer는 [`npu.cc`](../../scheduler/npu.cc#L366)에서 제외한다. 그러므로 Transformer mapping 파일의 존재나 CPython include는 PyTorch·Transformer·vLLM 지원의 근거가 아니다.

## 필요한 frontend 경계

다음 경계를 명시적으로 도입한다.

```text
PyTorch model ─┐
               ├─> capture/export adapter ─> NPUsim graph IR ─> lowering + mapping ─> core simulator
vLLM workload ─┘                              └> serving trace ────────────────────────┘
```

### NPUsim graph IR

frontend IR은 framework object나 Python C API pointer가 아니라 직렬화 가능한 독립 형식이어야 한다. 최소 계약은 다음과 같다.

- graph/node ID, op kind, data/control dependency, source graph location
- tensor ID, logical shape와 symbolic/dynamic dimension binding, layout, dtype, accumulator/quantization metadata
- parameter/activation/KV-cache의 logical byte size와 read/write/alias/in-place semantic
- GEMM transpose, attention head 수·head dimension·causal mask, normalization axis·epsilon, RoPE, activation 같은 op attribute
- model/frontend/adapter/IR schema version, export option, unsupported-op report

shape 또는 dtype를 임의 값으로 채우거나 지원하지 않는 op를 0-cycle로 생략해서는 안 된다. IR validation 단계에서 명시적으로 거부하거나, 사용자가 선택한 conservative fallback cost model로만 처리해야 한다.

### PyTorch adapter

PyTorch adapter는 eager execution 결과를 파싱하기보다 export/capture graph를 입력으로 받아야 한다. `torch.export` 또는 FX 계열 capture를 사용할 수 있으나, 특정 Python ABI나 PyTorch 내부 API가 core simulator ABI에 새면 안 된다.

- Python package/CLI가 model 또는 export artifact와 example/symbolic input contract를 받아 graph IR을 만든다.
- timing-only 실행에서는 weight 값 전체보다 shape/dtype/quantization/parameter identity를 우선 보존하고, value-aware mode에서만 값을 선택적으로 제공한다.
- adapter manifest는 PyTorch/Python version, export mode, model revision/hash, input binding을 기록한다.
- C++ core는 graph IR file만 읽는다. Python embedding은 필수 runtime dependency가 아니어야 한다.

### vLLM serving trace adapter

vLLM adapter의 산출물은 framework 내부 object가 아니라 logical serving trace다. 각 scheduling event에는 최소한 다음이 필요하다.

- request ID, arrival time, prompt length, generated-token target/stop 상태
- prefill 또는 decode phase, event batch의 request 목록, 각 sequence의 position/KV length
- scheduler policy와 option: continuous batching, chunked prefill, priority/admission, prefix-cache reuse 등의 활성 여부와 version
- KV block size, cache dtype/quantization, allocation/free/eviction/reuse event
- tensor/pipeline parallel rank 및 collective가 있다면 logical group과 byte count

NPUsim은 이 trace를 replay하여 event별 shape와 traffic을 simulator에 전달한다. trace는 vLLM Python control-plane cycle을 추정하지 않으며, scheduler가 내린 batch/admission 결정을 workload input으로 받아 accelerator-side cost를 평가한다.

## 단계별 구현 과제와 완료 조건

### P0-A: frontend 계약 및 PyTorch graph ingestion

- versioned graph IR schema와 안정적인 serialization을 정의한다.
- PyTorch adapter CLI/package를 제공해 작은 CNN과 decoder block을 graph IR로 export한다.
- hard-coded Python 3.8 include/link option과 반환값을 버리는 `Pytorch` 경로를 제거하거나 legacy experimental path로 격리한다.
- 모든 graph node에 shape/dtype/layout/source 정보를 기록하고, IR validator가 dangling tensor, unresolved shape, unsupported dtype/op를 거부한다.
- **검증:** PyTorch reference graph와 export IR의 node dependency, tensor shape/dtype, parameter byte count가 golden manifest와 일치한다.

### P0-B: LLM operator lowering과 mapping 연결

- 최소 decoder inference 집합을 명시한다: embedding, linear/GEMM, batched matmul, view/transpose/layout transform, residual elementwise, RMSNorm/LayerNorm, RoPE, causal masked softmax, attention, GELU/SiLU 및 sampling boundary.
- 각 op를 NPUsim primitive/compound workload로 lowering하고, 지원 상태를 `modeled`, `fallback`, `rejected`로 결과에 표시한다.
- GEMM/attention의 M/N/K 및 batch/head/sequence dimension이 mapping/config와 자동 검증되게 한다. mapping이 없으면 추정 결과를 내지 말고 필요한 mapping을 보고한다.
- view/reshape처럼 arithmetic cost가 없는 op도 layout 변경 또는 materialization 여부를 명시한다.
- **검증:** synthetic decoder block과 공개 가능한 소형 decoder model에서 op coverage가 100%이거나, 미지원 op가 명시적으로 실패한다. 0-cycle silent skip은 허용하지 않는다.

### P0-C: KV cache 및 prefill/decode workload

- KV cache를 별도 logical memory object로 모델링하며 layer/head/sequence/block layout, dtype/quantization, read/write/append byte를 기록한다.
- prefill과 decode를 별도 event type으로 모델링하고, decode step마다 sequence length 증가가 attention shape·KV traffic에 반영되게 한다.
- cache capacity, residency tier(GLB/DRAM 등), allocation/free/eviction/prefix reuse의 모델 범위를 config에 명시한다.
- **검증:** 1 request prefill→여러 decode, 서로 다른 prompt length의 2 request continuous batching, capacity overflow/eviction fixture에서 KV byte와 event 순서가 기대 trace와 일치한다.

### P0-D: vLLM serving trace replay와 metric

- vLLM version과 scheduler option을 manifest에 기록하는 trace exporter를 제공한다.
- NPUsim trace replay가 batch composition, phase, token position, KV event를 손실 없이 읽고 shape/mapping validation을 수행한다.
- request/token 단위로 arrival, queue, dispatch, prefill completion, 각 decode completion을 기록한다.
- TTFT, TPOT/ITL, end-to-end request latency percentile, token throughput/goodput, batch size, KV capacity/high-water mark를 산출한다.
- **검증:** 고정 trace 재실행의 deterministic output, request 순서·token count 보존, batch별 shape/KV byte golden comparison을 CI regression으로 추가한다.

### P1: serving-system 및 multi-device 확장

- tensor/pipeline parallel collective, MoE routing/expert dispatch, LoRA adapter residency, speculative decoding의 trace/IR semantics를 추가한다.
- router/link contention과 multi-chip timing이 검증된 뒤 collective communication을 hardware model에 연결한다.
- scheduler/control-plane overhead를 별도 host model로 추가할지, external trace timestamp에 포함할지 결과 문서에서 선택 가능하게 한다.

## 공통 완료 기준

- 결과 파일에는 graph IR schema, PyTorch/vLLM/adapter revision, model hash, trace hash, accelerator config 및 mapping revision을 기록한다.
- graph IR validation, operator lowering, KV cache accounting, serving trace replay 각각에 독립 micro-test와 golden regression이 존재한다.
- shape/dtype/op/mapping 불일치는 명시적 오류 또는 명시적 fallback으로 처리되며, timing total에서 묵시적으로 제외되지 않는다.
- functional mode에서는 선택한 소형 fixture에 대해 PyTorch reference와 layer/KV-cache boundary의 값 또는 허용 오차를 비교한다. timing-only mode에서는 값 없이도 동일한 shape/traffic trace를 재현한다.
- 대표 workload는 synthetic decoder, small open decoder prefill/decode, variable-length continuous-batching trace를 포함한다. 실제 모델 가중치와 dataset은 라이선스·배포 제약에 따라 외부 download fixture로 분리한다.

## 비목표와 해석 주의

- vLLM의 GPU kernel 성능 또는 Python scheduler의 구현 시간을 NPUsim이 자동으로 대체·예측한다는 주장은 이 범위에 포함하지 않는다.
- value-aware inference와 exact token output은 별도 functional/precision 검증이 선행되어야 하며, timing trace replay만으로 vLLM output bit-exactness를 보장하지 않는다.
- PyTorch graph export가 성공했다는 사실만으로 모든 op가 hardware timing에 반영됐다는 뜻은 아니다. 결과의 operator coverage와 fallback/rejection report가 함께 제공돼야 한다.
