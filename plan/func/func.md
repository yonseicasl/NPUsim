# NPUsim Compression/Decompression 실험 계획

## 1. 실험 목표

본 실험은 LLM weight compression의 이득과 decompression engine의 비용을 함께 모델링하여, compression이 실질적인 성능 향상으로 이어지는 조건을 분석하는 것을 목표로 한다.

### 연구 질문

- Compression ratio가 증가하면 성능은 얼마나 향상되는가?
- 어느 지점부터 decompression engine이 새로운 병목이 되는가?
- Compression의 이득을 유지하기 위해 필요한 decompression throughput은 얼마인가?
- Decode와 prefill에서 compression의 효과는 어떻게 달라지는가?
- DRAM bandwidth에 따라 compression의 효과가 어떻게 달라지는가?

### 핵심 메시지

> LLM weight compression은 메모리 용량과 off-chip bandwidth 요구량을 줄이지만, 실제 성능은 compression ratio뿐 아니라 decompression throughput과 compute overlap에 의해 결정된다.

---

## 2. 모델링 범위

현재 NPUsim은 CSR, CSC, SparseMap과 같은 compressed data format과 zero skipping을 지원한다. 이번 실험에서는 compressed weight를 복원하는 별도의 hardware component를 추가한다.

~~~text
DRAM
  ↓ Compressed weights + metadata
Global/Weight Buffer
  ↓ Compressed tile
DECA-inspired Decompression Engine
  ↓ Dense tile
Local Buffer / PE Array
~~~

### Decompression engine 파라미터

- Startup latency
- Decode throughput
- Decompression engine 개수
- Input/output queue 크기
- Metadata 및 scale 크기
- DRAM transfer와 decompression의 overlap
- 압축 효율이 낮을 때 dense format을 사용하는 bypass 기능

초기 구현에서는 weight compression만 대상으로 하며 activation compression은 제외한다.

---

## 3. Compression Ratio 정의

Compression ratio는 metadata와 quantization scale을 포함하여 계산한다.

$$
CR =
\frac{\text{Dense weight bytes}}
{\text{Compressed value bytes}
+\text{Metadata bytes}
+\text{Scale bytes}}
$$

### 데이터 크기 적용 범위

- DRAM → Global Buffer: compressed size
- Global Buffer 저장 공간: compressed size
- Decompression Engine → Local Buffer: dense size
- Local Buffer → PE: dense size
- MAC operation count: dense baseline과 동일

DECA-inspired engine이 dense tile을 출력하므로, 초기 실험에서는 zero-skipping에 의한 연산량 감소를 적용하지 않는다. 이를 통해 compression에 따른 memory traffic 감소와 decompression overhead를 독립적으로 분석한다.

---

## 4. 구현 계획

### 4.1 Decompression Engine Interface

공통 `DecompressionEngine` base class를 구현한다.

주요 기능은 다음과 같다.

- `compressed_size()`
- `metadata_size()`
- `decode_cycles()`
- `can_accept()`
- `is_output_ready()`

### 4.2 Decompression Engine 종류

- `None`: 압축하지 않은 dense baseline
- `Ideal`: decompression overhead가 없는 upper bound
- `DECA`: latency와 throughput을 반영한 실제 모델
- `HuffmanDecoder`: 선택적 추가 구현

최소 실험은 `None`, `Ideal`, `DECA` 세 가지로 구성한다.

### 4.3 Pipeline 연결

1. Compressed weight tile을 DRAM에서 Global Buffer로 전송한다.
2. Global Buffer의 compressed tile을 decoder input queue에 삽입한다.
3. Decoder startup latency와 throughput을 적용한다.
4. Decompression이 완료되면 dense tile을 Local Buffer로 전송한다.
5. Decoder queue가 가득 찬 경우 DRAM 또는 scheduler를 stall시킨다.
6. Dense tile이 준비되지 않은 경우 PE 실행을 대기시킨다.

### 4.4 추가 통계

- Compressed DRAM traffic
- Metadata traffic
- Decode busy cycles
- Decoder input/output queue stall cycles
- PE waiting-for-decoder cycles
- Dense bypass 횟수
- Decoder utilization
- Global/Local Buffer occupancy

---

## 5. 모델 검증 계획

### 5.1 Byte Accounting 검증

- Dense weight 크기 검증
- Compressed payload 크기 검증
- Metadata 및 scale 크기 검증
- Decoder 출력 크기가 원래 dense 크기와 동일한지 검증

### 5.2 Timing 검증

1. `CR=1`이고 decode latency가 0이면 dense baseline과 cycle이 동일해야 한다.
2. Decoder throughput이 무한대이면 Ideal compression 결과와 동일해야 한다.
3. DRAM bandwidth가 충분히 크면 compression speedup이 거의 사라져야 한다.
4. Overlap을 사용하지 않으면 다음 관계를 따라야 한다.

   $$
   T_{\text{total}}
   \approx T_{\text{DRAM}} + T_{\text{decode}}
   $$

5. Streaming overlap을 적용하면 다음 관계에 가까워져야 한다.

   $$
   T_{\text{total}}
   \approx L_{\text{decode}}
   + \max(T_{\text{DRAM}}, T_{\text{decode}})
   $$

6. Decoder throughput이 낮아질수록 PE waiting cycle이 증가해야 한다.

### 5.3 외부 결과 비교

가능하면 DECA 논문의 compressed GEMM 결과 또는 analytical model과 경향을 비교한다.

DECA의 구조를 완전히 재현하지 않는다면 논문에서는 `DECA-inspired decompression engine`으로 표현한다.

---

## 6. Workload 구성

전체 LLM을 처음부터 functional simulation하기보다 대표 GEMM layer를 먼저 평가한다.

### 대상 연산

- QKV projection
- Attention output projection
- FFN up projection
- FFN gate projection
- FFN down projection

### 실행 형태

| 실행 단계 | GEMM M dimension |
|---|---:|
| Decode | 1, 4, 8 |
| Small-batch inference | 16, 32 |
| Prefill | 128, 512 |

Llama 계열 7B 및 70B 모델의 representative layer dimension을 사용하는 것을 우선 고려한다.

---

## 7. Sensitivity 실험

### 7.1 실험 변수

| 변수 | 설정 |
|---|---|
| Compression ratio | 1×, 1.25×, 1.5×, 2×, 3×, 4× |
| Decoder throughput | PE 소비율 대비 0.5×, 1×, 2× |
| DRAM bandwidth | 0.5×, 1×, 2×, 4× |
| 실행 형태 | Decode, Prefill |
| Decoder model | None, Ideal, DECA |

### 7.2 실험 순서

#### Phase 1: Synthetic Sensitivity

Compression ratio를 독립적인 파라미터로 설정하여 전체적인 성능 경향과 break-even point를 분석한다.

이 결과는 architecture-level sensitivity 분석으로 사용하며, 실제 LLM compression 결과로 해석하지 않는다.

#### Phase 2: Realistic Compression Points

실제 quantized 또는 sparse LLM weight에서 다음 크기를 측정한다.

- Compressed values
- Index 또는 bitmap metadata
- Group quantization scale
- Padding 및 alignment overhead

측정된 compression ratio를 synthetic sensitivity 결과 위에 표시한다.

---

## 8. 평가 지표

### 주요 성능 지표

- End-to-end execution cycles
- Dense baseline 대비 speedup
- DRAM traffic reduction
- Decoder utilization
- PE utilization
- PE waiting-for-decoder cycles
- Decoder queue stall cycles

### 선택적 에너지 지표

- DRAM access energy
- Global Buffer access energy
- Decoder energy
- 전체 accelerator energy

Decoder energy 또는 area의 신뢰할 수 있는 모델을 확보하지 못한 경우, 초기 실험에서는 execution cycle과 traffic 분석에 집중한다.

---

## 9. 예상 결과

### 낮은 Compression Ratio

- Traffic 감소가 작음
- Metadata와 decoding overhead가 상대적으로 큼
- Dense baseline보다 성능이 낮아질 가능성 존재

### 중간 Compression Ratio

- DRAM traffic 감소가 decoder overhead를 상쇄
- Compression의 break-even point 발생
- Memory-bound workload에서 성능 향상

### 높은 Compression Ratio

- DRAM latency는 지속적으로 감소
- Decompression engine이 새로운 병목으로 전환
- Compression ratio 증가에 따른 speedup 포화

### Decode와 Prefill 비교

- Decode는 weight reuse가 적고 memory-bound하므로 compression 효과가 큼
- Prefill은 compute intensity와 weight reuse가 증가하여 compression 효과가 감소
- Batch 또는 GEMM의 M dimension이 증가할수록 decoder overhead를 숨기기 쉬워짐

### DRAM Bandwidth 영향

- 낮은 DRAM bandwidth에서는 compression 효과가 큼
- 높은 DRAM bandwidth에서는 compression 효과가 감소
- 충분한 bandwidth에서는 decoder overhead가 오히려 성능을 제한할 수 있음

---

## 10. 논문 Figure 구성

하나의 Figure를 세 개의 subplot으로 구성한다.

### Figure (a): Compression Ratio Sensitivity

- X축: Compression ratio
- Y축: Normalized execution cycles
- 비교: Dense, Ideal, DECA
- 패널: Decode와 Prefill

### Figure (b): Decoder Throughput Sensitivity

- X축: Compression ratio
- Y축: Decoder throughput
- 색상: Dense baseline 대비 speedup

### Figure (c): Execution Cycle Breakdown

- DRAM transfer
- Decompression
- On-chip data movement
- PE computation
- Stall

대표 compression ratio인 `1×`, `2×`, `4×`를 비교한다.

---

## 11. Configurability 평가

성능 결과와 함께 다음 구현 지표를 보고한다.

- Decompression engine 추가 LOC
- 새로 추가한 파일 수
- 수정한 기존 파일 수
- 기존 PE model 수정 여부
- 기존 scheduler 수정 여부
- Configuration만으로 decoder를 활성화하거나 교체할 수 있는지
- 기존 accelerator regression test 통과 여부

### 목표

- 기존 PE와 scheduler의 핵심 구현은 변경하지 않는다.
- Decompression engine을 독립적인 component로 추가한다.
- Configuration을 통해 `None`, `Ideal`, `DECA`를 선택할 수 있도록 한다.

---

## 12. 논문 서술 방향

### Motivation

> As LLMs continue to grow in size, weight compression has become an important technique for reducing memory capacity and bandwidth requirements. However, compressed weights must be reconstructed before computation, and the decompression overhead can become a new performance bottleneck.

### Configurability Case Study

> To demonstrate NPUsim's configurability, we integrate a DECA-inspired weight-decompression engine as an optional hardware component between the compressed weight buffer and the PE array.

### Evaluation Question

> Using the extended model, we investigate how compression ratio, decompression throughput, memory bandwidth, and workload characteristics jointly affect LLM inference performance.

### 예상 결론

> The results reveal that increasing the compression ratio does not always translate into proportional performance gains. Once off-chip memory pressure is sufficiently reduced, decompression throughput becomes the primary bottleneck, resulting in diminishing performance returns.

---

## 13. 진행 일정

| 기간 | 작업 |
|---|---|
| 1–2일 | 모델 범위 및 decompression parameter 확정 |
| 3–5일 | DecompressionEngine 및 queue 구현 |
| 6–7일 | Microbenchmark와 byte/timing 검증 |
| 8–10일 | LLM representative layer 실험 |
| 11–12일 | Sensitivity 분석 및 이상 결과 점검 |
| 13–14일 | Figure 생성 및 논문 문장 작성 |

---

## 14. 최종 산출물

- DECA-inspired decompression component
- Decoder configuration file
- Byte 및 timing unit test
- LLM representative layer configuration
- Compression-ratio sensitivity 결과
- Decoder-throughput sensitivity 결과
- Cycle breakdown 결과
- Configurability 구현 지표
- 논문용 Figure와 case-study 설명
