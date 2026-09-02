# Simulation issue 분류

이 디렉터리는 이슈를 검증 대상에 따라 분리한다. 같은 구현 결함이 두 결과에 영향을 줄 수 있지만, 각 문서에는 해당 검증 계약에 직접 필요한 문제와 완료 조건만 둔다.

- **Functional simulation**: 값, tensor 좌표, 데이터 이동 및 reduction 결과가 reference와 일치하는가.
- **Timing simulation**: cycle, traffic, contention, utilization, dynamic/static energy가 정의된 하드웨어 모델과 일치하는가.

| 구성요소 | Functional 정확성 | Timing/energy 정확성 |
|---|---|---|
| Runtime datatype | [functional/runtime_datatype.md](functional/runtime_datatype.md) | [timing/runtime_datatype.md](timing/runtime_datatype.md) |
| PE (MAC) | [functional/pe.md](functional/pe.md) | [timing/pe.md](timing/pe.md) |
| Local buffer | — | [timing/local_buffer.md](timing/local_buffer.md) |
| PE array | — | [timing/pe_array.md](timing/pe_array.md) |
| Spatial architecture | [functional/spatial_architecture.md](functional/spatial_architecture.md) | [timing/spatial_architecture.md](timing/spatial_architecture.md) |
| Systolic array | [functional/systolic_array.md](functional/systolic_array.md) | [timing/systolic_array.md](timing/systolic_array.md) |
| Adder tree | [functional/adder_tree.md](functional/adder_tree.md) | [timing/adder_tree.md](timing/adder_tree.md) |
| PyTorch frontend / LLM serving | — | [frontend/pytorch_vllm.md](frontend/pytorch_vllm.md) |
| Global buffer | — | [timing/global_buffer.md](timing/global_buffer.md) |
| Multi-chip / NoP | — | [timing/multi_chip.md](timing/multi_chip.md) |
| DRAM | — | [timing/dram.md](timing/dram.md) |
| Static energy (횡단) | — | [timing/static_energy.md](timing/static_energy.md) |
| Global cycle / compute-memory overlap (횡단) | — | [timing/global_cycle_overlap.md](timing/global_cycle_overlap.md) |

PyTorch frontend와 LLM serving 연동은 functional/timing 양쪽의 입력 계약을 바꾸는 교차 이슈이므로 [`frontend/pytorch_vllm.md`](frontend/pytorch_vllm.md)에서 별도로 관리한다. 이 문서는 framework graph·serving trace의 보존과 simulator lowering 범위를 정의하며, core component 이슈를 대체하지 않는다.

공통 설정 파싱·상태 초기화 결함은 functional 결과를 무효화하면 `functional/`에, 결과는 유지하지만 비용 추정을 왜곡하면 `timing/`에 기록한다. 둘 모두에 영향을 주는 경우에는 각 문서에서 자기 관점의 영향과 검증을 명시한다.
