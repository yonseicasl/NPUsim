# NPUsim 현재 버전: 기능 및 엔지니어링 완성도 평가

## 평가 범위와 결론

이 문서는 저장소의 정적 코드 검토와 설정·빌드 경로 점검에 기반한다. 외부 의존성(Nebula, DRAMSim3, OpenCV 등)이 현재 작업공간에 없어 전체 simulator 빌드, end-to-end golden 결과, 전체 simulator sanitizer 및 calibration은 수행하지 못했다. 다만 아래 결함 수정 후 dependency-free regression과 config/mapping parser ASan·UBSan은 수행했다.

NPUsim은 PE, PE array, global buffer, multi-chip, DRAM을 분리하고 mapping table로 데이터플로를 지정하는 **연구용 CNN accelerator timing/energy simulator**로서 좋은 구조를 갖고 있다. 최초 평가에서 확인한 multi-chip, functional buffer, configuration validation 결함은 아래와 같이 수정했다. 다만 convolution/FC 외 레이어 비용, end-to-end golden 결과와 calibration이 없으므로 범용·검증된 full-model simulator로 볼 수 없다는 범위 판단은 유지한다.

## 구현 조치 현황 (2026-08-07)

| 우선순위·진단 | 조치 | 검증 상태 |
| --- | --- | --- |
| P0 multi-chip 구성 불일치 | `[accelerator] num_chips`를 component 생성 전에 읽고, 물리 chip 수 및 active mapping 상한 검증 추가 | 제공 accelerator 7개 정적 계약 통과; end-to-end multi-chip은 외부 의존성 부재로 미실행 |
| P0 PE-array temporal buffer | 공통 helper로 통합, `separate/shared` 매핑과 legacy `*_size` fallback 수정, `num_pes` 한 번만 적용 | 세 PE-array 구현 및 제공 설정의 계산 검증 |
| P1 numeric initialization | `memset(..., 1.0, ...)` 제거, element 수 ceiling 및 `std::fill_n(data_t{})` 적용 | default·user integer·user float functional 구문/초기화 검증 |
| P1 수명 관리 | scheduler clear-before-delete 제거, NPU 소유 객체 전체 해제, 포인터 초기화 | 구문 및 정적 소유권 검토 통과 |
| P1 config/mapping validation | section/key/숫자/vector/mapping 11항목, 필수 component·개수·active capacity 검증 추가 | 182개 설정과 8개 음성 fixture, parser ASan·UBSan 통과 |
| full-model 불완전 | 빈 pooling API 제거; PE는 예상 밖 layer를 명시적으로 거부; NPU는 비지원 layer를 경고·제외하고 total을 partial로 표시 | PE 3개 data type mode 구문 검증; 비용 모델 구현 자체는 남음 |
| 설정 실행 경로 | launcher와 model을 루트 `configs/` 기준으로 통일, 파일 존재 및 안전한 identifier 검사 | shell syntax·path rejection 검증 |
| 중복·dead code | temporal-buffer helper 통합, 미사용 mapping overload·thread 선언·systolic workspace·offset TODO 제거 | `rg` 기반 재발 검사 및 관련 translation unit 구문 검증 |
| 추가: PE local buffer key 불일치 | PE도 *_size/*_buffer를 호환 처리하고 non-zero를 강제 | 세 functional data type 전체 변경 translation unit 구문 검증 |
| 추가: systolic dead flatten workspace | null write 또는 초대형 할당 후 결과를 사용하지 않던 flatten/workspace 경로 제거 | 세 functional data type systolic 구문 검증 |
| 추가: mapping zero·overflow 및 USER include 순서 | zero stride·누적 overflow 거부, energy2.map 136개 stride 수정, scheduler의 data.h 직접 include | 8개 음성 fixture 및 integer/float 전체 변경 translation unit 검증 |
| build 재현성 | option 변수/문법·Makefile 오타 수정, Nebula/DRAMSim3 revision 고정, PyTorch 자동 clone 제거 | `bash -n`, `git diff --check`, validation suite 통과 |
| 자동 검증 | `tests/run_validation.sh`와 GitHub Actions 추가 | 일반 및 parser ASan·UBSan 로컬 통과 |
| golden/reference/calibration | 구현하지 않음 | 외부 실행 환경과 기준 결과가 필요 |

> 아래 결함 설명은 최초 원인 분석을 보존하기 위한 수정 전 서술이다. 현재 조치 여부는 위 표를 우선한다.

### 추가 진단 상세

- **PE local buffer key 계약:** PE-array helper만 input_buffer/weight_buffer/output_buffer를 읽고 개별 PE는 input_size/weight_size/output_size만 읽었다. TPU, MAERI, Simba 등에서는 local buffer가 0 element로 할당될 수 있었다.
- **systolic functional workspace:** workspace_input은 기본 convolution 경로에서 null인 채 flatten() 대상으로 전달될 수 있었다. GEMM 설정에서는 byte 크기 두 개를 곱해 매우 큰 배열을 만들었지만, flatten 결과는 실제 transfer source로 사용되지 않았다.
- **mapping 산술 계약:** mapping 값 자체는 unsigned 범위여도 component 누적 곱에서 overflow할 수 있었고 GROUP/STRIDE=0을 허용했다. 실제 simba/mobilenetv3/energy2.map에서 주석은 stride 1인데 136개 component 행이 stride 0이었다.
- **USER data type include 순서:** scheduler.h가 macro 정의만 포함하고 data_t 정의는 직접 포함하지 않아, include 순서에 따라 USER_INTEGER/USER_FLOAT 빌드가 실패했다.

위 항목은 모두 수정됐으며 일반 regression, parser ASan·UBSan, 세 functional data type의 경고-as-error 구문 검증을 통과했다.

| 관점 | 평가 | 판단 |
| --- | ---: | --- |
| 아키텍처 모델링·DSE 구조 | 8/10 | 계층형 컴포넌트와 mapping 분리가 강점 |
| 단일 칩 timing 연구 코드 | 5/10 | 제한된 설정에서만 사용 권장 |
| functional simulation 신뢰도 | 2/10 | 버퍼 초기화·설정 불일치 문제 존재 |
| multi-chip simulation 신뢰도 | 1/10 | 객체 수와 active-chip 접근이 정합하지 않음 |
| 공개 프레임워크/소프트웨어 완성도 | 3/10 | 검증, 재현성, 오류 처리 부족 |

## 이전 기능 평가에서 바로잡아야 할 점

### 실제 실행 지원 범위는 설정 자산보다 좁다

- 루트 `configs/`에는 Eyeriss 계열, MAERI, Simba, TPU 계열과 BERT/GPT/ViT용 mapping을 포함한 많은 설정 파일이 있다.
- 그러나 실행 파일은 `models/accelerators`, `models/networks`, `models/mappings`만 탐색한다. 현재 이 경로에는 Eyeriss와 AlexNet/ResNet-50 중심의 일부 자산만 있다.
- 실행 루프에서 NPU timing simulation 대상은 convolution과 fully-connected뿐이다. max/avg pooling은 조건문에서 주석 처리되어 있고, shortcut, concat, softmax, upsample, excitation은 accelerator timing에 반영되지 않는다.
- 따라서 BERT/GPT/ViT mapping 파일의 존재만으로 Transformer end-to-end 지원을 주장할 수 없다. 네트워크 로더도 CNN 중심 `nebula::convolutional_t`를 생성한다.

### 기본 모드는 functional/DRAM-accurate 검증 모드가 아니다

- 빌드 스크립트의 기본값은 `FUNCTIONAL=0`, `DRAMSIM3=0`이다.
- 기본 실행 결과는 실제 tensor 값의 정확도를 검증한 결과라기보다, 설정 기반 timing/energy 추정 결과로 해석해야 한다.
- unit test, regression test, CI, golden output, 논문/실측 수치와의 calibration report가 저장소에 없다.

## 기능 및 구조 측면의 장점

- **아키텍처 분리:** PE, spatial/adder-tree/systolic PE array, global buffer, multi-chip, DRAM이 별도 컴포넌트다.
- **cycle-driven 진행:** 요청, 전송, 연산을 단계적으로 진행하며 compute와 memory hierarchy의 상호작용을 모델링한다.
- **세부 통계:** layer/network 단위로 computation, access, transfer, static energy와 cycle을 수집한다.
- **mapping 기반 DSE:** MAC부터 DRAM까지 매핑을 지정하고 input/weight/output-stationary dataflow를 지원한다.
- **확장 지점:** dense 및 여러 sparse format, functional mode, DRAMSim3 연동을 위한 코드 경로가 존재한다.

## 코드 리뷰: 확인된 결함

### P0 — multi-chip 구성과 실제 객체 생성 수가 다르다

Accelerator config는 `num_chips`를 사용하지만, `npu_t::init()`은 존재하지 않는 `num_processors`와 `num_pes`만 읽는다. 따라서 `num_processors`는 기본값 1로 남고, PE array와 global buffer도 하나만 생성된다.

반면 실행 루프는 mapping으로부터 계산한 active chip 수까지 `pe_arrays[i]`와 `global_buffers[i]`에 접근한다. 제공 mapping을 정적으로 집계하면 active `CHIPS_X=64`, active `CHIPS_Y=36`인 사례가 있다. 다중 칩 mapping에서는 범위 밖 접근, 크래시, 또는 잘못된 결과가 발생할 수 있다.

영향: Simba, Eyeriss-v2, TPUv3 등 multi-chip accelerator 결과는 수정 전 신뢰할 수 없다.

### P0 — PE-array temporal buffer 설정이 잘못 해석된다

`spatial_arch`, `adder_tree`, `systolic_array`의 공통적인 버퍼 초기화 코드에는 다음 문제가 있다.

- `memory_type=shared`일 때 내부 enum을 `SEPARATE`로, `separate`일 때 `SHARED`로 반대로 기록한다.
- 현재 `models/accelerators/eyeriss.cfg`는 `input_size/weight_size/output_size`를 쓰지만 구현은 `input_buffer/weight_buffer/output_buffer`만 읽는다. 이 설정에서는 PE-array 버퍼 크기가 0으로 남는다.
- 크기를 `num_pes`배 한 뒤 allocation 수에서 다시 `num_pes`배 한다. PE당 버퍼 크기를 의도한 구현이라면 PE 수만큼 과다 할당이다.
- reset은 allocation된 전체 메모리가 아니라 `input_size/weight_size/output_size`만 초기화한다.

영향: functional path의 데이터 주소와 값, 버퍼 용량, access model이 서로 불일치할 수 있다.

### P1 — `memset(..., 1.0, ...)`는 숫자 1 초기화가 아니다

`memset`의 값 인자는 byte 값이다. `1.0`은 `1`로 변환되어 각 byte를 `0x01`로 채우며, float/integer `data_t`의 값 `1`을 만들지 않는다. PE local buffer, MAC register, PE-array buffer, global buffer에서 반복된다.

영향: functional simulation의 초기 tensor 값이 의도와 달라진다. `std::fill_n` 또는 실제 데이터 로더를 사용해야 한다.

### P1 — 수명 관리 오류

`npu_t` destructor는 `schedulers.clear()`를 먼저 호출한 후 scheduler를 delete하려고 한다. clear 이후 크기는 0이므로 scheduler 객체가 해제되지 않는다.

영향: 프로세스 종료 시 memory leak이며, 반복 시뮬레이션 또는 라이브러리 임베딩 시 문제가 커진다.

### P1 — configuration/mapping 검증이 없다

설정 parser와 mapping parser는 section 형식, key의 중복, section 이전의 key-value, 숫자 변환 실패, vector 길이, component 수, mapping 수와 network layer 수의 정합성을 검증하지 않는다. 오류 설정이 즉시 실패하지 않고 기본값 또는 불완전한 값으로 실행될 수 있다.

특히 vector parser는 입력 길이가 부족해도 성공을 반환하고, mapping parser는 11개 항목이 있다고 가정한다.

### P1 — 필수 accelerator component 계약을 검증하지 않는다

`npu_t::init()`은 설정 section을 순회하며 PE array, global buffer, multi-chip, DRAM을 선택적으로 생성한다. 그러나 초기화 직후 `pe_arrays[0]`, `multi_chip`, `dram`을 조건 없이 역참조하고 `connect()`에서도 같은 객체들을 사용한다. 즉 필수 section이 누락되었거나 section 이름이 잘못된 config는 명시적인 validation 오류 대신 null dereference 또는 vector 범위 밖 접근으로 실패할 수 있다.

영향: 설정 오류의 원인을 사용자가 진단하기 어렵고, 새로운 accelerator config를 추가할 때 실패 방식이 비결정적이다. 초기화 전에 필수 component의 정확히 한 개 존재 여부, PE array/global buffer 수의 일치, `num_processors`와 `num_chips`의 계약을 검사해야 한다.

### P1 — full-model 연산 구현이 불완전하다

- `pe_t::max_pooling()`과 `pe_t::avg_pooling()`은 빈 함수다.
- NPU run loop는 convolution/connected 외 레이어를 timing simulation하지 않는다.
- systolic array의 tile-offset update에는 TODO가 남아 있다.

영향: pooling, residual graph, concat, activation 등 비-convolution 연산이 포함된 모델의 layerwise/total timing은 전체 모델 비용이 아니다.

## 코드 리뷰: 중복 및 유지보수 위험

- `spatial_arch`, `adder_tree`, `systolic_array`의 temporal-buffer 초기화가 대량 복제되어 있다. 실제로 동일한 memory-type 반전과 `memset` 오류가 세 구현 모두에 반복된다. base class helper로 통합해야 한다.
- `execute_thread()`는 선언만 있고 구현·호출이 없다. config의 `num_threads` 역시 코드에서 사용되지 않아 실제 병렬 실행 기능은 제공되지 않는다.
- parameter-order를 받는 `mapping_table_t::calculate_parameter_size` overload는 호출되지 않는다. 고정 길이 배열 `strcpy`, 잘못된 대소문자 조건(`P || P`) 등 잠재 결함이 있어 제거 또는 테스트 후 재구현이 필요하다.
- `systolic_array_t`의 GEMM workspace는 input만 할당·사용하며 weight/output workspace는 선언·해제만 된다.
- 빌드 스크립트는 `USER_INTEGER`를 정의하고 `USE_INTEGER`를 검사하며, float option에 쉘 문법 오류가 있다.
- library Makefile의 `schduler` 오타로 scheduler header 변경이 dependency에 잡히지 않는다.

## 권장 수정 순서

1. **multi-chip 정합성 복구:** config의 `num_chips` 또는 multi-chip dimensions를 기준으로 PE array/global buffer를 생성하고, active-chip 수가 vector 크기를 넘지 않도록 검증한다.
2. **버퍼 모델 단일화:** 세 PE-array 구현에서 공통 buffer allocation helper를 만들고, config key/단위/ownership을 하나로 정의한다. memory type enum 반전과 이중 `num_pes` 곱셈을 수정한다.
3. **functional correctness 복구:** 모든 `memset(..., 1.0, ...)`을 `std::fill`로 대체하고, PyTorch/Nebula reference와 layerwise output을 비교한다.
4. **입력 validation 추가:** 필수 component의 존재·개수와 연결 가능 여부를 먼저 검사한 뒤, 필수 key, enum 값, vector 길이, mapping section 수, active component 수와 physical capacity를 검증하고 오류를 명시적으로 반환한다.
5. **기능 범위 명확화 또는 확장:** pooling/elementwise/shortcut/concat 등 full-model layer의 timing 및 memory traffic을 구현하거나, 공식 지원 범위를 convolution/FC-only로 제한한다.
6. **검증 기반 구축:** 대표 단일 칩 예제부터 build/run/golden-result regression test를 만들고, AddressSanitizer/UndefinedBehaviorSanitizer와 CI를 추가한다.
7. **정리:** RAII(`std::vector`, `std::unique_ptr`)로 소유권을 정리하고, dead API와 미사용 코드, 중복 구현을 제거한다.

## 사용 권고

수정 전에는 단일 칩 Eyeriss 계열의 제한된 timing/energy 연구에만 사용하고, 결과에는 지원 레이어와 제외된 비용을 명시해야 한다. multi-chip 결과, sparse functional 결과, 전체 CNN/Transformer 정확도 결과는 regression 및 sanitizer 검증 전까지 논문의 정량 근거로 사용하지 않는 것이 안전하다. README의 full-model/functional 표현도 이 제한된 지원 범위와 검증 상태를 반영하도록 수정해야 한다.
