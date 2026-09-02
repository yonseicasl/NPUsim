# NPUsim 현재 버전 평가

## 요약

NPUsim은 CNN 기반 DNN accelerator의 설계공간 탐색을 위한 연구용 cycle-level simulator로서, 아키텍처 계층과 데이터플로를 분리한 구조가 강점이다. 다만 현 버전은 기능 범위, 자동 검증, 의존성 재현성 측면에서 아직 연구 프로토타입 성격이 강하다.

> **평가 근거와 한계:** 외부 의존성이 작업공간에 없어 전체 simulator 빌드, end-to-end golden 결과, 전체 simulator sanitizer 및 calibration은 아직 수행하지 못했다. 대신 dependency-free regression과 config/mapping parser ASan·UBSan 검증을 추가·수행했다.

## 수정 반영 상태 (2026-08-07)

| 초기 진단 | 현재 상태 | 검증 |
| --- | --- | --- |
| multi-chip 객체 수·active chip 불일치 | 수정 완료 | 7개 accelerator의 `num_chips == height*width`, active component 상한, `npu.cc` 구문 검증 |
| PE-array temporal buffer 해석·과다 할당 | 수정 완료 | 세 PE-array 구현을 공통 helper로 통합하고 제공 설정의 byte 수·할당 수 검증 |
| 잘못된 numeric `memset`과 functional 분기 오류 | 수정 완료 | default/`USER_INTEGER`/`USER_FLOAT` functional 구문 및 zero initialization 검증 |
| scheduler 및 component 수명 관리 | 수정 완료 | destructor 순서·전체 owned object 해제 코드 검토와 구문 검증 |
| config/mapping/필수 component 검증 부재 | 수정 완료 | 182개 설정 자산 통과, malformed·duplicate·short/long/non-numeric 음성 테스트 통과 |
| `models/`와 루트 `configs/` 경로 불일치 | 수정 완료 | launcher가 루트 `configs/`만 사용하며 누락 파일을 초기화 전에 거부 |
| build option·Makefile·dead API·중복 초기화 | 수정 완료 | shell syntax, 세 functional type, 정적 회귀 검사 통과 |
| PE local buffer key 불일치·systolic dead workspace | 수정 완료 | 0-byte/null/OOM 경로 제거, 세 functional data type의 전체 변경 translation unit 구문 검증 |
| mapping zero/overflow 허용 | 수정 완료 | zero stride·누적 overflow 거부 fixture와 energy2.map stride 오타 136개 교정 |
| USER mode include-order 의존 | 수정 완료 | scheduler가 data.h를 직접 포함하며 integer/float 전체 변경 translation unit 검증 통과 |
| 의존성 버전 미고정 | 수정 완료 | Nebula·DRAMSim3 commit 고정 및 override 지원, PyTorch source clone 제거 |
| unit/regression/CI 부재 | 부분 해소 | dependency-free regression과 GitHub Actions parser sanitizer 추가 |
| full-model layer timing 부재 | 범위 제한으로 처리 | convolution/FC 이외 레이어를 경고·제외하고 total을 partial로 표시, README에 제한 명시 |
| golden output·reference 비교·calibration 부재 | 미완료 | 외부 의존성과 기준 데이터가 필요 |

> 아래 장단점과 우선 과제는 최초 진단의 근거를 보존한다. 이미 수정된 항목의 현재 판단은 위 상태표가 우선하며, full-model 확장·golden regression·calibration은 남은 과제다.


| 관점 | 평가 |
| --- | ---: |
| 연구용 프로토타입 완성도 | 7/10 |
| 공개 도구·범용 프레임워크 완성도 | 4/10 |

## 장점

- **명시적 아키텍처 모델링:** PE, PE array, global buffer, multi-chip, DRAM을 독립 컴포넌트로 구성한다. loop-centric 모델보다 하드웨어 구조와 데이터플로 변경을 분리해 DSE에 적합하다.
- **cycle-level 실행 흐름:** 연산, 데이터 전송, 요청을 cycle 단위로 진행해 compute와 memory hierarchy 간 상호작용을 모델링한다.
- **풍부한 비용 분석:** 계층별 및 네트워크 전체에 대해 computation/access/transfer/static energy와 cycle 통계를 산출한다.
- **구성 가능한 데이터플로:** MAC부터 DRAM까지의 mapping table과 input/weight/output-stationary 방식을 지원한다.
- **아키텍처 자산:** Eyeriss, Eyeriss v2, MAERI, Simba, TPU, TPU v3 등 여러 accelerator 설정과 다수의 mapping 파일이 있다.
- **희소성·기능 시뮬레이션 설계:** dense와 여러 sparse format, 실제 데이터 기반 functional simulation, DRAMSim3 연동을 위한 코드 경로가 준비되어 있다.

## 단점 및 제약

- **실제 timing 실행 레이어가 제한적:** `npu_t::run()`은 convolution과 fully-connected 레이어만 accelerator에서 실행한다. pooling 처리는 주석 처리되어 있고 shortcut, concat, softmax, upsample, excitation 등의 비용·의존성은 전체 timing에 반영되지 않는다.
- **Transformer 지원은 실험적 수준:** BERT/GPT/ViT용 mapping 파일은 있으나, 네트워크 로더가 CNN 중심이며 end-to-end Transformer 실행을 보장하는 구현·예제가 부족하다.
- **설정 자산과 실행 경로의 불일치:** 루트 `configs/`에는 다양한 accelerator/network/mapping 설정이 있지만, 실행 프로그램은 `models/` 아래 설정만 탐색한다. 현재 `models/`에는 일부 Eyeriss, AlexNet, ResNet-50 자산만 있다.
- **multi-chip 결과는 현재 신뢰할 수 없음:** accelerator config의 `num_chips`와 객체 생성에 쓰는 `num_processors`가 정합하지 않다. 그 결과 PE array/global buffer는 하나만 생성될 수 있지만 실행 루프는 active chip 수만큼 접근한다. multi-chip mapping은 범위 밖 접근·크래시·오답의 위험이 있다. 세부 근거는 `NPUSIM_ENGINEERING_ASSESSMENT.md`의 P0 항목을 따른다.
- **기능 시뮬레이션과 DRAMSim3가 기본 비활성:** 기본 빌드는 timing/energy estimator 성격이 강하며, 실제 값 기반 정확도 검증 경로는 별도 빌드 옵션을 켜야 한다.
- **재현성 부족:** Nebula, DRAMSim3, PyTorch, OpenCV 등의 외부 의존성을 실행 스크립트에서 동적으로 처리하고 버전을 고정하지 않는다. 현 작업공간에도 `ext/`가 없어 end-to-end 실행을 바로 검증할 수 없다.
- **자동 검증 부재:** unit test, regression test, CI, golden output, 논문/실측 결과와의 calibration 문서가 없다. cycle/energy 모델의 정확도를 지속해서 보장하기 어렵다.
- **구현 정리 필요:** 빌드 스크립트의 `USER_INTEGER`/`USE_INTEGER` 변수 불일치, float 옵션 처리 문제, scheduler 객체 해제 순서 오류 등 유지보수 위험 신호가 있다.
- **공개 문서와 실제 지원 범위의 차이:** README는 full-model·functional simulator로 소개하지만, 현재 실행 경로는 convolution/fully-connected 중심이고 functional mode에도 초기화·버퍼 설정 결함이 있다. 현 상태의 공식 지원 범위는 “검증 전 단일 칩 CNN timing/energy 연구 코드”로 제한해 표기하는 것이 정확하다.

## 우선 개선 과제

1. multi-chip 객체 수와 active-chip 수의 정합성을 우선 복구하고, 범위를 벗어난 접근을 검증 단계에서 차단한다.
2. `configs/`와 `models/` 설정 경로를 통합하고, 공식적으로 지원하는 accelerator/network/mapping 조합을 명시한다.
3. 의존성 버전 고정과 Dockerfile 또는 재현 가능한 설치 스크립트를 제공한다.
4. Eyeriss + AlexNet 등 대표 조합의 build/run/golden-result 비교를 CI regression test로 만든다.
5. pooling, residual/shortcut, concat, elementwise, activation의 latency와 memory traffic 모델을 추가한다.
6. functional mode 결과를 PyTorch 또는 Nebula 기준 결과와 계층별로 비교하는 정확도 검증을 추가한다.
7. 대표 accelerator의 latency/energy를 논문 또는 실측 수치와 비교한 calibration 보고서를 제공한다.

## 결론

현재 NPUsim은 **하드웨어·데이터플로 탐색을 위한 강한 연구 코드**다. 다만 범용 simulator로 신뢰성 있게 사용하려면 레이어 커버리지, 실행 가능 구성, 자동 검증, 재현성 패키징을 보완해야 한다.
