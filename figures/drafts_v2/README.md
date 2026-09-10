# NPUsim — Fig. 2–4 초안

작성 기준: 2026-09-10 작업 폴더의 `abstract.tex`, `1_introduction.tex`, `2_background.tex`, `3_related_work.tex`, `4_npusim.tex`의 활성 본문 전체.
스타일 참고: `paper_reference/DATE2027__NPUsim_.pdf`의 Fig. 2–4와 `paper_reference/DATE2025__RoTA____Accepted.pdf`의 구조 도식. 흰 배경, 검은 선, 점선 그룹, 회색 음영, 직교 연결선, sans-serif 글꼴을 사용했다.

`index.html`을 브라우저에서 열면 세 그림, 설계 의도, 영문 캡션을 함께 볼 수 있다. ‘논문 폭’은 7인치/3.5인치 CSS 폭을 보여준다. 화면상의 물리적 크기는 OS/브라우저 배율에 따라 다를 수 있다.

| 그림 | 편집 원본 | 권장 논문 폭 | 핵심 내용 |
|---|---|---|---|
| Fig. 2 | `fig2_npusim_overview.svg` | 두 단, 7 in | 분리된 입력, 실행 스케줄러, 하드웨어 모델, 비용 및 결과 |
| Fig. 3 | `fig3_configurable_architecture.svg` | 한 단, 3.5 in | 설정 블록 → 공통 모델 라이브러리 → 세 아키텍처 예시 |
| Fig. 4 | `fig4_decoupled_scheduling.svg` | 두 단, 7 in | (a) 타일 스케줄과 요청 전파, (b) 버퍼링·지연·에너지 |

동일한 이름의 PDF는 벡터 출력이고, PNG는 2배 해상도의 미리보기다. Fig. 4의 패널은 `fig4a_execution_scheduling.*`, `fig4b_timing_energy.*`로도 제공한다. 원본 LaTeX와 기존 그림은 변경하지 않았다.

## 그림에 반영한 판단

- Fig. 2: hardware specification과 scheduling scheme을 서로 독립적인 입력으로 나타냈다. PyTorch의 layer dimensions는 스케줄러로, optional tensor values는 functional execution으로 연결했다. 비용 모델은 타이밍과 활동별 비용 계산에 사용된다.
- Fig. 3: configuration block의 의미를 나타내는 개념도이며, 실제 설정 파일 문법을 주장하지 않는다. Eyeriss-like / TPU-like / Simba-like는 본문의 구성 예시다. PE와 chip 개수는 도식을 단순화한 것으로 실제 제품의 개수/치수를 뜻하지 않는다. chip 내부의 작은 직사각형들은 global buffer와 PEs를 나타낸다. 이미 구현된 모델의 설정 변경과 초기화 시 조립을 표현한다.
- Fig. 4(a): symbolic offsets `w₀`, `i₀`, `o₀`를 사용했다. WS 예시에서 weight tile은 네 번 재사용되고 input/output tile은 교체된다. 이는 설명용 매핑의 예시이며 모든 WS 매핑에 보편적인 input/output reuse count를 주장하지 않는다. 각 경계는 서로 다른 stationary dataflow를 가질 수 있다. GB = global buffer, LB = local buffer, RF = register file.
- Fig. 4(b): 실험 결과가 아닌 **설명용 타이밍**이다. 네 work tile의 DRAM→GB / GB→LB / MAC 비용을 2 / 2 / 4 cycles로 가정한다. 두 경계에 각각 2개 slot이 있고 consumer가 끝나야 slot이 해제된다. GB→LB는 [6,8), [10,12)에서 LB 점유로 대기하고 compute는 [4,20)에서 연속 실행한다. (a)의 symbolic tensor offsets와 (b)의 T0–T3는 서로 다른 표기이며, 후자는 파이프라인 work tile을 뜻한다. 공유 port 충돌 등의 추가 제약이 없는 예시다.
- Fig. 4는 현재 Section 4가 execution scheduler와 timing subfigure를 모두 참조하므로 두 패널을 유지했다. 가독성을 위해 두 단 폭을 권장한다.

## 영문 캡션 초안

**Fig. 2.** Overview of NPUsim. Separate hardware specifications, scheduling schemes, and DNN workloads drive modular hardware models and execution schedulers. Schedulers determine tile transfers and execution order, while hardware models determine cycle costs and component availability. Traced activities and component cost models produce runtime statistics; functional execution optionally evaluates actual tensor values.

**Fig. 3.** Configurable architecture modeling in NPUsim. Hardware configuration blocks select and parameterize a shared component library, which is instantiated and connected at initialization to compose spatial (Eyeriss-like), systolic (TPU-like), and multi-chip (Simba-like) architectures. Optional configuration blocks enable specialized optimization modules. The depicted array and chip counts are illustrative.

**Fig. 4.** Decoupled scheduling and cycle-level execution in NPUsim. (a) Per-boundary schedulers derive tensor tile offsets and reuse counts from the commanded dataflow and mapping. A weight-stationary example holds one weight tile across four input tiles; runtime requests propagate up the memory hierarchy as needed. (b) Hardware costs and buffering compose scheduled work into a cycle-level pipeline. In this illustrative double-buffered example, the compute stage is the bottleneck and upstream transfers stall while local buffer slots remain occupied. The same activities determine energy consumption.

## 재생성

```bash
python3 figures/drafts_v2/generate_figures.py
```

시스템의 `librsvg`, `libcairo`, `libgobject`와 Python 표준 라이브러리를 사용한다. 외부 이미지, 폰트 다운로드, 네트워크 의존성은 없다. SVG의 텍스트와 도형은 개별 편집 가능하다. 브라우저/편집기에서는 Arial / Helvetica / Liberation Sans 순으로 글꼴을 선택한다. PDF는 생성 환경에서 사용한 글꼴을 내장한다.

LaTeX 삽입 예시는 `latex_snippets.tex`에 있다. 현재 본문에 이미 있는 figure 환경을 대체할 때 사용하며 중복 삽입하지 않는다.
