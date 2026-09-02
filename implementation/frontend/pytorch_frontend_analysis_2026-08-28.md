# 2026-08-28 PyTorch frontend 전환 분석 + KN9 정리 + conv fixture 게이트

## 분석 결론

전환(커밋 b147068)은 **잘 이루어졌다**. 자체 테스트를 넘어 독립 검증으로 확인한 것:

1. **Nebula parity (결정적 증거)** — linear_relu fixture는 곧 RTL-검증된
   gemm_64x64x64다. 같은 mapping으로 두 경로를 A/B하면 compute-schedule
   3,518 cycle과 모든 DRAM traffic이 **비트 일치**, 차이는 정확히
   {frontend 출처 주석, SFU 블록, critical path +256 = SFU busy}뿐.
2. **conv2d e2e** (fixture 부재로 직접 제작) — MAC 7,776 = 6·6·8·3·3·3 정확,
   C=64 전체 커버리지에서 165,888 정확, 미달 매핑엔 명확한 경고.
3. **nebula loader 벽 구조적 회피** — IR 모드는 `load_data`를 건너뛴다
   (`if(!executable_ir_mode)`): 채널 1/3 제약이 IR 경로에는 없다.
4. **실패 경로 전부 시끄럽다** — SFU 없는 activation, 잘못된/중복/혼합 op_id,
   FUNCTIONAL 빌드, 기하-형상 불일치(15개 fail 지점) 모두 구체적 에러로 중단.
5. **회귀 무손상** — 외부 baseline 비트 불변, NVDLA holdout 불변.

유의점: ① torch 미설치 환경이라 실제 `torch.export` 캡처 경로는 이 환경에서
미검증(합성 graph 경유만). ② b147068이 당시 tree의 **다른 작업(NVDLA/timing)
미커밋 변경을 함께 커밋**해 이력이 섞였다. ③ IR은 fp32를 기록하지만 런타임
datatype은 config 소유 — report에 명시되는 의도된 계약.

## 조치 1 — KN9 4건 분류 (knobs 게이트 그린 복귀)

| knob | 분류 | 근거 |
| --- | --- | --- |
| `op_id` | OPTIONAL(도메인 밖) | mapping-파일 키 — 이 스캔의 declared-set(accelerator config)에 구조적으로 들어올 수 없음. binding 검증은 npu_frontend.cc가 수행 |
| `strict_profiles` | OPTIONAL | 검증 강도 토글(기본 off), 모델 수치를 안 움직여 liveness sweep이 볼 수 없음 |
| `profile_precision` | OPTIONAL | 프로파일 출처 주석; 미선언 = UNCALIBRATED 표기(정직한 기본값) |
| `softmax_operand_residency` | **UNCOVERED backlog** | dram·glb 분기가 비용에 영향 — glb 분기를 어떤 config도 실행하지 않음. 닫으려면 glb 선언 config + traffic delta 게이트 필요 (OPEN으로 기록) |

## 조치 2 — conv fixture + frontend 게이트

* `frontend/fixtures/conv_relu_{graph,exec}.json` + `conv_relu.map` — conv 3→8
  3×3 + fused ReLU (SFU conv 경로까지 커버).
* `workload_graph_test.cc`를 model_name 분기로 일반화 — 두 fixture 모두
  loader/adapter 계약을 고정. `run_pytorch_frontend_validation.sh`가 둘 다 실행.
* **`validation/frontend/check.py` 신설 (FE1–FE4)**: FE1 compile 결정성(fixture
  비트 재현), FE2 nebula parity(위 분석을 게이트로 영구화 — traffic 동일 +
  crit 차 == SFU busy 항등), FE3 conv MAC/SFU/commit 기하 항등, FE4 fail-loud
  부정 경로.

## 회귀

**26 스위트 26/26 전체 통과** (frontend 신설 포함, knobs 복귀) + unit test 2종 +
`check_timing --check-baseline --check-traffic` exit 0. 이 세션 최초의 전량 그린.

---

# 후속 — 실제 torch.export 캡처 검증 (미검증 항목 종결)

PyTorch **2.13.0+cpu**를 venv에 설치하고, fixture와 같은 기하의 실제 모델
(`nn.Linear(64,64)+ReLU`, `nn.Conv2d(3,8,3,bias=False)+ReLU`)을
`torch.export` 경유로 캡처했다.

| 검증 | 결과 |
| --- | --- |
| export → validate → compile | 양 모델 통과 |
| export 결정성 | 같은 factory 재실행 시 graph SHA 동일 (adb62bff…) |
| lowering 의미 동일성 | 실제 캡처의 operations(kind/activation/geometry)가 **fixture와 완전 동일** (linear·conv 모두) |
| run-ir 타이밍 | 실제 캡처 실행 결과가 fixture 실행과 **frontend 출처 주석 제외 비트 일치** |
| 미지원 op fail-loud | MaxPool1d 포함 모델 → `unsupported PyTorch operation aten.max_pool1d.default at node max_pool1d`로 compile 거부 |

발견한 사용상 유의점 1건: 실제 export의 conv op id는 torch.export 노드명을 따라
**"conv2d"**가 된다(합성 fixture는 "conv"). mapping의 `op_id`는 exported 노드명과
일치해야 하며, 불일치 시 binding 검증이 정확한 에러로 거부함을 확인했다.

재현: `pip install --index-url https://download.pytorch.org/whl/cpu torch numpy`
후 `PYTHONPATH=. python -m frontend.pytorch.cli export --factory <module>:<fn> ...`.

---

# 2026-08-31 후속 — pool fixture + SFU 계약을 FE 게이트에 추가

pool2d 작업(working tree, 08-28 18시 이후 안정)이 검증 가능한 상태임을 확인하고
요청대로 게이트에 편입했다. 그 사이 병행 세션이 FE 게이트를 residual-DAG(FE5,
BN/add/pool/concat 항등)로 확장해 둔 상태였으므로, 겹치지 않는 격차 — **전용
pool fixture(실제 캡처)와 pool의 SFU 요구 계약** — 만 추가했다.

## 추가된 것

* **`frontend/fixtures/pool_chain_{graph,exec}.json` + `.map`** — **실제
  torch.export(2.13.0+cpu)로 캡처**한 Conv2d(3→8,3×3) → MaxPool2d(2) →
  AvgPool2d(2). 합성이 아닌 실캡처를 fixture 소스로 사용(pool의 torch 경로까지
  함께 고정). op_id는 exported 노드명(conv2d)을 따름.
* **FE6** — pool chain 항등식: conv MAC 13,824(=8·8·8·3³), max pool scalar
  384(=128 출력×3 VMAX), avg pool 128(=96 VADD+32 VMUL), 층별 kind/mode 표기.
* **FE7** — pool SFU 계약: ① [sfu] 없는 config →
  `requires an [sfu] section for non-MAC timing` 거부, ② `supported_ops = relu`
  (vmax 제외) 변형 config → `needs SFU primitive 'vmax', outside this
  architecture's supported_ops contract` 거부 — 실행 유닛 부재는 zero-cost
  기본값이 아니라 아키텍처 사실이라는 계약의 검증.
* FE1 결정성 목록에 pool fixture 편입(4종), loader unit test에
  pool-chain 분기(3 op, pool은 MAC mapping 없음 확인), sh 스크립트에 등재.

## 회귀

26 스위트 **26/26 전체 통과**, unit test(4 fixture 전부 PASS), `check_timing
--check-baseline --check-traffic` exit 0. 유의: 이 검증은 pool 작업이 포함된
**미커밋 working tree** 위에서의 결과다 — 해당 작업이 커밋되면 FE6/FE7이 그
상태를 계속 고정한다.

---

# 2026-08-31: flatten alias 소거 — 분류망 end-to-end 개통

## 배경

LeNet급 실모델(conv→pool→…→flatten→FC→softmax) 실캡처가
`unsupported PyTorch operation aten.flatten.using_ints`에서 컴파일 거부됨을
실측으로 확인(08-30). flatten은 contiguous 텐서의 원소 보존 재해석이라
데이터 이동·산술·트래픽이 물리적으로 0 — "비용 0짜리 연산 추가"가 아니라
**연산 자체를 소거하고 텐서를 storage의 별명(alias)으로 처리**하는 것이 정직한
모델이다. 연산으로 넣으면 lifetime 이중 계상과 가짜 DRAM/GLB 트래픽이 생긴다.

## 설계 — alias 주석 방식 (재배선 없음)

* **lowering** (`_RESHAPE_OPS` = flatten.using_ints / view.default /
  reshape.default): 노드를 소거하고 exec IR에 operation을 만들지 않는다. 출력
  텐서는 유지하되 `alias_of: <storage>` 주석 추가(사슬은 lowering에서 단일
  storage로 해소). 검증: dtype·원소수 보존, layout contiguous,
  storage_offset 0, strides가 row-major(크기 1 차원 제외) — 아니면
  fail-loud ("eliding it would hide a real copy"). coverage는
  captured==lowered 항등 유지(소거 노드도 lowered로 계상, 스키마 무변경 →
  기존 fixture 비트 불변).
* **소비자는 alias id를 그대로 참조** — C++ geometry 검증과 nebula 브리지가
  view의 shape(예: Linear input_features=256)를 자동으로 올바르게 사용.
* **lifetime canonical 통일** (`workload_lifetime_t::canonical`): 소비자 수·
  GLB 상주·바이트 부기를 전부 storage id 기준으로 — alias를 통해 읽는 소비자가
  producer 버퍼를 살려둔다. graph output 판정도 storage 해소 집합으로.
* **npu.cc 분류 2곳** (`runtime_bytes`, `graph_operand_stream`)은
  `storage_tensor()` 기준(파라미터/graph-input 분류가 view가 아닌 storage를
  따름).
* **검증 3중화**: Python `_validate_aliases`(미선언 대상/alias-of-alias/
  dtype·bytes 불일치/graph input alias/op-produces-alias 거부) + C++ 로더
  동형 검증 + 토폴로지 produced 전파(storage 해소).

## 실증 — LeNet end-to-end

실제 torch.export 캡처 `lenet-fixture`(Conv 1→6 5×5 +ReLU → MaxPool2 →
Conv 6→16 5×5 +ReLU → MaxPool2 → [flatten 소거] → Linear 256→120 +ReLU →
Linear 120→10 → Softmax), gemmini_sfu.cfg로 7개 연산 전 층 타이밍 실행:

| 층 | 항등식 | 결과 |
| --- | --- | --- |
| conv1 MACs | 6·24·24·25 = 86,400 | 정확 |
| pool1 scalar | 6·12·12·3 = 2,592 | 정확 |
| conv2 MACs | 16·8·8·150 = 153,600 | 정확 |
| pool2 scalar | 16·4·4·3 = 768 | 정확 |
| fc1 MACs | 120·256 = 30,720 | 정확 |
| fc2 MACs | 10·120 = 1,200 | 정확 |
| softmax scalar | 9+10+10+9+1+10 = 49 | 정확 |

**핵심 증거**: fc1의 residency = `input GLB-resident`, GLB 점유 256→120 —
Linear가 flatten alias를 **통과해** pool2 storage를 GLB에서 직접 소비(가짜
DRAM 왕복 없음). 전체 GLB 점유 사슬 0→3456→864→1024→256→120→10→0이 각
단계에서 정확한 텐서 바이트와 일치.

## 편입

* fixture: `lenet_{graph,exec}.json` + `lenet.map`(4 MAC 섹션, op_id =
  exported 노드명), 재컴파일 결정성 확인.
* Python unit 6종 추가(소거/사슬 해소/비연속 거부/원소수 거부/alias-of-alias/
  미선언 대상) — 16/16 통과.
* loader unit lenet 분기: alias 필드 + **lifetime이 alias 경계를 넘어 상주를
  유지·해제**함을 직접 단언 — 5 fixture 전부 PASS.
* FE 게이트 FE8 신설 + FE1 5종 확장: 7층 정확(8번째 층 부재 = 소거 증명),
  MAC/scalar 항등 7건, `input GLB-resident` + `256 -> 120` 경계 검사.

## 남은 것 (정직한 한계)

비연속/offset view의 materialization(진짜 복사 비용 모델)은 의도적으로 거부
유지. `adaptive_avg_pool2d`(ResNet류)와 MaxPool1d는 여전히 미지원 —
fail-loud로 남음.

## 부록: lenet-fixture 재현 factory (torch 2.13.0+cpu)

```python
class _LeNetLike(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.c1 = torch.nn.Conv2d(1, 6, 5)
        self.c2 = torch.nn.Conv2d(6, 16, 5)
        self.pool = torch.nn.MaxPool2d(2)
        self.f1 = torch.nn.Linear(16 * 4 * 4, 120)
        self.f2 = torch.nn.Linear(120, 10)
        self.relu = torch.nn.ReLU()

    def forward(self, x):
        x = self.pool(self.relu(self.c1(x)))
        x = self.pool(self.relu(self.c2(x)))
        x = torch.flatten(x, 1)
        x = self.relu(self.f1(x))
        return torch.softmax(self.f2(x), dim=-1)

def lenet():
    return _LeNetLike().eval(), (torch.randn(1, 1, 28, 28),), {}
```

`python3 -m frontend.pytorch.cli export --factory <module>:lenet --output
frontend/fixtures/lenet_graph.json --model-name lenet-fixture` 후 compile.
