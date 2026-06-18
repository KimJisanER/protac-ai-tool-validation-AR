# 현재 PROTAC In-silico 예측 도구는 "의미 있는 예측"이 가능한가?
### — 안드로겐 수용체(AR)를 사례로 한 분해 활성·투과성·삼원복합체(co-fold) 예측 모델의 비판적 평가

화학정보학·분자설계 기말 프로젝트. PROTAC 신약 설계용 공개 AI 도구(분해 활성·막 투과성·삼원복합체 co-fold 예측)가 **학습에 없던 신규 화합물에 대해 실제로 신뢰할 만한 예측을 하는지**를, 데이터가 풍부하고 임상 검증된 **AR(UniProt P10275)** 을 시금석으로 비판적으로 검증한다.

## 핵심 결론
도구들은 **자신이 학습한 분포(벤치마크)에서는 우수**해 보이나, **학습 미관측 화합물의 전향적(prospective) 예측에서는 성능이 크게 하락**하며 순위조차 평가 설계에 따라 뒤집힐 만큼 불안정하다.

| 도구 부류 | 보고/in-distribution | 본 검증(전향적·AR) | 판정 |
|---|---|---|---|
| 분해 활성 (PROTAC-STAN/Ribes/DegradeMaster) | AUROC 0.85–0.91 | pooled 0.57·타겟 내 ≈0.45(무작위 부근)·도구 간 κ=0.13 | 미성숙 |
| 투과성 (Potts-Guy/PROTAC-TS) | (소분자/펩타이드 검증) | AD 외삽·예측 압축, 두 투과지표 무상관 | 제한적 |
| Co-fold (Boltz-2) | (AF3급) | ligand_ipTM 0.937 / protein_ipTM 0.671(단일시드) | 부분적 |

→ **단독 예측엔 아직 미성숙. 실험 검증 전제의 약한 사전확률(weak prior)·다중 도구 합의로만 사용 권장.**

## 저장소 구성
- `AR PROTAC 탈모 치료제 연구 조사 (본문 완성).md` — **제출용 본 보고서**(①~⑥)
- `REPORT.md` — 동일 분석의 압축 보고서
- `PLAN.md` / `ENV_SETUP.md` / `SETUP.md` — 실행 계획·환경 메모
- `카멜레온성 예측 방법론 정리.md` — 분자 카멜레온성 예측 방법론 서베이
- `scripts/` — 전 분석 코드(데이터 동결, STAN/Ribes 추론·평가, Potts-Guy/PROTAC-TS 투과, 게재연도, 화학공간 HTML 등)
- `outputs/` — 결과 표(CSV)·그림(PNG)·로그
- `*.html` — 인터랙티브 화학공간(`chemspace_view1_AR`, `chemspace_view2_split`)

## 데이터(미포함) — 직접 내려받기
저작권·라이선스 때문에 원본 데이터는 repo에 포함하지 않았다(`.gitignore`).
- **PROTAC-DB 3.0**(`protac.xlsx`) 및 이전 버전(`ProtacDB_MTL_CLS.csv`): http://cadd.zju.edu.cn/protacdb/ (사용 정책 준수)
- **case study**(`jm4c02226_si_002.csv`, 논문 PDF): *J. Med. Chem.* 2024, 67(24), 22218–22244, DOI 10.1021/acs.jmedchem.4c02226 의 Supporting Information
- 평가한 도구는 각 공식 repo에서 설치: PROTAC-STAN, PROTAC-Degradation-Predictor(Ribes), DegradeMaster, PROTAC-TS, Boltz-2(co-fold)

## 환경
- 메인 분석 `protac` (Python 3.10, torch 2.8+cu128, RDKit 2026.03.1) · 투과 `protac_ts`(TabPFN) · co-fold `protac_boltz`(Boltz-2) · 도구 교차검증 `protac_ribes`

## 재현 개요
원본 데이터를 위 출처에서 받아 배치한 뒤 `scripts/`의 단계별 스크립트를 순서대로 실행한다(데이터 동결 → STAN 전향 평가 → 다중 도구 비교 → 투과 → co-fold → 게재연도/화학공간). 세부 경로·산출물 매핑은 본 보고서 부록 참조.

---
> 본 저장소의 표/그림 수치는 `scripts/` 출력(`outputs/`)을 단일 출처로 한다. AI 도구의 보고 수치를 그대로 신뢰하지 말고 누수·타겟·라벨을 통제해 재검증하라는 것이 본 프로젝트의 메타-교훈이다.
