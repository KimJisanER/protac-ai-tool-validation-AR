# **안드로겐 수용체 표적 PROTAC의 탈모 치료제 가능성 평가 및 AI 기반 신규 분자 설계**

## 

## **1\. 서론: 타겟 선정 이유 및 AR PROTAC의 임상적 가능성**

안드로겐성 탈모(Androgenetic Alopecia, AGA)는 5-알파 환원효소에 의해 생성된 디하이드로테스토스테론(DHT)이 모낭의 안드로겐 수용체(Androgen Receptor, AR)와 결합하여 모발의 소형화와 탈락을 유발하는 만성 질환이다.1 기존의 탈모 치료제인 경구용 피나스테리드 등은 수용체 활성을 억제하기 위해 약물이 타겟에 지속적으로 결합해야 하는 '점유 구동(Occupancy-driven)' 방식을 따르며, 이는 필연적으로 전신적인 안드로겐 억제와 성기능 장애 등의 심각한 부작용을 동반한다.3  
이를 극복하기 위한 혁신적인 대안으로, 표적 단백질 자체를 유비퀴틴-프로테아좀 시스템(UPS)을 통해 원천적으로 분해해 버리는 AR 표적 PROTAC(Proteolysis-Targeting Chimeras) 기술이 주목받고 있다. 점유가 아닌 '사건 구동(Event-driven)' 방식으로 촉매처럼 작용하는 이 기술은 국소 투여용 치료제 개발에 최적화되어 있다.5 실제로 Kintor Pharmaceutical이 개발한 국소 도포용 AR PROTAC인 GT20029는 최근 중국에서 진행된 임상 2상 시험에서 그 우수성이 증명되었다. GT20029를 투여한 환자군(특히 0.5% 매일 도포 및 1.0% 주 2회 도포 군)은 위약군 대비 12주 차에 통계적으로 유의미한 모발 수(TAHC) 및 모발 두께(TAHW) 증가를 달성하였다.6 무엇보다 전신 혈중 흡수량이 정량 한계 미만에 머물러 심각한 부작용 없이 두피에서만 안전하게 작용함이 입증되며 차세대 탈모 치료제로서의 실질적인 가능성을 열었다.8

## 

## **2\. AR PROTAC 개발의 한계: 분해 효율과 피부 투과성의 딜레마**

임상적 타당성이 입증되었음에도 불구하고, 새로운 국소 도포용 AR PROTAC을 설계하는 것은 화학적으로 매우 난해한 미션이다. 탈모약으로서 기능하기 위해서는 표적 단백질의 분해(Degradation) 효율이 충분히 높아야 할 뿐만 아니라, 견고한 두피의 각질층을 원활히 통과하되 전신 혈류로는 빠져나가지 않고 모낭 깊숙이 머무르는 '피부 체류성(Skin retention)'을 동시에 만족해야 하기 때문이다.10  
문제는 대부분의 PROTAC 분자가 분자량 700\~1200 Da, 높은 위상학적 극성 표면적(TPSA)을 지녀 전통적인 신약 개발 기준인 '리핀스키의 5법칙'을 훌쩍 벗어나는 bRo5(beyond Rule of Five) 화학 공간에 존재한다는 점이다.12 이로 인해 단순한 수동 투과도를 측정하는 PAMPA 분석 등에서 PROTAC은 대부분 정량 한계 미만의 저조한 투과성을 보인다.14  
분해 효율을 유지하면서 피부 투과성을 최적화하는 과정은 극도로 까다롭다. 예를 들어, 침습적인 마이크로니들 투여가 필요했던 초기 AR PROTAC인 TJA-107을 최적화한 화합물 C6의 개발 사례를 보면, 링커의 강성을 조절하고 구조 내 TPSA를 줄이며 지질친화도(clogP)를 높이는 정교한 조율을 통해 나노몰(nM) 수준의 강력한 AR 분해능과 1700배가 넘는 우수한 피부 대 혈장 농도 비(Skin-to-plasma ratio)를 확보해 냈다.10 또한, PROTAC과 같이 유연한 분자는 주변 용매 환경에 따라 분자 내 수소 결합(IMHB)을 형성하여 스스로 웅크리며 극성을 숨기는 '분자 카멜레온성(Molecular chameleonicity)'을 띠는데, 이러한 노출된 극성 표면적(EPSA)의 동적인 궤적 등 수많은 변수들이 복합적으로 작용하여 투과성이 결정된다.15 즉, 분자 구조의 미세한 조합 변화가 형태와 투과성에 미치는 영향을 인간의 직관과 힘만으로 예상하는 것은 사실상 불가능에 가깝다.

## **3\. 화학정보학 및 AI 기반 정량적 예측의 필요성**

인간의 힘으로 설계하기 까다로운 bRo5 화합물의 특성과, 낮은 용해도 및 막 내 갇힘 현상으로 인해 기존 체외(in vitro) 세포 투과도 분석법(PAMPA, Caco-2)이 빈번하게 한계를 보이는 점을 고려할 때14, 축적된 실험 데이터를 바탕으로 분자의 물성을 파악하고 인공지능(AI)을 통해 정량적으로 예측해보는 화학정보학적 접근은 점차 그 의미와 중요성이 커지고 있다.  
1,600개 이상의 PROTAC 구조와 DC50 등의 분해 활성 데이터를 보유한 PROTAC-DB나 방대한 생물활성 데이터를 제공하는 ChEMBL 등의 공개 데이터베이스를 활용하면18, 방대한 화학 공간 내에서 직관을 뛰어넘는 구조-활성 패턴을 도출할 수 있다. 구조가 공개되지 않은 GT20029에 의존하는 대신, 이러한 기존 데이터를 학습한 AI 예측 모델(단백질-리간드 복합체 구조 예측, 투과도 기반 ML 모델 등)을 활용하면 약물의 효능과 피부 체류성이라는 상충하는 다중 파라미터를 인실리코(In silico) 상에서 효과적으로 동시 최적화할 수 있다.

## **4\. 본 보고서의 연구 계획 및 목표**

이상의 배경을 바탕으로, 본 보고서에서는 다음과 같은 연구 흐름을 통해 AR PROTAC 탈모 치료제의 새로운 가능성을 탐구하고자 한다.

1. **화학적 분포 공간 탐색:** PROTAC-DB 및 ChEMBL 데이터베이스를 활용하여, 활성이 이미 알려진 AR 표적 PROTAC 분자들의 구조적 특징(링커 및 E3 리가아제 리간드 조합)과 분자량, TPSA 등의 물리화학적 분포 공간을 살펴본다.  
2. **AI 도구의 예측 정확도 평가:** 확보된 AR PROTAC 데이터를 바탕으로 결합 포즈(Binding pose)와 결합 친화도(Affinity)를 예측하는 최신 AI 구조 예측 및 도킹 도구들을 적용해 보고, 그 예측 성능과 정확도를 기존 실험 데이터와 비교하여 비판적으로 평가한다.  
3. **투과성(Permeability) 예측 및 검증:** 분자의 이론적 극성 표면적(TPSA) 및 카멜레온성을 대변하는 분자 기술자를 바탕으로 피부 투과성 및 체류성에 대한 기계 학습 예측을 수행하여, 본 인실리코 예측 방법론의 타당성(Validity)을 확인한다.  
4. **신규 분자 구조 제시:** 구축된 예측 및 평가 파이프라인을 기반으로, AR 분해(Degradation) 효능이 충분히 높을 것으로 예측되면서도 동시에 우수한 피부 투과성을 지녀 국소 탈모약으로서의 물성이 극대화된 새로운 PROTAC 후보 화합물을 최종적으로 도출하여 제시할 것이다.

## **5\. 데이터/구조 확보**

PROTAC DB 3.0

Jingxuan Ge, Shimeng Li, Gaoqi Weng, Huating Wang, Meijing Fang, Huiyong Sun, Yafeng Deng, Chang- Yu Hsieh, Dan Li, Tingjun Hou. PROTAC-DB 3.0: an updated database of PROTACs with extended pharmacokinetic parameters. *Nucleic Acids Research*, 2025 Jan 6;53(D1):D1510-D1515. doi: [10.1093/nar/gkae768](https://doi.org/10.1093/nar/gkae768). 

데이터셋 파일 : "jm4c02226\_si\_002.csv"  
DOI 에 담겨있는 column을 기반으로 논문 공개 날짜도 수집해줘  
ProtacDB\_MTL\_CLS.csv 이것은 이전버전 데이터 셋인데, 이 데이터 셋에 없는 jm4c02226\_si\_002.csv 내의 자료를 추려줘. 그리고 이 차이를 기반으로 time split을 해서 AI 모델들의 성능을 확인하고 싶어.

## **6\. 방법과 선택 이유**

Chemical Space 탐색

Degradation AI model  
[https://github.com/ABILiLab/DegradeMaster](https://github.com/ABILiLab/DegradeMaster)    
[https://github.com/fenglei104/DeepPROTACs](https://github.com/fenglei104/DeepPROTACs)  
[https://github.com/LiZhang30/AiPROTAC](https://github.com/LiZhang30/AiPROTAC)  
[https://github.com/PROTACs/PROTAC-STAN](https://github.com/PROTACs/PROTAC-STAN)  
[https://github.com/bryanc5864/DegradoMap](https://github.com/bryanc5864/DegradoMap)

Cell membrane permeability  
[https://github.com/ycu-iil/PROTAC-TS](https://github.com/ycu-iil/PROTAC-TS)  
Potts and Guy 다중 선형 회귀 모델 (Skin)  
Mitragotri 결정론적 예측 모델 (Deterministic Model)

co-fold 모델  
[https://github.com/jwohlwend/boltz](https://github.com/jwohlwend/boltz)

case-study  
[https://pubs.acs.org/doi/10.1021/acs.jmedchem.4c02226](https://pubs.acs.org/doi/10.1021/acs.jmedchem.4c02226)  
각 분자 별 결과는 별도로 파일: "jm4c02226\_si\_002.csv"

#### **참고 자료**

1. Androgenetic Alopecia \- StatPearls \- NCBI Bookshelf \- NIH, 6월 16, 2026에 액세스, [https://www.ncbi.nlm.nih.gov/books/NBK430924/](https://www.ncbi.nlm.nih.gov/books/NBK430924/)  
2. Treatment options for androgenetic alopecia: Efficacy, side effects, compliance, financial considerations, and ethics \- PMC, 6월 16, 2026에 액세스, [https://pmc.ncbi.nlm.nih.gov/articles/PMC9298335/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9298335/)  
3. The diagnosis and treatment of androgenetic alopecia: a review of the most current management | Śliwa | Forum Dermatologicum, 6월 16, 2026에 액세스, [https://journals.viamedica.pl/forum\_dermatologicum/article/view/95391/76623](https://journals.viamedica.pl/forum_dermatologicum/article/view/95391/76623)  
4. PROTACs: Past, Present and Future \- PMC \- NIH, 6월 16, 2026에 액세스, [https://pmc.ncbi.nlm.nih.gov/articles/PMC10237031/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10237031/)  
5. Technical characteristics of PROTAC technology in drug research and development, 6월 16, 2026에 액세스, [https://axispharm.com/technical-characteristics-of-protac-technology-in-drug-research-and-development/](https://axispharm.com/technical-characteristics-of-protac-technology-in-drug-research-and-development/)  
6. Exciting Breakthrough in Hair Loss Treatment: GT20029 Emerges as a Promising Topical Option for Androgenetic Alopecia \- Bauman Medical, 6월 16, 2026에 액세스, [https://www.baumanmedical.com/gt20029-hair-loss-treatment-phase-2-results/](https://www.baumanmedical.com/gt20029-hair-loss-treatment-phase-2-results/)  
7. Full article: Efficacy and safety of topical GT20029 in male patients with androgenetic alopecia: a multicenter, randomized, double-blind, placebo-controlled phase 2 study \- Taylor & Francis, 6월 16, 2026에 액세스, [https://www.tandfonline.com/doi/full/10.1080/09546634.2025.2574304](https://www.tandfonline.com/doi/full/10.1080/09546634.2025.2574304)  
8. Kintor Pharma Announces Its Trial for Androgenetic Alopecia Reaches Primary Endpoint, 6월 16, 2026에 액세스, [https://clival.com/news/kintor-pharma-announces-its-trial-for-androgenetic-alopecia-reaches-primary-endpoint](https://clival.com/news/kintor-pharma-announces-its-trial-for-androgenetic-alopecia-reaches-primary-endpoint)  
9. First-in-Class Topical GT20029 Demonstrates Promising Phase 2 Efficacy and Tolerability for AGA | Dermatology Times, 6월 16, 2026에 액세스, [https://www.dermatologytimes.com/view/first-in-class-topical-gt20029-demonstrates-promising-phase-2-efficacy-and-tolerability-for-aga](https://www.dermatologytimes.com/view/first-in-class-topical-gt20029-demonstrates-promising-phase-2-efficacy-and-tolerability-for-aga)  
10. Expanding the Scope of PROTACs: Opportunities and Challenges in ..., 6월 16, 2026에 액세스, [https://pmc.ncbi.nlm.nih.gov/articles/PMC12670424/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12670424/)  
11. Discovery of a Novel Non-invasive AR PROTAC Degrader for the Topical Treatment of Androgenetic Alopecia \- PubMed, 6월 16, 2026에 액세스, [https://pubmed.ncbi.nlm.nih.gov/39641607/](https://pubmed.ncbi.nlm.nih.gov/39641607/)  
12. Refinement of Computational Access to Molecular Physicochemical Properties: From Ro5 to bRo5 \- PMC, 6월 16, 2026에 액세스, [https://pmc.ncbi.nlm.nih.gov/articles/PMC9511483/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9511483/)  
13. Targeted Protein Degradation with PROTACs and Molecular Glues \- Blog, 6월 16, 2026에 액세스, [https://blog.crownbio.com/targeted-protein-degradation-with-protacs-and-molecular-glues](https://blog.crownbio.com/targeted-protein-degradation-with-protacs-and-molecular-glues)  
14. Systematic Investigation of the Permeability of Androgen Receptor ..., 6월 16, 2026에 액세스, [https://pmc.ncbi.nlm.nih.gov/articles/PMC7429968/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7429968/)  
15. A rational control of molecular properties to discover new oral drugs in the beyond-Rule-of-5 (bRo5) chemical space, 6월 16, 2026에 액세스, [https://iris.unito.it/retrieve/1ba8b290-4ce5-448f-8fff-833f9a0d754c/Thesis\_Diego\_Garcia\_Jimenez\_final\_version.pdf](https://iris.unito.it/retrieve/1ba8b290-4ce5-448f-8fff-833f9a0d754c/Thesis_Diego_Garcia_Jimenez_final_version.pdf)  
16. Enhancing Permeability Through Exposed Polar Surface Area (EPSA) for Beyond Rule of Five (bRo5) Drug Candidates \- WuXi AppTec DMPK, 6월 16, 2026에 액세스, [https://dmpkservice.wuxiapptec.com/articles/373-enhancing-permeability-through-exposed-polar-surface-area-epsa-for-beyond-rule-of-five-bro5-drug-candidates/](https://dmpkservice.wuxiapptec.com/articles/373-enhancing-permeability-through-exposed-polar-surface-area-epsa-for-beyond-rule-of-five-bro5-drug-candidates/)  
17. Conquering the beyond Rule of Five Space with an Optimized High-Throughput Caco-2 Assay to Close Gaps in Absorption Prediction \- PMC, 6월 16, 2026에 액세스, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11280027/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11280027/)  
18. PROTAC-DB, 6월 16, 2026에 액세스, [https://cadd.zju.edu.cn/protacdb/about](https://cadd.zju.edu.cn/protacdb/about)  
19. PROTAC-DB \- Database Commons, 6월 16, 2026에 액세스, [https://ngdc.cncb.ac.cn/databasecommons/database/id/7313](https://ngdc.cncb.ac.cn/databasecommons/database/id/7313)