
# EGR Knockin-Knockout RNA-seq Analysis
이 프로젝트는 EGR2/EGR3 넉인 및 넉아웃 마우스에서 CD4+/CD8+ T세포의 전사체 데이터를 분석하는 코드입니다.

## 📂 폴더 구조
- `1_Preprocessing_Norm_TCC.R` : 데이터 전처리 및 정규화
- `2_PCA_Dendrogram_EGR_Tcells.R` : PCA 및 덴드로그램 분석
- `3_DEG_Analysis_EGR2_EGR3.R` : 차등 발현 유전자(DEG) 분석
- `4_GO_KEGG_Enrichment_EGR2_EGR3.R` : 기능적 분석 (GO/KEGG)
- `5_Visualization_Plots_EGR.R` : 최종 시각화 및 보고서 생성
- `data/` : 원본 및 전처리된 데이터 저장
- `results/` : 분석 결과 저장
- `figures/` : 생성된 플롯 저장
- `scripts/` : 보조 함수 및 서브 스크립트 저장