# transcriptome과 RNAseq 용어 범위 차이
- Transcriptome analysis(전사체 분석)과 RNA-seq(RNA sequencing)은 서로 관련이 깊지만 동의어는 아닙니다.
- RNA-seq은 전사체 분석을 수행하는 방법 중 하나이며, 전사체 분석은 RNA-seq뿐만 아니라 다른 기술을 포함하는 더 넓은 개념입니다.

## 🔹 Transcriptome Analysis (전사체 분석)란?
📌 **전사체(transcriptome)**는 세포나 조직 내에서 발현되는 모든 RNA(전사된 유전자)를 총칭하는 개념입니다. <br>
📌 **전사체 분석(Transcriptome Analysis)**은 이 RNA들을 분석하여 유전자 발현 수준, 조절 기작, 차등 발현(DEG), 기능적 연관성(GO, KEGG) 등을 연구하는 과정을 의미합니다. <br>

### 전사체 분석의 주요 목표
- 특정 조건(예: 질병 vs 건강한 상태, 넉인 vs 넉아웃)에서 어떤 유전자들이 발현되는지 분석
- 차등 발현 유전자(DEG, Differentially Expressed Genes) 분석
- 기능적 분석 (Gene Ontology, KEGG pathway 등)
- 전사체 수준에서 유전자 네트워크 분석

### 전사체 분석을 수행하는 다양한 방법
전사체 분석에는 여러 기술이 있으며, 그중 RNA-seq은 가장 널리 사용되는 방법입니다.

| 분석 방법                     | 설명                                      |
|------------------------------|-------------------------------------------|
| **RNA-seq (RNA sequencing)** | 차세대 염기서열 분석(NGS) 기반 RNA 발현 분석 |
| **Microarray**               | 마이크로어레이 칩을 이용한 유전자 발현 프로파일링 |
| **qRT-PCR**                  | 특정 유전자의 발현량을 정량적으로 분석        |
| **NanoString**               | 특정 유전자 패널의 발현을 직접 측정            |
| **Ribo-seq**                 | 리보솜-연관 RNA를 분석하여 번역 중인 유전자 연구 |


## 🔹 RNA-seq (RNA Sequencing)이란?
📌 **RNA-seq(RNA sequencing, RNA 시퀀싱)**은 NGS(차세대 염기서열 분석, Next-Generation Sequencing)를 이용하여 RNA 서열을 읽고 유전자 발현을 정량화하는 기술입니다. <br>
📌 전사체 분석을 수행하는 대표적인 방법으로, 거의 모든 전사체 연구에서 사용됩니다. <br>

### RNA-seq의 주요 특징
✅ 고해상도: 기존 마이크로어레이보다 정밀한 유전자 발현 분석 가능 <br>
✅ 데이터 스펙트럼이 넓음: 잘 알려진 mRNA뿐만 아니라, lncRNA, miRNA, circRNA 등 다양한 RNA를 분석 가능 <br>
✅ 새로운 유전자/전사체 탐색 가능: 기존에 알려지지 않은 유전자, 스플라이싱 변이 탐색 가능 <br>
✅ NGS 기반으로 빠르고 대량 데이터 처리 가능 <br>

### RNA-seq의 주요 분석 과정
1. RNA 추출 및 라이브러리 제작
2. NGS 시퀀싱 수행 (Illumina, PacBio 등)
3. Raw 데이터 처리 (FASTQ 파일)
4. Read mapping (유전체 정렬)
5. 유전자 발현량 정량화 (Count, TPM, FPKM)
6. 차등 발현 분석(DEG)
7. 기능적 분석 (GO, KEGG, GSEA 등)

## 🔹 결론: Transcriptome Analysis vs RNA-seq
✅ RNA-seq은 전사체 분석을 수행하는 대표적인 방법이다. <br>
✅ Transcriptome analysis(전사체 분석)는 RNA-seq뿐만 아니라 마이크로어레이, qRT-PCR 등 여러 방법을 포함하는 더 넓은 개념이다. <br>
✅ 전사체 분석은 생물학적 질문을 해결하기 위한 연구 방법이며, RNA-seq은 이를 수행하는 강력한 도구이다. <br>
✅ 즉, RNA-seq을 사용해서 전사체 분석을 수행할 수 있다! 🚀 <br>