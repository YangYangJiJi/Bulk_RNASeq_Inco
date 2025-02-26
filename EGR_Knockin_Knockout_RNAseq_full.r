### -------------------------###
### R 경로설정               ###
### -------------------------###

#경로설정 방법1
getwd() #현재 working directory 위치 확인
setwd("C:/myGithub/Transcriptome_Analysis/EGR_Knockin_Knockout_RNAseq") #working directory 직접 설정

#방법2
#또는 이렇게 해도 됨.
#현재 스크립트가 저장된 위치를 mainDir 변수에 넣음
mainDir=dirname(rstudioapi::getSourceEditorContext()$path) 
mainDir #경로 맞는지 확인
setwd(mainDir) #현재 스크립트가 저장된 위치로 working directory 설정
getwd() #경로 맞는지 확인


### -------------------------###
### 데이터 input             ###
### -------------------------###

#BiocManager 
library(BiocManager)
BiocManager::install("ArrayExpress")
library(ArrayExpress)

#ArrayExpress 데이터베이스에서 데이터 다운로드
#1) ArrayExpress::getAE 함수 이용하기
dat = ArrayExpress::getAE("E-MTAB-5338",type="processed") # 분석 데이터 다운로드
#Downloading file : symonds_counts_UCSC_refGene_mm10.txt

#다운로드 정보 확인
dat$dataFiles$type
dat$dataFiles$file

#2) Url download
dat$dataFiles$url
url <- "ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/338/E-MTAB-5338/Files/symonds_counts_UCSC_refGene_mm10.txt"
destfile <- file.path(getwd(), "symonds_counts_UCSC_refGene_mm10.txt")

#download.file(다운로드할 url, 어디에 저장할지 변수명, method="auto") #이런식으로도 다운로드 받을 수 있음
download.file(url, destfile, method="auto") 
#사실 보통은 그냥 사이트에서 다운받음. R에서 다운받는 일은 적다.

counFile = dat$processedFiles
counFile = 'symonds_counts_UCSC_refGene_mm10.txt'

#여기서 잘 해야함. 데이터를 잘 가져와야 나중에 꼬이지 않을것이다. 정말 중요한 부분.
#이 부분은 raw data file을 직접 열어서 어떤 명령을 내릴지 분석가의 판단이 중요함.
dataDF = read.table(counFile,
                    sep = '\t', #seperate임. 구분자를 넣어준다. 이거는 탭으로 구분한다는 뜻
                    header = T, #내 샘플의 헤더를 그대로 가져온다. (True)
                    row.names=1, #첫번째 열을 row name으로 불러옴.
                    check.names = F) #컬럼 이름을 체크할거냐? #만약 T하면, 슬래시를 점으로 바꿔줌.

metadata = read.table(paste0(getwd(), '/metadata.txt'),
                      sep = '\t',
                      header = T,
                      quote='"', #안에 쌍따옴표 넣으면, 각 데이터에 씌워진 쌍따옴표 다 없애줌.
                      stringsAsFactors = F, #R 내부에서 factor를 따로 지정해서 sorting 하고 싶으면 True. #멘토님은 애초에 metadata부터 잘 정리해서 넣는걸 추천 
                      check.names = T)

head(dataDF)
head(metadata)


### ------------------------------###
### 발현된 마우스 유전자 필터링   ###
### ------------------------------###

#dataDF안에 metadata의 SampleName 순으로 정렬 
dataDF=dataDF[,metadata$Sample.Name]

# 이제 카운트파일, 메타데이터 다 불러오고, 두개의 열을 잘 정렬함. 
dim(dataDF)
colnames(dataDF)
class(dataDF) #샘플이름을 잘 선정하는 것도 중요하다. 
#전체 유전자 중, 발현된 유전자로 필터링 할거임.
#12개 샘플 중, 유전자가 발현되었다고 정의할 때, 한 유전자가 5이상의 카운트 값을 가지는데, 
#9개 샘플 이상에서 유전자가 발현되면 발현된걸로 치겠다. 이렇게 필터링


#연산할 떄 R이 이해하기 쉽도록 데이터 프레임을 매트릭스화 시킴.
dataDF = data.matrix(dataDF)
dataDF = dataDF[rowSums(dataDF>=5)>=9,] #뜯어볼거임.


# a=dataDF>=5 #일단 5이상의 카운트 값을 확인 #5이상이면 True로 반환
# b=rowSums(dataDF>=5) #rowSums는 행을 다 합함. 리스트가 유전자 갯수만큼 나오게 됨. 각 행마다 답이 1개씩 생기기에 
# c=rowSums(dataDF>=5)>=9 #9개 이상 발현된 것 True로 반환
# dataDF=dataDF[c,]
#이것을 한줄로 쓰면 dataDF = dataDF[rowSums(dataDF>=5)>=9,] 가 된다.
#지금까지 한게 발현된 마우스 유전자 필터링한거임.

dim(dataDF) #2만4천에서 1만3천개로 유전자가 줄었다. 필터링 되었다.


### --------------------###
### 정규화 - TCC 사용   ###
### --------------------###

### count파일을 Normalization할거임

#데이터의 정규화 및 차등발현 (DE)을 볼 수 있는 패키지 설치.
BiocManager::install("limma")
BiocManager::install("TCC, force=true") #전사체 분석에서 가장 많이 사용하는 패키지. DEG분석, 데이터 normalization 등 
library(limma)
library(TCC)

#new라는 함수안에 TCC라는 클래스르 만들겠다. 
#dataDF라는 카운트 넣어주고, metadata의 Group컬럼을 불러와서 factor함수에 넣음. 
#그러면 레벨이 지정됨. 
#as.numeric 함수 안에 넣으면 숫자로 변함. 
tcc=new('TCC', dataDF, group = as.numeric(factor(metadata$Group)))
class(tcc) #클래스 형태로 된 것을 확인.
show(tcc)
# TCC는 차등발현유전자(DEG) 분석을 수행하기 위한 BIOconductor 패키지 클래스
# RNA-seq 데이터를 입력으로 받아, 차등발현유전자를 식별하고 시각화 및 통계분석을 수행하는데 사용


## TCC 객체 생성 후 정규화 및 분석 수행
#tmm의 정규화 방법 이용해서 우리의 카운트 파일을 정규화 시킴
tcc=calcNormFactors(tcc, norm.method='tmm', test.method='edger',FDR=0.05)
show(tcc) #샘플 안에 norm.factor값이 계산되어 넣어짐.
# STEP 1 : 데이터 정규화
#TMM (trimmed Mean of M-values) 정규화 방법 : edgeR 패키지에서 제공하는 방법
#라이브러리 크기와 발현 값의 편향을 보정하는 데 자주 사용. 전사체 분석에서 많이씀.

#정규화된 데이터 추출
normDF=getNormalizedData(tcc)
head(normDF)
#지금까지 카운트 값을 정규화 완료함.
#카운트 파일에서 각 라이브러리 사이즈에 맞게, 각 분포에 맞게 정규화된 파일
#그러면 샘플간 비교가 가능해진다.


### --------------------###
### data melting        ###
### --------------------###

## 이제 샘플 QC를 진행할거임. 그림그리기
BiocManager::install('reshape2')
library(reshape2)

#wide 형식의 데이터를 long 형식으로 변환
dataDF_melt = melt(dataDF) #열의 정보들을 행으로 구조변환 시킴. melt함수로 
head(dataDF_melt)

dim(dataDF_melt) #아까 1만3천개에서 1만3천 x 12 = 16만2천개.

colnames(dataDF_melt) = c('gene','sample','count') #컬럼이름 설정.
#행 이름이었던 유전자 이름을 값으로 가지고 옴. (엄연히 다르다. 인덱싱 정보를 값으로 가져옴)


dataDF_melt = merge(dataDF_melt, metadata[,c('Sample.Name','Group')], by.x = 'sample', by.y='Sample.Name')
head(dataDF_melt)
#마치 엑셀의 vlookup과 같은 기능을 하는 merge함수
#샘플이름을 기준으로 합쳐짐.

### 정규화 한 데이터를 멜트시키고, 메타데이터와 합치는 과정
normDF_melt=melt(normDF)
colnames(normDF_melt)=c('gene','sample','count')

normDF_melt=merge(normDF_melt, metadata[,c('Sample.Name','Group')], by.x = 'sample', by.y='Sample.Name')
head(dataDF_melt, 3) #노말리 안된것
head(normDF_melt, 3) #노말리 된것

library(ggplot2) #데이터 시각화를 위한 강력한 도구로, 데이터를 이해하고 전달하는 데 큰 도움을 주는 패키지


### --------------------###
### box plot            ###
### --------------------###
library(RColorBrewer)
display.brewer.pal(7, "BrBG")
brewer.pal(7,"Pastel1")
#"#FBB4AE" "#B3CDE3" "#CCEBC5" "#DECBE4" "#FED9A6" "#FFFFCC" "#E5D8BD"

colnames(dataDF_melt)
ggplot(dataDF_melt) +
  geom_boxplot(aes(x=sample, y=log2(count+1), fill = Group)) + # 그룹별로 색깔 나옴
  theme_bw() + #테마
  labs(title = 'counts before normalization') + #라벨. 플랏 제목
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_fill_manual(values=c("Egr2_Kin_CD4" = "#FBB4AE",
                             "Egr2_Kin_CD8" = "#B3CDE3",
                             "Egr2/3_DKO_CD4" = "#DECBE4",
                             "Egr2/3_DKO_CD8" = "#CCEBC5"))

dev.off() #그림 삭제
ggsave('Fig.Boxplot-1.png', width=5, height=4) #플랏 그림 저장

### 이제 정규화 한 box plot 그림
ggplot(normDF_melt) +
  geom_boxplot(aes(x=sample, y=log2(count+1), fill = Group)) + # 그룹별로 색깔 나옴
  theme_bw() + #테마
  labs(title = 'counts after normalization') + #라벨. 플랏 제목
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave('Fig.Boxplot-2.png', width=5, height=4)


### --------------------###
### PCA plot            ###
### --------------------###

### PCA : 데이터의 고차원 정보를 저차원으로 축소하여 주요 경향성을 확인
#데이터의 분산을 최대한 보존하는 방향으로 새로운 축 (주성분, PC1, PC2)을 정의
#PC3까지 써서 공간으로, 3차원으로 만들수 있음. R에서 돌려볼 수 있다.
#PCA는 노말리 한 카운트 값으로 PCA 그릴 것

log2DF = log2(normDF+1)
PCA = prcomp(t(log2DF)) #t는 행에 샘플이 갈 수 있게끔 전치시킨다.

View(PCA)
names(PCA)

# 분산은 각 주성분(PC)의 기여도
# 값이 클수록 해당 주성분이 데이터의 변동성을 더 많이 설명
variances <- PCA$sdev^2 #표준편차에 제곱해서 분산구하기.
#총 분산 계산
total_variance <- sum(variances)

#각 주성분의 기여도 계산 
proportion_variance <- variances / total_variance

#누적 기여도 계산
cumulative_variance <- cumsum(proportion_variance)

pc1 = round(proportion_variance[1]*100,2) #소수점 둘째자리까지만 보겠다.
pc2 = round(proportion_variance[2]*100,2) #소수점 둘째자리까지만 보겠다.


# 각 변수가 주성분(PC)을 형성하는 데 기여하는 선형 결합 계수
head(PCA$rotation)
head(sort(abs(PCA$rotation[, "PC1"]), decreasing = TRUE))
#유전자의 절댓값이 클수록 PC1에 기여하는 부분이 많은 유전자.
#기여도가 큰 순서대로 분류함. sort

#메타데이터를 같이 쓰기 위해 merge
head(PCA$x) #주성분 점수 
pcaDF = merge(PCA$x, metadata[c('Sample.Name','Group')], by.x=0, by.y = 'Sample.Name') 
head(pcaDF)

unique(metadata$Group)
ggplot(pcaDF) +
  geom_point(aes(x=PC1, y=PC2, color=Group), size=2)+
  geom_hline(yintercept =0, color='red', linetype = 'dashed')+ #y축 그리기
  geom_vline(xintercept =0, color='red', linetype = 'dashed')+ #x축 그리기
  labs(title='PCA plot',
       x=paste0('PC1 (', round(pc1, 2), '%)', sep=''),
       y=paste0('PC2 (', round(pc2, 2), '%)', sep=''))+
  theme_bw() + #겉에 검정색 네모박스.+
  geom_text_repel(aes(x=PC1, y=PC2, label = Row.names), size=2) +
  scale_color_manual(values=c("Egr2_Kin_CD4" = "#FBB4AE",
                              "Egr2_Kin_CD8" = "#B3CDE3",
                              "Egr2/3_DKO_CD4" = "#DECBE4",
                              "Egr2/3_DKO_CD8" = "#CCEBC5"))

#PCA 정보로 
#PC1이 PC2보다 분산력이 더 크다. 
#파란계열은 CD8 / 빨간계열은 CD4 -> PC1 기준으로 보면 셀 타입별로 구분
#PC2에 따라서는 넉인 넉아웃이 구분됨. 
#셀 타입에 따른 유전 거리가 더 멀다.
#넉인 넉아웃에 따라 생명력이 크게 차이났다면 PC1에 의해 구분되었을 것(?)
#괄호 안에는 설명력을 나타냄. PC1은 61.74%의 설명력을 가짐. PC2는 21.32%의 설명력을 가짐

#이름이 겹치지 않게끔 만들어주는 패키지 함수
#만약 아웃라이어가 생겼어도 무슨 샘플인지 알 수 있음
library(ggrepel)

ggsave('Fig3.PCAplot.png', width=6, height=4)


### --------------------###
### dendrogram          ###
### --------------------###

#샘플 간의 거리를 계산하여, 각 샘플을 유사도에 따라 군집화
hc=hclust(dist(t(log2DF)), method="average")
head(log2DF)
#행이 유전자, 열이 샘플
#dist는 거리계산을 행 기준으로 함. 
#유전자 간의 각각 거리 보면 너무 많음..
#그래서 트랜스포즈로 행열 바꾸고 샘플 간의 거리를 확인.

#BiocManager::install('ggdendro')
library(ggdendro)

ggdendrogram(hc,rotate=TRUE, leaf_labels = TRUE)
ggsave('Fig4.dendrogram.png', width=4, height = 4) 
#경향성 분석.
#1. 한 조건당 반복수끼리 잘 모임
#2. 같은 셀 타입끼리 잘 묶임. 셀 타입에 의해 구분된 조건들이 더 잘 묶임.

#또 다른 덴드로 그램 함수
BiocManager::install('factoextra')
library(factoextra)
fvix_dend(hc, k=4, rect = TRUE, show_labels=TRUE)

#PCA, dendrogram 모두 다 샘플 QC 과정이다.


### -------------------------###
### expression correlation   ###
### -------------------------###
#### expression correlation 상관관계 분석. -> 후 heatmap 그리기
BiocManager::install('pheatmap')
library(pheatmap)

#pheatmap(cor(log2DF)) 은 아규먼트가 많다. 
#상관계수 관계에 따라 heatmap이 그려짐.
#트리도 그려짐. hcluster(계층적 군집)됨.
#빨강 : 자기 자신과의 상관관계는 1이다. 
#상관계수는 0에서 1.

#색깔 : 저 3개의 색깔을 기준으로 100개의 색을 만들어줌.

pheatmap(cor(log2DF),
         main = 'expression correlation',
         color = colorRampPalette(c('blue','white','red'))(100), #c는 벡터를 만들기 위해.
         treeheight_row = 5,
         treeheight_col = 5,
         cutree_rows = 5, #히트맵을 구역별로 자를 수 있음.
         cutree_cols = 4,  #히트맵을 구역별로 자를 수 있음.
         filename = 'Fig5.correlation.png',  #저장
         width = 5,
         height = 4) 
#상관관계 해석
#CD8에서 상관관계가 더 가까움.
#가장 진한 파란색 : 상관관계가 가장 멈. 셀 타입도 다르고, 넉아웃 넉인 뮤테이션이 들어갔기 때문
#CD8넉인, CD4넉아웃 -> 빨강
#CD8넉아웃, CD4넉인 -> 파랑 -> 더 생물학적으로 많은 차이가 난다는 것을 파악 가능.


### -------------------###
### DEG by cell type   ###
### -------------------###

# 필요한 라이브러리 로드
library(TCC)  # TCC 패키지 로드

# 스텝 1: 데이터 정규화 (TMM 방법 사용)
# 정규화 방법으로 TMM (trimmed mean of M-values)을 사용하여 샘플 간의 변동성을 보정
cellTCC = new('TCC', dataDF, group = as.numeric(factor(metadata$Cell.Type)))
cellTCC = calcNormFactors(cellTCC, norm.method = 'tmm', test.method = 'edger',FDR = 0.05)

# 스텝 2: 차등발현유전자 (DEG) 식별
# edgeR를 사용하여 차등발현유전자 분석
cellTCC = estimateDE(cellTCC, test.method = 'edger', FDR = 0.05)

# DEG 결과를 추출 (로그2 폴드 변화가 2 이상이고, p-value가 0.05 이하인 유전자들)
cellDEG = cellTCC[which(abs(cellTCC$m.value) >= 2 & cellTCC$q.value < 0.05),]
dim(cellDEG)  # DEG로 선정된 유전자 개수 확인

# DEG 유전자들만 추출하여 데이터 준비
log2DF_cellDEG = log2DF[rownames(log2DF) %in% cellDEG$gene_id,]

# 스케일링: 각 유전자의 발현 값을 평균 0, 표준편차 1로 변환
log2DF_cellDEG_scale = t(scale(t(log2DF_cellDEG)))

# pheatmap 그리기
pheatmap(log2DF_cellDEG_scale,
         main = 'DEGs by Cell Type',  # 히트맵 타이틀 설정
         color = colorRampPalette(c('blue', 'white', 'red'))(100),  # 색상(파랑-빨강)
         show_rownames = FALSE,  # 너무 많은 유전자가 있어 행 이름을 숨김
         cluster_cols = FALSE,  # 열(샘플)에 대한 군집화 비활성화
         clustering_method = 'ward.D',  # 유전자 군집화 방법 설정
         treeheight_row = 20,  # 유전자 트리의 높이 설정
         cuttree_rows = 2,  # 유전자 발현 패턴에 따라 2개의 군집으로 나눔
         gaps_col = c(3, 6, 9),  # 샘플별 구분선 위치 설정
         filename = 'Fig6.DEGs_by_Cell_Type.png',  # 파일로 저장
         width = 5, height = 5)


###### DEG - cell type
c(1,1,1,1,2,2,2)
library(TCC)
#스텝1 : 데이터 정규화 : TMM 정규화 방법 : edgeR 패키지에서 제공하는 방법.
#스텝2 : 차등발현유전자(DEG) 식별 
#edgeR : 음이항 분포 기반으로 유전자 발현의 변동성을 모델링.
#일단 셀타입에 따라 CD4면 4, CD8이면 8로 변환됨.
cellTCC = new('TCC', dataDF, group = as.numeric(factor(metadata$Cell.Type)))
cellTCC = calcNormFactors(cellTCC, norm.method = 'tmm', test.method = 'edger', FDR = 0.05) #여기서 method를 deseq, voom으로 바꾸면 해당 계산법대로 계산해줌.
cellTCC = estimateDE(cellTCC, test.method ='edger', FDR = 0.05) #정규화된 카운트 값 넣어, 엣지알로 차등발현 계산해주세요.
#일단은 셀타입기준으로 CD4, CD8의 발현 차이 분석

cellTCC = getResult(cellTCC, sort=F) #sort 안하고, 내가 넣은 유전자 순대로 뱉어달라.
#이게 DEG 구하는 방법. 

cellDEG = cellTCC[abs(cellTCC$m.value) >=2 & cellTCC$q.value <0.05,]
#같은 결과임.
cellDEG = cellTCC[which(abs(cellTCC$m.value) >=2 & cellTCC$q.value <0.05),]
dim(cellDEG) #574   7 발현된것 574개. 
#CD4와 CD8 셀 타입 비교했을 때 차이나는 유전자가 574개이다.
#로그 2 폴드 체인지 값이기에 1 이상이면 2배, 2 이상이면 4배차이.

log2DF_cellDEG = log2DF[rownames(log2DF) %in% cellDEG$gene_id,]

#a.value는 평균
#m.value는 로그 2 폴드 체인지. 우리가 흔히 말하는 차등발현유전자 어쩌구..
#estimate DEG는 1 -> DEG 한것 표시.
#%in% 은 벡터 값이 있는 애들만 true로 반환.
#총 발혇ㄴ된 유전자에서 DEG면 true, 아니면 flase
#이 벡터값을 다시 count 값, 행의 위치에 인덱스로 넣어주면 DEG 값을 데이터로 넣을 수 있음.
dim(log2DF_cellDEG) #574  12 # 왜 7에서 12로 증가하는거지..?

log2DF_cellDEG_scale = t(scale(t(log2DF_cellDEG)))
head(log2DF_cellDEG_scale)
#스케일 맞추는 이유는 pheatmap을 그릴 때 차이를 직관적으로 보기 위해서.
#그 행의 평균값과 표준편차값으로 각각을 스케일 맞춰줌.
#평균을 0으로 맞춤, 표준편차로 맞춤.
#각 유전자 별로 / 행에 있던 유전자를 컬럼으로 위치시키고 스케일링

dev.off()
pheatmap(log2DF_cellDEG_scale,
         main = 'DEGs by Cell Type',
         color = colorRampPalette(c('blue','white','red'))(100),
         show_rownames = FALSE, #row가 너무 많아서 안보이게 함.
         cluster_cols = FALSE,
         clustering_method = 'ward.D', #유전자 끼리 클러스터링.
         treeheight_row = 20,
         cuttree_rows = 2, #유전자 발현 패턴이 비슷한 애들끼리 4개씩 끊김.
         gaps_col = c(3,6,9),
         filename = 'Fig6.DEGs by Cell Type.png',  #저장
         width = 5,
         height = 5) 
#gaps_col : 컬럼 기준으로 3번째, 6번재 9번째에서 끊어줌. #보기편하게.
#행이 유전자라서 더 발현차이가 잘 보임. 그래서 스케일링을 했다. 
dev.off()


### -------------------###
### DEG by mutation    ###
### -------------------###
# 필요한 라이브러리 로드
library(TCC)  # TCC 패키지 로드

# 스텝 1: 데이터 정규화 (TMM 방법 사용)
# 정규화 방법으로 TMM (trimmed mean of M-values)을 사용하여 샘플 간의 변동성을 보정
cellTCC = new('TCC', dataDF, group = as.numeric(factor(metadata$Cell.Type)))
cellTCC = calcNormFactors(cellTCC, norm.method = 'tmm', test.method = 'edger',FDR = 0.05)

# 스텝 2: 차등발현유전자 (DEG) 식별
# edgeR를 사용하여 차등발현유전자 분석
cellTCC = estimateDE(cellTCC, test.method = 'edger', FDR = 0.05)

# DEG 결과를 추출 (로그2 폴드 변화가 2 이상이고, p-value가 0.05 이하인 유전자들)
cellDEG = cellTCC[which(abs(cellTCC$m.value) >= 2 & cellTCC$q.value < 0.05),]
dim(cellDEG)  # DEG로 선정된 유전자 개수 확인

# DEG 유전자들만 추출하여 데이터 준비
log2DF_cellDEG = log2DF[rownames(log2DF) %in% cellDEG$gene_id,]

# 스케일링: 각 유전자의 발현 값을 평균 0, 표준편차 1로 변환
log2DF_cellDEG_scale = t(scale(t(log2DF_cellDEG)))

# pheatmap 그리기
pheatmap(log2DF_cellDEG_scale,
         main = 'DEGs by Cell Type',  # 히트맵 타이틀 설정
         color = colorRampPalette(c('blue', 'white', 'red'))(100),  # 색상 설정 (파랑-빨강)
         show_rownames = FALSE,  # 너무 많은 유전자가 있어 행 이름을 숨김
         cluster_cols = FALSE,  # 열(샘플)에 대한 군집화 비활성화
         clustering_method = 'ward.D',  # 유전자 군집화 방법 설정
         treeheight_row = 20,  # 유전자 트리의 높이 설정
         cuttree_rows = 2,  # 유전자 발현 패턴에 따라 2개의 군집으로 나눔
         gaps_col = c(3, 6, 9),  # 샘플별 구분선 위치 설정
         filename = 'Fig6.DEGs_by_Cell_Type.png',  # 파일로 저장
         width = 5, height = 5)

#DEG - mutation
#메타데이터 그룹 정보를 뮤테이션으로 줬을 뿐이다. 
colnames(dataDF)
metadata$Sample.Name

mutTCC = new('TCC', dataDF, group = as.numeric(factor(metadata$Mutation)))
mutTCC = calcNormFactors(mutTCC, norm.method = 'tmm', test.method = 'edger', FDR = 0.05) 
mutTCC = estimateDE(mutTCC, test.method ='edger', FDR = 0.05) 
mutTCC = getResult(mutTCC, sort=F)
#이게 DEG 구하는 방법. 
mutDEG = mutTCC[abs(mutTCC$m.value) >=2 & mutTCC$q.value <0.05,]
mutDEG = mutTCC[which(abs(mutTCC$m.value) >=2 & mutTCC$q.value <0.05),]
dim(mutDEG) #204   7
log2DF_mutDEG = log2DF[rownames(log2DF) %in% mutDEG$gene_id,]
log2DF_mutDEG_scale = t(scale(t(log2DF_mutDEG)))
head(log2DF_cellDEG_scale)

pheatmap(log2DF_mutDEG_scale,
         main = 'DEGs by Mutation Type',
         color = colorRampPalette(c('green','black','red'))(100),
         show_rownames = FALSE, #row가 너무 많아서 안보이게 함.
         cluster_cols = FALSE,
         clustering_method = 'ward.D', #유전자 끼리 클러스터링.
         treeheight_row = 20,
         cuttree_rows = 2, #유전자 발현 패턴이 비슷한 애들끼리 4개씩 끊김.
         gaps_col = c(3,6,9))
#filename = 'Fig7.DEGs by Cell Type.png',  #저장
#width = 5,
#height = 5) 

dev.off()

### ---------------------------------------------###
### functional enrichment - cell type - CD4_up   ###
### ---------------------------------------------###

#######functional enrichment - cell type

BiocManager::install('clusterProfiler')
BiocManager::install('org.Mm.eg.db')
library(clusterProfiler)
library(org.Mm.eg.db)
#유전자 기능 주석을 활용하여 오믹스 데이터를 해석할 수 있도록 하는 풍부도 분석 도구
#만약에 clusterProfiler에 없는 정보이면 IPA 등을 활용해서 gene annotation해서 기능분석 할 수 있다. 

CD4_up= as.numeric(cellTCC$gene_id[cellTCC$m.value <= -2 & cellTCC$q.value <=0.05])
#CD4에서 up인 것을 뽑음.
#왜 -2 이하일까? 
#컨트롤 대비 케이스에서 발현이 높은것을 로그 폴드 체인지라고 함. 
#logFC 절댓값 <= 1 -> 컨트롤 vs 케이스 에서 차이나는 것
#로그 폴드 체인지 값이 1이상 : 케이스에서 up regulate 된 DEG 보기 
#logFC 값이 -1이하 : 케이스에서 down regulate 된 DEG 보기
#CD8 vs CD4에서음수인것은 CD8에서 발현이 낮다 = CD4에서 발현이 up되었다. 
#CD8 vs CD4에서양수인것은 CD8에서 발현이 높다
#CD4를 컨트롤로 잡은 것임. 이걸 어떻게 알 수 있냐? 어디서 조건 설정? 
#아까 cell type heatmap 그릴 때 메타데이터에 있는 cell.type을 팩터로, 숫자형으로 줌.
#메타데이터 상 순서(레벨)가 CD4먼저임. CD4가 1이됨. 1이 컨트롤이됨. / CD8은 2로 케이스가 됨.
#로그폴드 양수이면 케이스에서 up regulate된것이기에 케이스=CD8에서 up된것.
#즉, 컨트롤, 케이스를 바꾸려면 메타데이터 순서를 바꾸거나, 조건을 바꾸거나. 정렬을 바꾸거나 하면 됨.

length(CD4_up) #267개

CD4_enrich = enrichGO(CD4_up, 'org.Mm.eg.db', ont='BP', keyType='ENTREZID', minGSSize=5, maxGSSize=200) #과발현된 유전자
class(CD4_enrich)
CD4_enrich = CD4_enrich@result #결과 파일 불러오기
CD4_enrich = CD4_enrich[CD4_enrich$qvalue <= 0.05,] #유의수준 0.05 #q-value를 기준으로 데이터 QC
CD4_enrich = CD4_enrich[order(CD4_enrich$Count, decreasing = c(T)),] #내림차순으로 유전자 정렬 #count값으로 sorting 함. 
head(CD4_enrich)

##### 그래프 그리기
CD4_enrich$ID = factor(CD4_enrich$ID, levels = CD4_enrich$ID) #이거 해야 내가 넣은 순서대로 볼 수 있음. 
head(CD4_enrich)

ggplot(CD4_enrich[1:10,])+ #10개에 대한 GO에 대해서 값이 나오게 함. 
  geom_point(aes(y=ID, x=GeneRatio, size=Count, color=qvalue))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

dev.off()

ggsave('Fig8.CD4-GOenrich.png', width=5, height=4) #그래프 저장
write.table(CD4_enrich, 'Table.CD4-GOenrich.txt', sep='\t', quote = FALSE, col.names = TRUE)

### -----------------------------------------------###
### functional enrichment - cell type - CD4_down   ###
### -----------------------------------------------###

#######functional enrichment - cell type - CD4_down


CD4_down= as.numeric(cellTCC$gene_id[cellTCC$m.value >= 2 & cellTCC$q.value <=0.05])
length(CD4_down) #307

CD4_enrich_down = enrichGO(CD4_down, 'org.Mm.eg.db', ont='BP', keyType='ENTREZID', minGSSize=5, maxGSSize=200) #과발현된 유전자
class(CD4_enrich_down)
CD4_enrich_down = CD4_enrich_down@result #결과 파일 불러오기
CD4_enrich_down = CD4_enrich_down[CD4_enrich_down$qvalue <= 0.05,] #유의수준 0.05 #q-value를 기준으로 데이터 QC
CD4_enrich_down = CD4_enrich_down[order(CD4_enrich_down$Count, decreasing = c(T)),] #내림차순으로 유전자 정렬 #count값으로 sorting 함. 
head(CD4_enrich_down)

##### 그래프 그리기
CD4_enrich_down$ID = factor(CD4_enrich_down$ID, levels = CD4_enrich_down$ID) #이거 해야 내가 넣은 순서대로 볼 수 있음. 
head(CD4_enrich_down)

ggplot(CD4_enrich_down[1:10,])+ #10개에 대한 GO에 대해서 값이 나오게 함. 
  geom_point(aes(y=ID, x=GeneRatio, size=Count, color=qvalue))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

dev.off()

ggsave('Fig11.CD4_down-GOenrich.png', width=5, height=4) #그래프 저장
write.table(CD4_enrich, 'Table.CD4_down-GOenrich.txt', sep='\t', quote = FALSE, col.names = TRUE)


### ---------------------------------------------###
### functional enrichment - DKO up regulate      ###
### ---------------------------------------------###
library(clusterProfiler)
library(org.Mm.eg.db)
DKO_up= as.numeric(mutTCC$gene_id[mutTCC$m.value <= -2 & mutTCC$q.value <=0.05])
#CD4 = DKO = 컨트롤 그룹
#CD8 = Kin = 케이스 그룹 임을 메타데이터 순서를 통해 확인
length(DKO_up) ##99개

DKO_enrich = enrichGO(DKO_up, 'org.Mm.eg.db', ont='BP', keyType='ENTREZID', minGSSize=5, maxGSSize=200) #과발현된 유전자
class(DKO_enrich)
DKO_enrich = DKO_enrich@result #결과 파일 불러오기
#유의수준 0.05 #q-value를 기준으로 데이터 QC
DKO_enrich = DKO_enrich[DKO_enrich$qvalue <= 0.05,] 
#내림차순으로 유전자 정렬 #count값으로 sorting 함. 
DKO_enrich = DKO_enrich[order(DKO_enrich$Count, decreasing = c(T)),] 
head(DKO_enrich)

##### 그래프 그리기
#이거 해야 내가 넣은 순서대로 볼 수 있음. 
DKO_enrich$ID = factor(DKO_enrich$ID, levels = DKO_enrich$ID) 
head(DKO_enrich)

ggplot(DKO_enrich[1:10,])+ #10개에 대한 GO에 대해서 값이 나오게 함. 
  geom_point(aes(y=ID, x=GeneRatio, size=Count, color=qvalue))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

dev.off()

ggsave('Fig9.DKO_up-GOenrich.png', width=5, height=4) #그래프 저장
write.table(CD4_enrich, 'Table.DKO_up-GOenrich.txt', sep='\t', quote = FALSE, col.names = TRUE)




### ---------------------------------------------###
### functional enrichment - DKO down regulate    ###
### ---------------------------------------------###
################## DKO에서 down regulated 유전자 추출
# DKO에서 down-regulated 유전자 추출
# m.value >= 2: 로그2 폴드 변화가 2 이상인 경우 -> 케이스에서 발현이 낮아짐 (down-regulated)
DKO_down = as.numeric(mutTCC$gene_id[mutTCC$m.value >= 2 & mutTCC$q.value <= 0.05])

length(DKO_down)  # down-regulated 유전자 개수 확인 ### 105

# GO enrichment 분석
DKO_enrich_down = enrichGO(DKO_down, 'org.Mm.eg.db', ont = 'BP', 
                           keyType = 'ENTREZID', minGSSize = 5, maxGSSize = 200)

# GO enrichment 결과 확인
class(DKO_enrich_down)
DKO_enrich_down = DKO_enrich_down@result  # 결과 파일 불러오기

# 유의수준 0.05로 데이터 QC
DKO_enrich_down = DKO_enrich_down[DKO_enrich_down$qvalue <= 0.05,]  # q-value 기준 필터링

# Count 값으로 내림차순 정렬
DKO_enrich_down = DKO_enrich_down[order(DKO_enrich_down$Count, decreasing = TRUE),] 
head(DKO_enrich_down)

##### 그래프 그리기
# GO ID 순서대로 표시
DKO_enrich_down$ID = factor(DKO_enrich_down$ID, levels = DKO_enrich_down$ID)
head(DKO_enrich_down)

# 상위 10개 GO 항목에 대한 산점도 그래프
ggplot(DKO_enrich_down[1:10,]) +
  geom_point(aes(y = ID, x = GeneRatio, size = Count, color = qvalue)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

dev.off()

# 그래프 저장
ggsave('Fig10.DKO_down-GOenrich.png', width = 5, height = 4)

# 결과를 텍스트 파일로 저장
write.table(DKO_enrich_down, 'Table.DKO_down-GOenrich.txt', sep = '\t', quote = FALSE, col.names = TRUE)

# q-value가 0.05 이하인 GO ID와 description 추출
DKO_enrich_down_significant = DKO_enrich_down[DKO_enrich_down$qvalue <= 0.05,]

# q-value가 작은 순서대로 정렬
DKO_enrich_down_significant_sorted = DKO_enrich_down_significant[order(DKO_enrich_down_significant$qvalue),]

# GO ID와 Description만 추출
GO_description = DKO_enrich_down_significant_sorted[, c("ID", "Description")]

# 결과 출력
print(GO_description)

# 추가로, GO ID와 관련된 유전자 수가 5개 이상인 항목만 필터링
DKO_enrich_down_filtered = DKO_enrich_down_significant_sorted[DKO_enrich_down_significant_sorted$Count >= 5,]

# 필터링된 GO ID와 Description 출력
GO_filtered_description = DKO_enrich_down_filtered[, c("ID", "Description")]
print(GO_filtered_description)
