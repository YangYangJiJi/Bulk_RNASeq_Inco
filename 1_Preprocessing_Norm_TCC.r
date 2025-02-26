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