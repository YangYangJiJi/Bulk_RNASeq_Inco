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