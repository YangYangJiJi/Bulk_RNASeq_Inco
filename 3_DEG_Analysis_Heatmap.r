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