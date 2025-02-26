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