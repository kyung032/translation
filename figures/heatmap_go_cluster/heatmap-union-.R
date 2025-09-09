

###############################################
# 1. 필요한 패키지 로드
###############################################
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(MASS)
library(cluster)

###############################################
# 2. RNA-seq & Ribo-seq 정규화 데이터 로드
###############################################

# RNA count 데이터
rna_g2_vs_g1 <- read.table("DESeq2_output_RNA/DEG_Norm_G2vsG1_RNA.txt", header=TRUE, sep="\t", row.names=1)
rna_g3_vs_g2 <- read.table("DESeq2_output_RNA/DEG_Norm_G3vsG2_RNA.txt", header=TRUE, sep="\t", row.names=1)

# RIBO count 데이터
ribo_g2_vs_g1 <- read.table("DESeq2_output_RIBO/DEG_Norm_G2vsG1_RIBO.txt", header=TRUE, sep="\t", row.names=1)
ribo_g3_vs_g2 <- read.table("DESeq2_output_RIBO/DEG_Norm_G3vsG2_RIBO.txt", header=TRUE, sep="\t", row.names=1)

# DEG 리스트 불러오기 (하나도 변는 유전)

up_G2vsG1 <- read.table("DESeq2_TE_upregulated_G2vsG1.txt", sep ="\t", quote="", fill =TRUE, header=TRUE)
up_G3vsG2 <- read.table("DESeq2_TE_upregulated_G3vsG2.txt", sep ="\t", quote="", fill =TRUE, header=TRUE)
down_G2vsG1 <- read.table("DESeq2_TE_downregulated_G2vsG1.txt", sep ="\t", quote="", fill =TRUE, header=TRUE)
down_G3vsG2 <- read.table("DESeq2_TE_downregulated_G3vsG2.txt", sep ="\t", quote="", fill =TRUE, header=TRUE)

colnames(up_G2vsG1)


union_genes <- Reduce(union, list(
  up_G2vsG1$GeneID,
  down_G2vsG1$GeneID,
  up_G3vsG2$GeneID,
  down_G3vsG2$GeneID
))

cat("📌 합집합 DEG 유전자 개수:", length(union_genes), "\n")
# 📌 합집합 DEG 유전자 개수: 2980

# 공통 유전자 필터링
rna_g2_vs_g1 <- rna_g2_vs_g1[union_genes, ]
rna_g3_vs_g2 <- rna_g3_vs_g2[union_genes, ]
ribo_g2_vs_g1 <- ribo_g2_vs_g1[union_genes, ]
ribo_g3_vs_g2 <- ribo_g3_vs_g2[union_genes, ]


###############################################
# 3. G2 샘플 평균값 계산
###############################################

# G2 샘플의 평균값 계산
rna_g2_avg <- (rna_g2_vs_g1[, c("G2_1_RNA", "G2_2_RNA")] + rna_g3_vs_g2[, c("G2_1_RNA", "G2_2_RNA")]) / 2
ribo_g2_avg <- (ribo_g2_vs_g1[, c("G2_1_RIBO", "G2_2_RIBO")] + ribo_g3_vs_g2[, c("G2_1_RIBO", "G2_2_RIBO")]) / 2

# G1, G3 샘플 포함하여 새로운 매트릭스 생성
rna_counts <- cbind(rna_g2_vs_g1[, c("G1_1_RNA", "G1_2_RNA")], rna_g2_avg, rna_g3_vs_g2[, c("G3_1_RNA", "G3_2_RNA")])
ribo_counts <- cbind(ribo_g2_vs_g1[, c("G1_1_RIBO", "G1_2_RIBO")], ribo_g2_avg, ribo_g3_vs_g2[, c("G3_1_RIBO", "G3_2_RIBO")])

# 새로운 샘플 이름 설정
colnames(rna_counts) <- c("G1_1", "G1_2", "G2_avg1", "G2_avg2", "G3_1", "G3_2")
colnames(ribo_counts) <- c("G1_1", "G1_2", "G2_avg1", "G2_avg2", "G3_1", "G3_2")


###############################################
# 4. Translational Efficiency (TE) 계산
###############################################

# TE = log2(Ribo-seq / RNA-seq)
te_matrix <- log2((ribo_counts + 1) / (rna_counts + 1))  # 1을 더해 0값 방지

# Z-score 변환
te_matrix_z <- t(scale(t(te_matrix), center = TRUE, scale = TRUE))

# 결과 확인
cat("📊 TE 매트릭스 크기:", dim(te_matrix_z), "\n")
#📊 TE 매트릭스 크기: 2980 6 

write.csv(te_matrix_z, file="TE_Z_Score_union_2980_DEG.csv")


###############################################
# 1️⃣ Gaussian Mixture Model (GMM) 클러스터링
###############################################
library(mclust)

# 📌 BIC 값 저장을 위해 초기화
bic_values <- numeric(length = 14)  # G = 2~15 (14개)
names(bic_values) <- 2:15  # 클러스터 개수(G)와 매칭

# 📌 BIC 값 계산 (G = 2 ~ 15)
for (k in 2:15) {
  gmm_model <- tryCatch(
    Mclust(te_matrix_z, G = k),  # GMM 실행
    error = function(e) return(NULL)  # 오류 발생 시 NULL 반환
  )
  
  if (!is.null(gmm_model)) {
    bic_values[as.character(k)] <- gmm_model$BIC  # BIC 저장
  } else {
    bic_values[as.character(k)] <- NA  # 오류 발생 시 NA 저장
  }
}

# 📌 유효한 BIC 값만 필터링
valid_idx <- which(!is.na(bic_values))
bic_values <- bic_values[valid_idx]  
valid_G <- as.numeric(names(bic_values))  # 클러스터 개수를 숫자로 변환

# 📌 BIC 값이 없을 경우 종료
if (length(valid_G) == 0) {
  stop("⚠️ GMM 실행 불가: 유효한 BIC 값이 없음. 데이터 확인 필요.")
}

# 📌 BIC 값 플로팅 (Elbow Point 찾기)
plot(valid_G, bic_values, type="b", pch=19, col="blue",
     xlab="Number of Clusters (G)", ylab="BIC", main="BIC vs G")

### 그래프 저장하여 elbow point: BIC값이 증가했다가 완만해지는 점 표현
#최적 G는 그래프항 11
#현재 제공된 그래프를 볼 때, 대략 11개에서 급격한 증가세가 완만해지는 경향이 뚜렷합니다. 따라서 최적의 클러스터 개수는 11개

# 📌 최적 G로 GMM 실행
gmm_res <- Mclust(te_matrix_z, G = 11)

# 📌 클러스터 개수 확인
cat("📌 GMM 클러스터 개수 확인:\n")
print(table(gmm_res$classification))

#1   2   3   4   5   6   7   8   9  10  11 
#284 224 219 558 160 254 271 230 224 398 158 

# GMM 클러스터를 TE matrix 정보 추가
# 행렬을 데이터 프레임으로 변환하면서 Cluster 정보를 추가


# GMM 클러스터 결과를 데이터 프레임으로 변환
te_matrix_z_clustered <- as.data.frame(te_matrix_z)  

# Cluster 열 추가
te_matrix_z_clustered$Cluster <- factor(gmm_res$classification)

# 📌 확인
head(te_matrix_z_clustered)
str(te_matrix_z_clustered)


# 다시 구조 확인
str(te_matrix_z_clustered)

#G1_1       G1_2    G2_avg1    G2_avg2       G3_1       G3_2
#ENSG00000000419 -0.9869972 -0.5616157  0.9270010  1.5539667 -0.4947970 -0.4375579
#ENSG00000003756 -0.7986002 -0.4118099  1.4239557  1.0988638 -0.8756243 -0.4367851
#ENSG00000004534 -0.7084493  0.1249279  0.6131771  1.6432470 -0.8953752 -0.7775275
#ENSG00000004779 -0.3302416  0.2214062  0.8261864  1.3887324 -0.9725288 -1.1335546
#ENSG00000004864  0.6613362  0.7728327 -1.6865442 -0.7483091  0.3267313  0.6739531
#ENSG00000004897  0.6238204  0.4760755 -1.7344035 -0.6924670  0.5737564  0.7532183
#Cluster
#ENSG00000000419       1
#ENSG00000003756       2
#ENSG00000004534       9
#ENSG00000004779       3
#ENSG00000004864       5
#ENSG00000004897       5


###############################################
# 2️⃣ 클러스터별 GO term enrichment 분석
###############################################

library(clusterProfiler)
library(org.Hs.eg.db)

# 모든 클러스터 저장
all_clusters <- unique(te_matrix_z_clustered$Cluster)

# 클러스터별 GO term 분석
enrichment_results <- list()
all_go_results <- list()

for (i in all_clusters) {
  genes_in_cluster <- rownames(te_matrix_z_clustered)[te_matrix_z_clustered$Cluster == i]
  
  # GO enrichment
  ego <- enrichGO(gene = genes_in_cluster, OrgDb = org.Hs.eg.db, 
                  keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05)
  
  enrichment_results[[as.character(i)]] <- ego
  
  # 원본 GO term 저장
  if (!is.null(ego) && nrow(ego@result) > 0) {
    all_go_results[[as.character(i)]] <- ego@result
  } else {
    all_go_results[[as.character(i)]] <- data.frame(
      ID = NA, Description = "No significant GO term",
      Cluster = i, pvalue = NA, padj = NA, GeneRatio = NA, BgRatio = NA
    )
  }
}

# ✅ 원본 GO terms 데이터 저장
all_go_results_df <- bind_rows(all_go_results, .id = "Cluster")
write.csv(all_go_results_df, "2980_DEG_All_GO_terms_heatmap.csv", row.names = FALSE)

# ✅ 최종 논문 스타일 요약본 저장 (p-value 기준 상위 5개 GO Term)
go_summary <- lapply(enrichment_results, function(ego) {
  if (!is.null(ego) && nrow(ego@result) > 0) {
    ego@result %>%
      filter(pvalue < 0.05) %>%
      arrange(pvalue) %>%
      slice_head(n = 5)
  } else {
    data.frame(Cluster = NA, Description = "No significant GO term")
  }
})

go_summary_df <- bind_rows(go_summary, .id = "Cluster")
write.csv(go_summary_df, "2980_DEG_All_GO_terms_heatmap_parentGO.csv", row.names = FALSE)

cat("✅ GO Term 분석 완료!\n")


###############################################
# 3 HEATMAP
###############################################

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 📌 클러스터별 색상 매핑 (Named Vector)
cluster_levels <- levels(te_matrix_z_clustered$Cluster)
cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(cluster_levels)), cluster_levels)

# 📌 클러스터별 간격을 조정한 row annotation 생성
annotation_row <- data.frame(Cluster = factor(te_matrix_z_clustered$Cluster))

row_ha <- rowAnnotation(
  Cluster = annotation_row$Cluster, 
  col = list(Cluster = cluster_colors)  # Named vector 사용
)

# ✅ 정상 작동 여부 확인
print(cluster_colors)
#1         2         3         4         5         6         7         8 
#"#8DD3C7" "#F8F8B6" "#CAAEC5" "#D68E8F" "#B2B2A5" "#D8C965" "#DED3B3" "#E3D5DC" 
#9        10        11 
#"#C191C2" "#CAE0C4" "#FFED6F" 

print(row_ha)


# 📌 논문 스타일 히트맵 생성
ht <- Heatmap(
  as.matrix(te_matrix_z_clustered[, 1:6]),  # 유전자 발현 값 (Z-score 변환됨)
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),  # 색상 설정 (논문 스타일)
  
  # 📌 클러스터별 그룹화
  row_split = annotation_row$Cluster,
  cluster_rows = TRUE,  
  cluster_columns = TRUE,  
  show_row_names = FALSE,  
  show_column_names = TRUE,  
  row_title_gp = gpar(fontsize = 10, fontfamily = "Arial"),  
  column_title_gp = gpar(fontsize = 10, fontfamily = "Arial"),  
  
  # 📌 클러스터별 색상 표시
  left_annotation = row_ha
)

# ✅ TIFF 형식으로 저장
tiff("2980_GO_Terms_Heatmap.tiff", width = 8, height = 10, units = "in", res = 600)
draw(ht)
dev.off()

cat("✅ 히트맵 생성 완료!\n")




