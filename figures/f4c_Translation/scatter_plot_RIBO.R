################################################################
#          1. scatter plot              #
################################################################
setwd("G:/riboseq_2025/riboseq_analysis/scatter_plot")

# scatter plot 을 그리기 위해 필요한 라이브러리
library(ggplot2)

################################################################
#          1) DATA INPUT              #
################################################################
g1_g2 = read.table("DEG_G1vsG2_RIBO.txt", header =T)
head(g1_g2)
g1_g2_fc=g1_g2$log2.FC.
g1_g2_pval=g1_g2$P.adj

g2_g3 = read.table("DEG_G2vsG3_RIBO.txt", header =T)
head(g2_g3)
g2_g3_fc=g2_g3$log2.FC.
g2_g3_pval=g2_g3$P.adj


# combined_data$color <- ifelse(combined_data$pvalue_G1_G2 < 0.05 & abs(combined_data$log2FoldChange_G1_G2) > 1, "red", "black")


# g1_g2와 g2_g3 데이터를 합침
volcano_data <- data.frame(
  G1_G2_fc = g1_g2_fc,
  G2_G3_fc = g2_g3_fc,
  G1_G2_padj = g1_g2_pval,
  G2_G3_padj = g2_g3_pval
)

################################################################
#          2) DATA ANALYSIS              #
################################################################

# 조건에 만족하는 데이터 필터링
significant_points <- subset(volcano_data, abs(G1_G2_fc) >= 1 & G1_G2_padj < 0.05)
# 조건에 만족하는 데이터 개수 출력
nrow(significant_points)
# 2597

# 조건에 만족하는 데이터 필터링
significant_points_HN <- subset(volcano_data, abs(G2_G3_fc) >= 1 & G2_G3_padj < 0.05)
# 조건에 만족하는 데이터 개수 출력
nrow(significant_points_HN)
# 2264


# 1 이상인 경우
significant_up <- subset(volcano_data, G1_G2_fc >= 1 & G1_G2_padj < 0.05)

# -1 이하인 경우
significant_down <- subset(volcano_data, G1_G2_fc <= -1 & G1_G2_padj < 0.05)

# 결과 출력
cat("1 이상인 데이터 개수:", nrow(significant_up), "\n")
# 934
cat("-1 이하인 데이터 개수:", nrow(significant_down), "\n")
# 1663

# 1 이상인 경우
significant_up_HN <- subset(volcano_data, G2_G3_fc >= 1 & G2_G3_padj < 0.05)

# -1 이하인 경우
significant_down_HN <- subset(volcano_data, G2_G3_fc <= -1 & G2_G3_padj < 0.05)

# 결과 출력
cat("1 이상인 데이터 개수:", nrow(significant_up_HN), "\n")
# 718
cat("-1 이하인 데이터 개수:", nrow(significant_down_HN), "\n")
# 1546 

# NA 값 제거
# NA를 1로 대체
volcano_data$G1_G2_padj[is.na(volcano_data$G1_G2_padj)] <- 1
volcano_data$G2_G3_padj[is.na(volcano_data$G2_G3_padj)] <- 1



################################################################
#          3) SCATTER PLOT 생성         #
################################################################

# Volcano plot 생성
plot <- ggplot(volcano_data, aes(x = G1_G2_fc, y = G2_G3_fc)) +
  # 빨간색 점
  geom_point(data = subset(volcano_data, abs(G1_G2_fc) >= 1 & G1_G2_padj < 0.05),
             color = "firebrick", size = 0.7) +
  # 검은색 점
  geom_point(data = subset(volcano_data, !(abs(G1_G2_fc) >= 1 & G1_G2_padj < 0.05)),
             color = "dimgray", size = 0.3) +
  ggtitle("Translation") +
  xlab("Nutrient deprivation") +
  ylab("Humanin treatment") +
  # 축 범위를 강제적으로 설정
  scale_x_continuous(limits = c(-10, 10)) +  # X축 범위 설정
  scale_y_continuous(limits = c(-10, 10)) +  # Y축 범위 설정
  coord_fixed(ratio = 1) +  # 축 비율 고정
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# 동일한 크기(6 cm × 6 cm)로 저장
ggsave("scatter_plot_RIBO.tiff", plot = plot, device = "tiff",
       width = 7, height = 7, units = "cm", dpi = 600)


