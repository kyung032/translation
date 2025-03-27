setwd("C:/Users/WIN/Desktop/science-2024-submit/5번")

###############################################
# 1️⃣ 필요한 라이브러리 로드
###############################################
library(ggplot2)
library(ggrepel)
library(dplyr)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")


library(org.Hs.eg.db)  # Gene Symbol 매핑

###############################################
# 2️⃣ 데이터 로드 (G2vsG1 및 G3vsG2)
###############################################

# 🔹 G2 vs G1
resDF_G2vsG1 <- read.table("DESeq2_TE_results_G2vsG1.txt", header=TRUE, sep="\t")

# 🔹 G3 vs G2
resDF_G3vsG2 <- read.table("DESeq2_TE_results_G3vsG2.txt", header=TRUE, sep="\t")

###############################################
# 3️⃣ -log10(padj) 계산 (FDR 변환)
###############################################
resDF_G2vsG1$negLog10Padj <- -log10(ifelse(resDF_G2vsG1$padj == 0, 1e-10, resDF_G2vsG1$padj))
resDF_G3vsG2$negLog10Padj <- -log10(ifelse(resDF_G3vsG2$padj == 0, 1e-10, resDF_G3vsG2$padj))

###############################################
# 4️⃣ 유의성 그룹 지정 (Up/Down/Not Sig)
###############################################

# 🔹 G2 vs G1
resDF_G2vsG1$Significance <- "Not Sig"
resDF_G2vsG1$Significance[resDF_G2vsG1$padj < 0.05 & resDF_G2vsG1$log2FoldChange >= 1] <- "Up"
resDF_G2vsG1$Significance[resDF_G2vsG1$padj < 0.05 & resDF_G2vsG1$log2FoldChange <= -1] <- "Down"

# 🔹 G3 vs G2
resDF_G3vsG2$Significance <- "Not Sig"
resDF_G3vsG2$Significance[resDF_G3vsG2$padj < 0.05 & resDF_G3vsG2$log2FoldChange >= 1] <- "Up"
resDF_G3vsG2$Significance[resDF_G3vsG2$padj < 0.05 & resDF_G3vsG2$log2FoldChange <= -1] <- "Down"

###############################################
# 5️⃣ Gene Symbol 매핑 (ENSG → 공식 심볼)
###############################################

# GeneID가 ENSG 형식이면 Symbol로 변환
if (!"Symbol" %in% colnames(resDF_G2vsG1)) {
  resDF_G2vsG1$Symbol <- mapIds(org.Hs.eg.db, 
                                keys = resDF_G2vsG1$GeneID, 
                                column = "SYMBOL", 
                                keytype = "ENSEMBL", 
                                multiVals = "first")
}

if (!"Symbol" %in% colnames(resDF_G3vsG2)) {
  resDF_G3vsG2$Symbol <- mapIds(org.Hs.eg.db, 
                                keys = resDF_G3vsG2$GeneID, 
                                column = "SYMBOL", 
                                keytype = "ENSEMBL", 
                                multiVals = "first")
}

###############################################
# 6️⃣ Volcano Plot 생성 함수 (Top 10 유전자 라벨 추가)
###############################################

plot_volcano <- function(data, title, filename) {
  # 상위 10개 Up 및 Down 유전자 (padj 기준으로 정렬)
  top_up <- data %>%
    filter(Significance == "Up") %>%
    arrange(padj) %>%
    head(10)
  
  top_down <- data %>%
    filter(Significance == "Down") %>%
    arrange(padj) %>%
    head(10)
  
  # Volcano Plot 생성
  volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = negLog10Padj, color = Significance)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("Up" = "#CC0033", "Down" = "#0066CC", "Not Sig" = "grey")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(
      title = title,
      x = expression(log[2]~"(Fold change)"),
      y = expression(-log[10]~"(FDR)")
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 10, margin = margin(t = 10)),
      axis.title.y = element_text(size = 10, margin = margin(r = 10)),
      axis.text = element_text(size = 10)
    ) +
    # 유전자 이름 라벨 추가 (상위 10개)
    geom_label_repel(data = top_up, aes(label = Symbol), 
                     color = "#CC0033", fill = "white", box.padding = 0.7, 
                     point.padding = 0.2, segment.color = "#CC0033",  segment.alpha = 1, size = 1.5, nudge_x = 2, max.overlaps = Inf, force = 15) +
    geom_label_repel(data = top_down, aes(label = Symbol), 
                     color = "#0066CC", fill = "white", box.padding = 0.7, 
                     point.padding = 0.2, segment.color = "#0066CC", segment.alpha = 1, size = 1.5, nudge_x = -5, max.overlaps = Inf, force = 15)
  
  #nudge > 축 평행이동, max.overlaps = Inf 겹쳐도 무조건 표시, box.padding > 박스 최소간격, point > 점 박스 최소간격
  
  # 플롯 저장
  ggsave(filename = filename, plot = volcano_plot, width = 12, height = 10, units = "cm", dpi = 600)
  
  # 플롯 출력
  print(volcano_plot)
}

###############################################
# 7️⃣ G2vsG1, G3vsG2 각각 Volcano Plot 생성
###############################################

# 🔹 G2 vs G1 TE 변화
plot_volcano(resDF_G2vsG1, "Translational Efficiency (G2 vs G1)", "volcano_plot_G2vsG1_practice7.png")

# 🔹 G3 vs G2 TE 변화
plot_volcano(resDF_G3vsG2, "Translational Efficiency (G3 vs G2)", "volcano_plot_G3vsG2_practice.png")

cat("✅ Volcano plot 저장 완료!\n")
