import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker
import numpy as np
from scipy.stats import ttest_ind, shapiro, mannwhitneyu
import itertools

# 색상 설정
colors = {
    "CONTROL": (0.18, 0.55, 0.34),
    "HBSS": (1.0, 0.4, 0.0),
    "HBSS+HN": (0.0, 0.3, 1.0)
}

# 데이터 불러오기 (이상치 제거 없음)
start_df = pd.read_csv("/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity/metagene_output/merged_start_codon_metagene.csv")  
stop_df = pd.read_csv("/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity/metagene_output/merged_stop_codon_metagene.csv")    

# Violin + Boxplot (이상치 제거 X)
def plot_violin_boxplot(df, title, y_label, filename, ylim=None):
    plt.figure(figsize=(6, 6), dpi=600)
    
    # Violin plot
    sns.violinplot(data=df, x="Group", y="RPM", hue="Group", palette=colors, inner=None, alpha=0.7, legend=False)
    # Boxplot (이상치 포함)
    sns.boxplot(data=df, x="Group", y="RPM", hue="Group", palette=colors, width=0.3, showcaps=False, showfliers=True, whiskerprops=dict(linestyle="-"), legend=False)
    # 개별 데이터 포인트 표시
    sns.stripplot(data=df, x="Group", y="RPM", color="black", size=1.5, jitter=True, alpha=0.6)

    plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.3e'))
    plt.title(title, fontsize=16, fontweight="bold")
    plt.ylabel(y_label, fontsize=10)
    plt.xlabel("Condition", fontsize=10)
    
    if ylim:
        plt.ylim(ylim)

    plt.savefig(filename, bbox_inches="tight")
    plt.show()

# 실행 및 저장 (이상치 제거 안함)
plot_violin_boxplot(start_df, "Around start codon", "P-site count", "metagene_output/start_violin_no_outlier.png")
plot_violin_boxplot(stop_df, "Around stop codon", "P-site count", "metagene_output/stop_violin_no_outlier.png")


# 통계 비교 수행 및 추가 분석 함수
def compare_groups(df, output_file):
    groups = df["Group"].unique()
    comparisons = list(itertools.combinations(groups, 2))
    results = []

    for g1, g2 in comparisons:
        data1 = df[df["Group"] == g1]["RPM"]
        data2 = df[df["Group"] == g2]["RPM"]

        p1 = shapiro(data1)[1]
        p2 = shapiro(data2)[1]

        if p1 > 0.05 and p2 > 0.05:
            stat, pval = ttest_ind(data1, data2)
            test_used = "t-test"
        else:
            stat, pval = mannwhitneyu(data1, data2)
            test_used = "Mann-Whitney"
        
        # 추가 통계량 계산 (IQR, 평균, 표준편차 등)
        iqr1 = np.percentile(data1, 75) - np.percentile(data1, 25)
        iqr2 = np.percentile(data2, 75) - np.percentile(data2, 25)
        mean1, std1 = data1.mean(), data1.std()
        mean2, std2 = data2.mean(), data2.std()

        results.append([
            g1, g2, f"{stat:.3e}", f"{pval:.3e}", test_used,
            f"{iqr1:.3e}", f"{mean1:.3e}", f"{std1:.3e}",
            f"{iqr2:.3e}", f"{mean2:.3e}", f"{std2:.3e}"
        ])
    
    # 데이터프레임 생성 및 CSV 저장
    results_df = pd.DataFrame(results, columns=[
        "Group1", "Group2", "Statistic", "P-value", "Test",
        "IQR_G1", "Mean_G1", "Std_G1", "IQR_G2", "Mean_G2", "Std_G2"
    ])
    results_df.to_csv(output_file, index=False)

    return results_df

# 실행 및 저장
start_stats_df = compare_groups(start_df, "metagene_output/start_violin_no_outlier.csv")
stop_stats_df = compare_groups(stop_df, "metagene_output/stop_violin_no_outlier.csv")

