import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr

# 데이터 로드
file1 = "Merged_polarity_stalling_850_G2vsG1_TE_Up_G3vsG2_TE_Down.csv"
file2 = "Merged_polarity_stalling_835_G2vsG1_TE_Down_G3vsG2_TE_Up.csv"

df_up_g2_down_g3 = pd.read_csv(file1)
df_down_g2_up_g3 = pd.read_csv(file2)


# RGB 색상 정규화
colors = {'up': np.array([255, 36, 0])/255, 'down': np.array([65, 105, 225])/255}

files_info = [
    ("Merged_polarity_stalling_850_G2vsG1_TE_Up_G3vsG2_TE_Down.csv", "G2 up, G3 down"),
    ("Merged_polarity_stalling_835_G2vsG1_TE_Down_G3vsG2_TE_Up.csv", "G2 down, G3 up")
]


# 상관관계 분석 및 저장 (과학적 표기법 적용)
for file, label in files_info:
    df = pd.read_csv(file)

    plt.figure(figsize=(12,9))
    bins = np.linspace(-1, 1, 60)

    if label == "G2 up, G3 down":
        plt.hist(df["Stalling_Change_G2vsG1"], bins=bins, color=colors['up'], alpha=0.9,
                 edgecolor='black', linewidth=0.7, label='G2 vs G1 (up)')
        plt.hist(df["Stalling_Change_G3vsG2"], bins=bins, color=colors['down'], alpha=0.9,
                 edgecolor='black', linewidth=0.7, label='G3 vs G2 (down)',
                 weights=-np.ones_like(df["Stalling_Change_G3vsG2"]))
    else:  # "G2 down, G3 up"
        plt.hist(df["Stalling_Change_G3vsG2"], bins=bins, color=colors['up'], alpha=0.9,
                 edgecolor='black', linewidth=0.7, label='G3 vs G2 (up)')
        plt.hist(df["Stalling_Change_G2vsG1"], bins=bins, color=colors['down'], alpha=0.9,
                 edgecolor='black', linewidth=0.7, label='G2 vs G1 (down)',
                 weights=-np.ones_like(df["Stalling_Change_G2vsG1"]))

    plt.xlabel("Stalling Index Change", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.title(f"Enhanced Visualization of Stalling Index Change - {label}", fontsize=16)
    plt.axhline(0, color='black', linewidth=1.2)
    plt.legend(fontsize=12)
    plt.grid(alpha=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"Enhanced_Stalling_Change_{label.replace(', ', '_').replace(' ', '')}.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Pearson 상관관계 및 p-value 계산
    corr, pval = pearsonr(df["Stalling_Change_G2vsG1"], df["Stalling_Change_G3vsG2"])

    # 과학적 표기법으로 데이터프레임 생성
    corr_df = pd.DataFrame({
        'Comparison': [f'{label} (Common genes)'],
        'Correlation_coefficient': [f"{corr:.6f}"],
        'p_value': [f"{pval:.2e}"]  # 과학적 표기법 적용
    })

    # CSV로 저장
    corr_df.to_csv(f"Stalling_Index_Correlation_{label.replace(', ', '_').replace(' ', '')}.csv", index=False)

    # 결과 출력
    print(f"Correlation of Stalling Index Changes ({label}):")
    print(corr_df)

    # 개별 Stalling Index 통계 분석 및 저장
    stats_df = df[["Stalling_Index_G1", "Stalling_Index_G2", "Stalling_Index_G3"]].describe().transpose()
    stats_df.to_csv(f"Stalling_Index_Summary_{label.replace(', ', '_').replace(' ', '')}.csv")

    print(f"\nSummary Statistics - {label}:")
    print(stats_df)

###############################################################################################
#Correlation of Stalling Index Changes (G2 up, G3 down):
#                      Comparison Correlation_coefficient    p_value
#0  G2 up, G3 down (Common genes)               -0.891114  3.10e-293

#Summary Statistics - G2 up, G3 down:
#                   count      mean       std  min       25%       50%       75%       max
#Stalling_Index_G1  850.0  0.449632  0.231373  0.0  0.276629  0.455507  0.607540  0.997655
#Stalling_Index_G2  850.0  0.545290  0.287985  0.0  0.323136  0.596020  0.772358  1.000000
#Stalling_Index_G3  850.0  0.434900  0.223695  0.0  0.271400  0.438846  0.595837  0.996275
#Correlation of Stalling Index Changes (G2 down, G3 up):
#                      Comparison Correlation_coefficient   p_value
#0  G2 down, G3 up (Common genes)               -0.972983  0.00e+00

#Summary Statistics - G2 down, G3 up:
#                   count      mean       std  min       25%       50%       75%       max
#Stalling_Index_G1  835.0  0.362008  0.172915  0.0  0.224434  0.360682  0.484260  0.908623
#Stalling_Index_G2  835.0  0.387377  0.233868  0.0  0.198949  0.384703  0.558700  1.000000
#Stalling_Index_G3  835.0  0.367769  0.169377  0.0  0.231559  0.369554  0.486852  0.896309