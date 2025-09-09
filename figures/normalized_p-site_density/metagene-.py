
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 📌 1️⃣ 각 그룹의 binned P-site density 파일 리스트
group_files = {
    "Control": ["G1_1_binned_p_site_density.csv", "G1_2_binned_p_site_density.csv"],
    "HBSS": ["G2_1_binned_p_site_density.csv", "G2_2_binned_p_site_density.csv"],
    "HBSS+HN": ["G3_1_binned_p_site_density.csv", "G3_2_binned_p_site_density.csv"]
}

# 📌 2️⃣ 그룹별 데이터 로드 및 병합
group_colors = {"Control": "#66CC66", "HBSS": "#E64A19", "HBSS+HN": "#0066FF"}
group_data = {}

for group, files in group_files.items():
    df_list = [pd.read_csv(file)[["Bin_center", "N"]].dropna() for file in files]
    merged_df = pd.concat(df_list)
    group_data[group] = merged_df

# 📌 3️⃣ 개선된 KDE 그래프 생성
fig, ax = plt.subplots(figsize=(8, 6), dpi=300)  # DPI 300 유지

for group, df in group_data.items():
    sns.kdeplot(
        x=df["Bin_center"], 
        weights=df["N"],  
        fill=True, 
        label=group, 
        alpha=0.3,  # 🔥 그룹별 투명도 설정 (기존 0.2 → 0.3)
        bw_adjust=0.2, 
        color=group_colors[group],
        linewidth=1.8  # 🔥 선 굵기 증가 (기존 대비 더 뚜렷하게)
    )

# 📌 4️⃣ Start/Stop Codon 강조 (검정 점선, 레이블 없음)
ax.axvline(x=0, color="black", linestyle="dashed", linewidth=1.5)
ax.axvline(x=1, color="black", linestyle="dashed", linewidth=1.5)

# 📌 5️⃣ Y축 간격을 0.5 단위로 조정
y_min, y_max = ax.get_ylim()
ax.set_yticks(np.arange(y_min, y_max + 0.5, 0.5))  

# 📌 6️⃣ 그래프 설정 (X축 기존 유지)
ax.set_xlabel("CDS position", fontsize=20, fontweight="bold")
ax.set_ylabel("Normalized P-site Density", fontsize=20, fontweight="bold")
ax.tick_params(axis="both", labelsize=16)

# 📌 7️⃣ 범례 조정 (기존 위치 유지, 좀 더 보기 좋게 조정)
ax.legend(
    fontsize=16,  # 🔥 그룹명 글씨 크기 조절 (기존 16 → 14)
    loc="upper left",  # 기존 위치 유지
    bbox_to_anchor=(0.1, 0.9),  # 🔥 범례를 그래프 내부 적절한 위치로 조정
    frameon=True,  # 🔲 범례 박스 테두리 추가
    facecolor="white",  # 배경 흰색 유지
    edgecolor="black",  # 🔥 범례 테두리 검정 추가
    labelspacing=0.6  # 🔥 범례 항목 간격 조정 (더 간결하게)
)

# 📌 8️⃣ 저장 (tight_layout 적용)
plt.tight_layout()
plt.savefig("Metagene_KDE_Science_Fixed.png", dpi=300)
plt.show()
