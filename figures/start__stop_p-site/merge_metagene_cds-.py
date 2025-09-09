import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns

# ✅ 파일 경로 설정
base_path = "/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity/metagene_output/"

# ✅ 파일 자동 검색 및 그룹화 (Start와 Stop 구분)
group_files = {
    "START": {
        "CONTROL": sorted(glob.glob(base_path + "G1_*_start_codon_metagene.csv")),
        "HBSS": sorted(glob.glob(base_path + "G2_*_start_codon_metagene.csv")),
        "HBSS+HN": sorted(glob.glob(base_path + "G3_*_start_codon_metagene.csv"))
    },
    "STOP": {
        "CONTROL": sorted(glob.glob(base_path + "G1_*_stop_codon_metagene.csv")),
        "HBSS": sorted(glob.glob(base_path + "G2_*_stop_codon_metagene.csv")),
        "HBSS+HN": sorted(glob.glob(base_path + "G3_*_stop_codon_metagene.csv"))
    }
}

# ✅ 함수: 그룹별 평균 계산
def load_and_average(file_list):
    df_list = [pd.read_csv(f) for f in file_list]
    merged_df = pd.concat(df_list)
    
    # ✅ 숫자형 데이터만 선택 (문자열 제외)
    numeric_cols = merged_df.select_dtypes(include=[np.number]).columns
    
    # ✅ 그룹화 후 평균 계산 (rel_pos 기준)
    avg_df = merged_df.groupby("rel_pos", as_index=False)[numeric_cols].mean()
    
    return avg_df

# ✅ Start Codon과 Stop Codon 각각 처리
for codon_type in ["START", "STOP"]:
    grouped_dfs = []
    for group, files in group_files[codon_type].items():
        avg_df = load_and_average(files)
        avg_df["Group"] = group  # 그룹 라벨 추가
        grouped_dfs.append(avg_df)

    # ✅ 모든 그룹 병합
    grouped_df = pd.concat(grouped_dfs)

    # ✅ 중복된 `rel_pos`가 그룹 내에서 발생했는지 확인
    print(f"🔍 중복된 (rel_pos, Group) 개수 ({codon_type}):", grouped_df.duplicated(subset=["rel_pos", "Group"]).sum())

    # ✅ 중복된 경우 제거
    grouped_df = grouped_df.drop_duplicates(subset=["rel_pos", "Group"])

    # ✅ NaN 값이 포함된 경우 0으로 대체
    grouped_df["RPM"] = grouped_df["RPM"].fillna(0)

    # ✅ `rel_pos` 기준으로 정렬
    grouped_df = grouped_df.sort_values(by=["rel_pos", "Group"])

    # ✅ CSV 저장
    output_file = base_path + f"G_GROUPED_{codon_type}_METAGENE.csv"
    grouped_df.to_csv(output_file, index=False)

    print(f"✅ {codon_type} 그룹 평균 계산 완료! 저장: {output_file}")

    # ✅ 색상 설정 (RGB 값을 0~1로 변환)
    colors = {
        "CONTROL": (46/255, 139/255, 87/255),  # 초록색
        "HBSS": (255/255, 36/255, 0/255),    # 주황색
        "HBSS+HN": (65/255, 105/255, 255/255)  # 파란색
    }

    # ✅ 선 두께 차별화
    linewidths = {
        "CONTROL": 6.0,  # 가장 굵게
        "HBSS": 2.0,  # 중간
        "HBSS+HN": 6.0  # 가장 얇게
    }

    # ✅ x, y 축 범위 설정
    min_x, max_x = grouped_df["rel_pos"].min(), grouped_df["rel_pos"].max()
    max_y = grouped_df["RPM"].max()

    # ✅ 시각화
    plt.figure(figsize=(8, 4), dpi=600)
    ax = plt.gca()

    # ✅ 그룹별로 개별적으로 선 그리기
    for group in ["CONTROL", "HBSS", "HBSS+HN"]:
        subset = grouped_df[grouped_df["Group"] == group]
        sns.lineplot(x=subset["rel_pos"], y=subset["RPM"], label=group,
                     color=colors[group], linewidth=linewidths[group], linestyle="-")  # 모두 실선

    # ✅ x축 간격 설정
    if codon_type == "START":
        plt.xticks(np.arange(min_x, max_x+1, 100), fontsize=16)  # 100 nt 간격
    else:  # STOP Codon
        plt.xticks(np.arange(min_x, max_x+1, 50), fontsize=16)  # 50 nt 간격

    # ✅ y축 간격 설정 (20,000 간격)
    plt.yticks(np.arange(0, max_y+1, 20000), fontsize=13)

    # ✅ 3자리 단위마다 콤마 표시
    ax.set_yticklabels([f"{int(x):,}" for x in ax.get_yticks()])

    # ✅ x, y축만 유지하고 내부 격자선 삭제
    plt.grid(False)  # 내부 격자 제거

    # ✅ 그래프 설정 (기준선 제거)
    plt.xlabel(f"Distance from {codon_type} codon", fontsize=20, fontweight="bold", labelpad=10)
    plt.ylabel("Ribosome Occupancy", fontsize=20, fontweight="bold", labelpad=10)


    # ✅ x, y축의 두께 2배 증가
    ax.spines['bottom'].set_linewidth(2.5)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['top'].set_visible(False)  # 위쪽 선 제거
    ax.spines['right'].set_visible(False)  # 오른쪽 선 제거


    # ✅ 범례를 완벽히 중앙에 정렬 + 그룹3 바깥으로 안 나가게 조정
    plt.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.17), 
               fontsize=20, ncol=3, handletextpad=0.5, columnspacing=0.5)


    # ✅ 그래프 저장
    plot_output = base_path + f"G_GROUPED_{codon_type}_METAGENE_CLEAN.png"
    plt.savefig(plot_output, dpi=600, bbox_inches="tight")  # 🔥 600 DPI 저장
    plt.show()

    print(f"✅ {codon_type} 그룹별 Metagene 그래프 저장 완료: {plot_output}")
