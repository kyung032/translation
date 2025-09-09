import os
import time
import matplotlib.pyplot as plt
import pandas as pd

# ✅ 논문용 그래프 설정 (폰트 크기 & 선 굵기 조정)
plt.rcParams.update({
    "font.size": 18,       # 폰트 크기 증가
    "axes.labelweight": "bold",
    "axes.titlesize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16
})

# ✅ 데이터 파일 위치
input_dir = "/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity_TE"
output_dir = "/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity_TE/metagene_output"
os.makedirs(output_dir, exist_ok=True)

# ✅ 그룹별 파일 목록 (직접 지정)
group_files = {
    "G1": ["G1_1_collapsed_transcripts_17368.csv", "G1_2_collapsed_transcripts_17368.csv"],
    "G2": ["G2_1_collapsed_transcripts_14769.csv", "G2_2_collapsed_transcripts_15039.csv"],
    "G3": ["G3_1_collapsed_transcripts_17473.csv", "G3_2_collapsed_transcripts_17610.csv"]
}

# ✅ 포함할 유전자 목록 (5개)
ensembl_gene_ids = {
    "CTDSP1": "ENSG00000144579",
    "ENTPD4": "ENSG00000197217",
    "EPHB2": "ENSG00000133216",
    "PPIA": "ENSG00000196262",
    "PPIL2": "ENSG00000112695"
}

# ✅ X축 범위 고정 (-50 ~ +50)
x_min, x_max = -50, 50

# ✅ 그룹별 색상 설정
sample_colors = {
    "G1": (46/255, 139/255, 87/255),   # Green
    "G2": (255/255, 36/255, 0/255),    # Red
    "G3": (65/255, 105/255, 225/255)   # Blue
}

# ✅ 데이터 타입 지정 (속도 향상)
dtypes = {
    "gene_id": "category",
    "Feature": "category",
    "rel_pos": "int32",
    "psite_center": "int8",
    "normalized_density": "float32"
}

# ✅ 유전자별 데이터를 로드 및 처리 함수 (속도 개선)
def load_and_process_gene_data(gene_name, gene_id):
    """특정 유전자에 대한 모든 그룹 데이터를 로드하고 CDS만 필터링하여 정리"""
    combined_data = []

    # ✅ 그룹별 데이터 로드 및 합산
    for group, file_list in group_files.items():
        all_files = []

        for file in file_list:
            file_path = os.path.join(input_dir, file)
            if not os.path.exists(file_path):
                print(f"⚠️ 파일 없음: {file_path}")
                continue

            # ✅ CSV 로드 (속도 최적화)
            df = pd.read_csv(file_path, sep=None, engine='python', skipinitialspace=True, header=0, dtype=dtypes)

            # ✅ `rel_pos`를 int 변환
            df["rel_pos"] = df["rel_pos"].astype("int32")

            # ✅ 해당 유전자의 CDS 필터링
            df_filtered = df[(df["gene_id"] == gene_id) & (df["Feature"] == "CDS")]

            # ✅ 그룹 컬럼 추가
            df_filtered = df_filtered[["rel_pos", "normalized_density"]]
            df_filtered.rename(columns={"normalized_density": group}, inplace=True)

            all_files.append(df_filtered)

        # ✅ 그룹 내 여러 파일 합산
        if all_files:
            group_df = pd.concat(all_files).groupby("rel_pos").sum().reset_index()
            combined_data.append(group_df)

    # ✅ 모든 그룹 데이터를 병합
    if combined_data:
        final_df = combined_data[0]
        for df in combined_data[1:]:
            final_df = final_df.merge(df, on="rel_pos", how="outer")

        # ✅ 누락된 좌표를 0으로 채우기
        final_df.fillna(0, inplace=True)

        return final_df
    return pd.DataFrame()  # 데이터 없으면 빈 데이터 반환

# ✅ 논문용 플롯 함수
def plot_rpf_density(df, gene_name, filename):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=600)  # 해상도 증가

    # ✅ X축 맞추기 (-50 ~ +100)
    ax.set_xlim(x_min, x_max)


    # ✅ Y축 변환 (PPIA만 `10^-4` 적용)
    if gene_name == "PPIA":
        df.iloc[:, 1:] /= 1e-4  # 모든 그룹(G1, G2, G3) 값을 10^-4로 변환
        y_label = "Normalized Density (×10⁻⁴)"
    else:
        y_label = "Normalized Density"

    # ✅ 그룹별 플롯 (선 굵기 3)
    for group in ["G1", "G2", "G3"]:
        if group in df.columns:
            ax.plot(df["rel_pos"], df[group], label=f"{group} - CDS", color=sample_colors[group], linewidth=3)

    # ✅ 축 라벨 및 제목
    ax.set_xlabel("Relative Position", fontsize=18, fontweight="bold")
    ax.set_ylabel("Normalized Density", fontsize=18, fontweight="bold")
    ax.set_title(f"RPF Reads Density - {gene_name} (CDS)", fontsize=20, fontweight="bold")

    # ✅ 스타일 적용
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)

    # ✅ Legend를 플롯 바깥으로 배치
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=16, frameon=True, edgecolor="black")


    # ✅ TIFF 및 SVG/PDF 저장
    plt.savefig(os.path.join(output_dir, f"{filename}.tiff"), format="tiff", dpi=600, bbox_inches="tight")
    plt.savefig(os.path.join(output_dir, f"{filename}.svg"), format="svg", bbox_inches="tight")
    plt.savefig(os.path.join(output_dir, f"{filename}.pdf"), format="pdf", bbox_inches="tight")
    plt.close()

# ✅ 실행 속도 테스트
start_time = time.time()

# ✅ 모든 유전자에 대해 처리
for gene_name, gene_id in ensembl_gene_ids.items():  
    print(f"🔄 {gene_name} 처리 중...")
    df_gene = load_and_process_gene_data(gene_name, gene_id)
    
    if not df_gene.empty:
        # ✅ CSV 파일 저장
        csv_filename = os.path.join(output_dir, f"{gene_name}_CDS_RPF.csv")
        df_gene.to_csv(csv_filename, index=False)

        # ✅ 그래프 생성
        plot_rpf_density(df_gene, gene_name, f"{gene_name}_CDS_RPF")

end_time = time.time()
print(f"✅ 실행 완료! 총 소요 시간: {end_time - start_time:.2f}초")

# ✅ 5개 유전자에 대해 그래프 다시 생성
for gene_name in ensembl_gene_ids.keys():
    file_path = os.path.join(output_dir, f"{gene_name}_CDS_RPF.csv")
    if os.path.exists(file_path):
        df_gene = pd.read_csv(file_path)
        plot_rpf_density(df_gene, gene_name, f"{gene_name}_CDS_RPF")
        print(f"✅ {gene_name} 그래프 저장 완료!")

print("🎯 논문용 고해상도 그래프 생성 완료!")
