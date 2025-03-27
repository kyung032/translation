import os
import pandas as pd
import numpy as np

# ------------------------------------------------
# 0) 파일 경로 정의 (업로드된 파일 기반)
# ------------------------------------------------
# TE 유전자 목록 파일 (Up/Down)
te_up_g3vs2_file = "DESeq2_TE_Upregulated_G3vsG2.txt"
te_down_g3vs2_file = "DESeq2_TE_Downregulated_G3vsG2.txt"
te_up_g2vs1_file = "DESeq2_TE_Upregulated_G2vsG1.txt"
te_down_g2vs1_file = "DESeq2_TE_Downregulated_G2vsG1.txt"

# P-site 데이터 (각 그룹별 Replicate 포함)
g1_file1 = "/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity_TE/G1_1_collapsed_transcripts_17368.csv"
g1_file2 = "/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity_TE/G1_2_collapsed_transcripts_17368.csv"
g2_file1 = "/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity_TE/G2_1_collapsed_transcripts_14769.csv"
g2_file2 = "/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity_TE/G2_2_collapsed_transcripts_15039.csv"
g3_file1 = "/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity_TE/G3_1_collapsed_transcripts_17473.csv"
g3_file2 = "/Users/kyungkim/Documents/riboseq_2025/riboseq_analysis/polarity_TE/G3_2_collapsed_transcripts_17610.csv"


# 출력 폴더 생성
output_dir = "output_stalling_index"
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------
# 1) Stalling Index 계산 함수
# ------------------------------------------------
def compute_stalling_index(df, init_range=(-50, 50)):
    """특정 구간의 Ribo 밀도를 전체 밀도로 정규화하여 Stalling Index 계산"""
    df_sum = df.groupby(["gene_id", "rel_pos"], as_index=False)["normalized_density"].sum()
    
    # 특정 구간(init_range)만 골라서 합
    region_df = df_sum[
        (df_sum["rel_pos"] >= init_range[0]) & 
        (df_sum["rel_pos"] <= init_range[1])
    ]
    region_sum = region_df.groupby("gene_id", as_index=False)["normalized_density"].sum()
    
    # 전체 합
    overall_sum = df_sum.groupby("gene_id", as_index=False)["normalized_density"].sum()

    # Stalling Index 계산 (NaN 값 처리)
    merged = overall_sum.merge(region_sum, on="gene_id", suffixes=("_total", "_region"), how="left")
    merged["stalling_index"] = merged["normalized_density_region"] / merged["normalized_density_total"]
    merged = merged[["gene_id", "stalling_index"]].fillna(0)
    
    return merged

# ------------------------------------------------
# 2) 각 Replicate 불러오기
# ------------------------------------------------
ribo_g1_1 = pd.read_csv(g1_file1)
ribo_g1_2 = pd.read_csv(g1_file2)
ribo_g2_1 = pd.read_csv(g2_file1)
ribo_g2_2 = pd.read_csv(g2_file2)
ribo_g3_1 = pd.read_csv(g3_file1)
ribo_g3_2 = pd.read_csv(g3_file2)

# ------------------------------------------------
# 3) Stalling Index 분석 및 저장 함수
# ------------------------------------------------
def analyze_stalling_index(te_file, scenario_name, group):
    """
    특정 TE 유전자 목록에 대한 각 Replicate별 Stalling Index 계산 및 저장
    """
    # 3-1) 유전자 목록 불러오기
    te_genes = pd.read_csv(te_file, sep="\t")["GeneID"].tolist()

    # 3-2) 그룹별 Replicate 선택
    if group == "G2vsG1":
        r1_1, r1_2 = ribo_g1_1, ribo_g1_2
        r2_1, r2_2 = ribo_g2_1, ribo_g2_2
        col_names = ["G1_1_SI", "G1_2_SI", "G2_1_SI", "G2_2_SI"]
    elif group == "G3vsG2":
        r1_1, r1_2 = ribo_g2_1, ribo_g2_2
        r2_1, r2_2 = ribo_g3_1, ribo_g3_2
        col_names = ["G2_1_SI", "G2_2_SI", "G3_1_SI", "G3_2_SI"]
    else:
        raise ValueError("group은 'G2vsG1' 또는 'G3vsG2'만 가능합니다.")

    # 3-3) 각 Replicate에서 해당 유전자 필터링 후 Stalling Index 계산
    si_r1_1 = compute_stalling_index(r1_1[r1_1["gene_id"].isin(te_genes)])
    si_r1_2 = compute_stalling_index(r1_2[r1_2["gene_id"].isin(te_genes)])
    si_r2_1 = compute_stalling_index(r2_1[r2_1["gene_id"].isin(te_genes)])
    si_r2_2 = compute_stalling_index(r2_2[r2_2["gene_id"].isin(te_genes)])

    # 3-4) 유전자 리스트 통합 (길이 맞추기)
    all_genes = set(si_r1_1["gene_id"]) | set(si_r1_2["gene_id"]) | set(si_r2_1["gene_id"]) | set(si_r2_2["gene_id"])
    all_genes = sorted(list(all_genes))

    # 3-5) 모든 데이터프레임 병합 (NaN 값은 0으로 대체)
    df_stalling = pd.DataFrame({"gene_id": all_genes})
    df_stalling = df_stalling.merge(si_r1_1, on="gene_id", how="left").rename(columns={"stalling_index": col_names[0]})
    df_stalling = df_stalling.merge(si_r1_2, on="gene_id", how="left").rename(columns={"stalling_index": col_names[1]})
    df_stalling = df_stalling.merge(si_r2_1, on="gene_id", how="left").rename(columns={"stalling_index": col_names[2]})
    df_stalling = df_stalling.merge(si_r2_2, on="gene_id", how="left").rename(columns={"stalling_index": col_names[3]})
    
    # NaN 값은 0으로 처리
    df_stalling.fillna(0, inplace=True)

    # 3-6) CSV 저장
    output_file = os.path.join(output_dir, f"Stalling_Index_{scenario_name}.csv")
    df_stalling.to_csv(output_file, index=False)
    print(f"✅ 저장 완료: {output_file}")

# ------------------------------------------------
# 4) G3vsG2 Up, G3vsG2 Down, G2vsG1 Up, G2vsG1 Down 분석 실행
# ------------------------------------------------
analyze_stalling_index(te_up_g3vs2_file, "G3vsG2_Up", "G3vsG2")
analyze_stalling_index(te_down_g3vs2_file, "G3vsG2_Down", "G3vsG2")
analyze_stalling_index(te_up_g2vs1_file, "G2vsG1_Up", "G2vsG1")
analyze_stalling_index(te_down_g2vs1_file, "G2vsG1_Down", "G2vsG1")
