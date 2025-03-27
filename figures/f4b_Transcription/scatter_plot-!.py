# Re-plotting Transcription scatter plot with legend including counts

plt.figure(figsize=(7, 7), dpi=300)

# Plot non-significant points (gray)
plt.scatter(
    non_significant_rna_points["G1_G2_fc"],
    non_significant_rna_points["G2_G3_fc"],
    s=2,
    c="dimgray",
    alpha=0.5,
    label=f"Non-Significant (n={num_non_significant_rna})"
)

# Plot significant points (red)
plt.scatter(
    significant_rna_points["G1_G2_fc"],
    significant_rna_points["G2_G3_fc"],
    s=3,
    c="firebrick",
    alpha=0.7,
    label=f"Significant (n={num_significant_rna})"
)

# Adding dash lines only at -1 and 1
plt.axhline(1, color="black", linewidth=0.7, linestyle="--")
plt.axhline(-1, color="black", linewidth=0.7, linestyle="--")
plt.axvline(1, color="black", linewidth=0.7, linestyle="--")
plt.axvline(-1, color="black", linewidth=0.7, linestyle="--")

# Adding a box around the plot
plt.gca().spines["top"].set_visible(True)
plt.gca().spines["right"].set_visible(True)
plt.gca().spines["left"].set_visible(True)
plt.gca().spines["bottom"].set_visible(True)
plt.gca().spines["top"].set_linewidth(1)
plt.gca().spines["right"].set_linewidth(1)
plt.gca().spines["left"].set_linewidth(1)
plt.gca().spines["bottom"].set_linewidth(1)

# Enlarging title, x-axis, y-axis labels, and ticks
plt.title("Transcription", fontsize=32, ha="center", pad=20)  # Title enlarged to 32pt
plt.xlabel("Nutrient deprivation", fontsize=30, labelpad=10)  # Axis labels closer to axes
plt.ylabel("Humanin treatment", fontsize=30, labelpad=10)  # Axis labels closer to axes
plt.xticks(ticks=[-8, -4, 0, 4, 8], fontsize=24)  # Adjusted tick intervals
plt.yticks(ticks=[-8, -4, 0, 4, 8], fontsize=24)  # Adjusted tick intervals
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.gca().set_aspect("equal", adjustable="box")

# Adjust legend font size and positioning
plt.legend(fontsize=24, loc="upper right", frameon=False)
plt.tight_layout()

# Save the updated RNA-seq plot (Transcription) with legend including counts
tiff_path_rna_boxed_with_legend = "/mnt/data/scatter_rna_boxed_with_legend.tiff"
plt.savefig(tiff_path_rna_boxed_with_legend, format="tiff", dpi=600, bbox_inches="tight")
plt.show()

# Return file path
tiff_path_rna_boxed_with_legend
