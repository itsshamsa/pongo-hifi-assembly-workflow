import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


summary = pd.read_csv(
    "results/assembly_summary.tsv",
    sep="\t"
)

summary["total_gb"] = summary["total_bp"] / 1_000_000_000
summary["n50_mb"] = summary["n50_bp"] / 1_000_000


figure, axes = plt.subplots(
    1,
    3,
    figsize=(14, 4)
)


axes[0].bar(
    summary["assembly"],
    summary["total_gb"]
)

axes[0].set_title("Assembly size")
axes[0].set_ylabel("Total length (Gb)")


axes[1].bar(
    summary["assembly"],
    summary["contigs"]
)

axes[1].set_title("Fragmentation")
axes[1].set_ylabel("Number of contigs")


axes[2].bar(
    summary["assembly"],
    summary["n50_mb"]
)

axes[2].set_title("Contiguity")
axes[2].set_ylabel("N50 (Mb)")


figure.suptitle("Pongo abelii hifiasm assembly comparison")

plt.tight_layout()

plt.savefig(
    "results/assembly_summary.png",
    dpi=300
)
