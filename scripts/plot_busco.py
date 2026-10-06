import pandas as pd
import matplotlib.pyplot as plt


# Read the tab-separated BUSCO results
df = pd.read_csv("results/busco_summary.tsv", sep="\t")

# Categories that together represent 100% of the BUSCO dataset
categories = ["single_copy", "duplicated", "fragmented", "missing"]

colors = {
    "single_copy": "#56B4E9",
    "duplicated": "#0072B2",
    "fragmented": "#E69F00",
    "missing": "#D55E00",
}

# This records where the next section of each stacked bar begins
bottom = [0] * len(df)

for category in categories:
    plt.bar(
        df["assembly"],
        df[category],
        bottom=bottom,
        label=category.replace("_", " ").title(),
        color=colors[category],
    )

    bottom = [
        previous + value
        for previous, value in zip(bottom, df[category])
    ]

plt.ylabel("BUSCO genes (%)")
plt.xlabel("Assembly")
plt.title("BUSCO completeness: Pongo abelii assemblies")
plt.ylim(0, 100)
plt.legend()
plt.tight_layout()

plt.savefig(
    "results/busco_comparison.png",
    dpi=300,
    bbox_inches="tight",
)

print(df)
print("Saved: results/busco_comparison.png")
