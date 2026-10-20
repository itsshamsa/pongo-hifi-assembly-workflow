import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Create a stacked BUSCO comparison plot."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="BUSCO summary table in TSV format."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output PNG filename."
    )

    parser.add_argument(
        "--title",
        default="BUSCO completeness",
        help="Title displayed above the figure."
    )

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    categories = [
        "single_copy",
        "duplicated",
        "fragmented",
        "missing"
    ]

    colors = {
        "single_copy": "#56B4E9",
        "duplicated": "#0072B2",
        "fragmented": "#E69F00",
        "missing": "#D55E00"
    }

    bottom = [0] * len(df)

    fig, ax = plt.subplots(figsize=(7, 5))

    for category in categories:
        ax.bar(
            df["assembly"],
            df[category],
            bottom=bottom,
            label=category.replace("_", " ").title(),
            color=colors[category]
        )

        bottom = [
            previous + value
            for previous, value in zip(bottom, df[category])
        ]

    ax.set_ylabel("BUSCO genes (%)")
    ax.set_xlabel("Assembly")
    ax.set_title(args.title)
    ax.set_ylim(0, 100)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.tight_layout()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(df.to_string(index=False))
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
