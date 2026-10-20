import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Combine BUSCO JSON results into one TSV table."
    )

    parser.add_argument(
        "--result",
        action="append",
        required=True,
        help="Assembly name and BUSCO JSON path: NAME=FILE"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV filename."
    )

    args = parser.parse_args()

    rows = []

    for result_argument in args.result:
        assembly_name, json_file = result_argument.split("=", 1)

        with open(json_file, "r") as busco_file:
            busco_data = json.load(busco_file)

        results = busco_data["results"]

        row = {
            "assembly": assembly_name,
            "complete": results["Complete"],
            "single_copy": results["Single copy"],
            "duplicated": results["Multi copy"],
            "fragmented": results["Fragmented"],
            "missing": results["Missing"],
            "total_buscos": results["n_markers"]
        }

        rows.append(row)

    results_table = pd.DataFrame(rows)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    results_table.to_csv(
        output_path,
        sep="\t",
        index=False
    )

    print(results_table.to_string(index=False))
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
