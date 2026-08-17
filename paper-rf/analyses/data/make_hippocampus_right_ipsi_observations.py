#!/usr/bin/env python3
"""Create a hippocampal observation table with raw right hemisphere as ipsilateral."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="paper-rf/data/hippocampus/observations.csv",
        help="Existing hippocampal observations table with i/c columns.",
    )
    parser.add_argument(
        "--network",
        default="paper-rf/data/hippocampus/network.csv",
        help="Network CSV used to determine model region labels.",
    )
    parser.add_argument(
        "--output",
        default="paper-rf/data/hippocampus/observations_right_ipsi.csv",
        help="Output observations table with right hemisphere mapped to i* labels.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    network_path = Path(args.network)
    output_path = Path(args.output)

    obs = pd.read_csv(input_path)
    network = pd.read_csv(network_path, nrows=0)
    labels = list(network.columns[1:])
    label_set = set(labels)

    out = obs.copy()
    swapped = []
    unpaired = []

    for label in labels:
        if not label.startswith("i"):
            continue
        mate = "c" + label[1:]
        if mate not in label_set:
            unpaired.append(label)
            continue
        if label not in out.columns or mate not in out.columns:
            raise ValueError(f"Missing matched observation columns: {label}, {mate}")
        left_values = obs[label].copy()
        right_values = obs[mate].copy()
        out[label] = right_values
        out[mate] = left_values
        swapped.append((label, mate))

    for label in labels:
        if label.startswith("c"):
            mate = "i" + label[1:]
            if mate not in label_set:
                unpaired.append(label)

    for label in sorted(set(unpaired)):
        if label in out.columns:
            out[label] = pd.NA

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    print(f"Wrote {output_path}")
    print(f"Swapped {len(swapped)} i/c label pairs.")
    if unpaired:
        print("Marked unpaired labels as missing: " + ", ".join(sorted(set(unpaired))))


if __name__ == "__main__":
    main()
