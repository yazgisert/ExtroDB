from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from common import BASE_DIR, TSV_DELIMITER, read_tsv

BASE = BASE_DIR
DEFAULT_INDEX = BASE / "searchable_resource" / "extremophile_master_index.tsv"


def parse_filter(values: List[str]) -> List[Tuple[str, str]]:
    # Filters use a lightweight column=value syntax so you can query the table
    # from the terminal without opening it in a spreadsheet first.
    filters = []
    for value in values:
        if "=" not in value:
            raise ValueError(f"Filter must use column=value syntax: {value}")
        column, expected = value.split("=", 1)
        filters.append((column, expected))
    return filters


def apply_filters(
    rows: List[Dict[str, str]],
    exact_filters: List[Tuple[str, str]],
    contains_filters: List[Tuple[str, str]],
) -> List[Dict[str, str]]:
    # Exact and substring filters can be combined in one command.
    filtered = rows

    for column, expected in exact_filters:
        filtered = [row for row in filtered if row.get(column, "") == expected]

    for column, expected in contains_filters:
        expected_lower = expected.lower()
        filtered = [
            row for row in filtered if expected_lower in row.get(column, "").lower()
        ]

    return filtered


def write_rows(rows: List[Dict[str, str]], columns: List[str]) -> None:
    # Writing TSV to stdout keeps this script easy to pipe into files or other
    # shell commands later.
    sys.stdout.write(TSV_DELIMITER.join(columns) + "\n")
    for row in rows:
        sys.stdout.write(TSV_DELIMITER.join(row.get(column, "") for column in columns) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter the searchable extremophile master index."
    )
    parser.add_argument(
        "--index",
        default=str(DEFAULT_INDEX),
        help="Path to the master index TSV.",
    )
    parser.add_argument(
        "--where",
        action="append",
        default=[],
        help="Exact filter using column=value syntax. Repeat as needed.",
    )
    parser.add_argument(
        "--contains",
        action="append",
        default=[],
        help="Substring filter using column=value syntax. Repeat as needed.",
    )
    parser.add_argument(
        "--columns",
        default="",
        help="Comma-separated output columns. Default is all columns.",
    )
    parser.add_argument(
        "--count",
        action="store_true",
        help="Print only the number of matched rows.",
    )
    parser.add_argument(
        "--list-columns",
        action="store_true",
        help="Print all available column names and exit.",
    )
    args = parser.parse_args()

    rows = read_tsv(Path(args.index))
    if not rows:
        return

    if args.list_columns:
        for column in rows[0].keys():
            print(column)
        return

    exact_filters = parse_filter(args.where)
    contains_filters = parse_filter(args.contains)
    filtered = apply_filters(rows, exact_filters, contains_filters)

    if args.count:
        print(len(filtered))
        return

    if args.columns:
        columns = [column.strip() for column in args.columns.split(",") if column.strip()]
    else:
        columns = list(rows[0].keys())

    write_rows(filtered, columns)


if __name__ == "__main__":
    main()
