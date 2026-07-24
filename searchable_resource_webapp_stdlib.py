from __future__ import annotations

import argparse
import csv
import json
import mimetypes
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse


BASE_DIR = Path(__file__).resolve().parent
APP_DIR = BASE_DIR / "searchable_resource_app"
INDEX_PATH = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"
SUMMARY_PATH = BASE_DIR / "searchable_resource" / "extremophile_master_summary.tsv"
TSV_DELIMITER = "\t"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=TSV_DELIMITER))


def iso_timestamp(path: Path) -> str:
    modified = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)
    return modified.isoformat().replace("+00:00", "Z")


def safe_join(base: Path, relative_path: str) -> Path | None:
    candidate = (base / relative_path).resolve()
    try:
        candidate.relative_to(base.resolve())
    except ValueError:
        return None
    return candidate


def dataset_payload() -> dict[str, object]:
    rows = read_tsv(INDEX_PATH)
    total_protein_count = sum(int(row["protein_count"]) for row in rows if row.get("protein_count"))
    return {
        "meta": {
            "source_file": str(INDEX_PATH.relative_to(BASE_DIR)),
            "updated_at": iso_timestamp(INDEX_PATH),
            "row_count": len(rows),
            "protein_count_sum": total_protein_count,
        },
        "rows": rows,
    }


def summary_payload() -> dict[str, object]:
    rows = read_tsv(SUMMARY_PATH)
    return {
        "meta": {
            "source_file": str(SUMMARY_PATH.relative_to(BASE_DIR)),
            "updated_at": iso_timestamp(SUMMARY_PATH),
            "row_count": len(rows),
        },
        "rows": rows,
    }


class SearchableResourceHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(APP_DIR), **kwargs)

    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        route = parsed.path

        if route == "/api/index":
            self.send_json(dataset_payload())
            return

        if route == "/api/summary":
            self.send_json(summary_payload())
            return

        if route.startswith("/files/"):
            relative_path = unquote(route[len("/files/") :])
            target = safe_join(BASE_DIR, relative_path)
            if target is None or not target.exists() or not target.is_file():
                self.send_error(HTTPStatus.NOT_FOUND, "File not found.")
                return
            self.send_local_file(target)
            return

        if route == "/":
            self.path = "/index.html"

        super().do_GET()

    def send_json(self, payload: dict[str, object]) -> None:
        body = json.dumps(payload, ensure_ascii=True).encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_local_file(self, path: Path) -> None:
        content_type = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        body = path.read_bytes()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Serve a small local web app for the searchable extremophile resource."
    )
    parser.add_argument("--host", default="127.0.0.1", help="Host interface to bind.")
    parser.add_argument("--port", type=int, default=8123, help="Port to serve on.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    server = ThreadingHTTPServer((args.host, args.port), SearchableResourceHandler)
    print(f"Serving searchable resource app at http://{args.host}:{args.port}")
    server.serve_forever()


if __name__ == "__main__":
    main()
