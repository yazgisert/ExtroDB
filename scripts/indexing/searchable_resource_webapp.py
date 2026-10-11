# FastAPI backend of ExtroDB: master index and summary as JSON, single and bulk (zip) FASTA downloads, and the static
# frontend. Runs in Docker (see README).

from __future__ import annotations

import csv
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from fastapi import BackgroundTasks, FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# In Docker the script sits next to searchable_resource/; locally it lives in scripts/indexing/
_SCRIPT_DIR  = Path(__file__).resolve().parent
BASE_DIR     = _SCRIPT_DIR if (_SCRIPT_DIR / "searchable_resource").exists() else _SCRIPT_DIR.parents[1]
APP_DIR      = BASE_DIR / "searchable_resource_app"
INDEX_PATH   = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"
SUMMARY_PATH = BASE_DIR / "searchable_resource" / "extremophile_master_summary.tsv"

# FASTA_DIR is where the *_raw_protein_fastas/ folders live.
# Defaults to BASE_DIR (local dev with host-mounted folders).
# In the Docker Zenodo setup, set FASTA_DIR=/app/fasta_data via environment.
FASTA_DIR = Path(os.environ.get("FASTA_DIR", str(BASE_DIR)))

app = FastAPI(
    title="Extremophile Proteome Resource",
    description="Browse and download curated extremophile proteomes.",
    version="2.0.0",
)

# CORS — allows the frontend (served by nginx on port 80) to call this backend
# (served on port 8000) without the browser blocking the request.
# In production, replace "*" with your actual domain.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def iso_timestamp(path: Path) -> str:
    modified = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)
    return modified.isoformat().replace("+00:00", "Z")


@app.get("/api/index")
def get_index():
    rows = read_tsv(INDEX_PATH)
    total_proteins = sum(
        int(row["protein_count"])
        for row in rows
        if row.get("protein_count", "").isdigit()
    )
    return JSONResponse({
        "meta": {
            "source_file": str(INDEX_PATH.relative_to(BASE_DIR)),
            "updated_at":  iso_timestamp(INDEX_PATH),
            "row_count":   len(rows),
            "protein_count_sum": total_proteins,
        },
        "rows": rows,
    })


@app.get("/api/summary")
def get_summary():
    rows = read_tsv(SUMMARY_PATH)
    return JSONResponse({
        "meta": {
            "source_file": str(SUMMARY_PATH.relative_to(BASE_DIR)),
            "updated_at":  iso_timestamp(SUMMARY_PATH),
            "row_count":   len(rows),
        },
        "rows": rows,
    })


@app.get("/api/download/{proteome_id}")
def download_fasta(proteome_id: str):
    rows = read_tsv(INDEX_PATH)
    row  = next((r for r in rows if r["proteome_id"] == proteome_id), None)

    if row is None:
        raise HTTPException(status_code=404, detail=f"Proteome '{proteome_id}' not found.")

    # Older index versions had a plural, pipe-separated column; take the first path.
    local_path_str = (row.get("local_fasta_path") or row.get("local_fasta_paths") or "").strip()
    local_path_str = local_path_str.split("|")[0].strip()

    if not local_path_str:
        raise HTTPException(status_code=404, detail="No local FASTA file recorded for this proteome.")

    fasta_path = FASTA_DIR / local_path_str
    if not fasta_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"FASTA file not found on disk: {local_path_str}",
        )

    return FileResponse(
        path=str(fasta_path),
        media_type="application/octet-stream",
        filename=fasta_path.name,
    )


class BulkDownloadRequest(BaseModel):
    ids: list[str]


@app.post("/api/download-bulk")
async def download_bulk(request: BulkDownloadRequest, background_tasks: BackgroundTasks):
    # IDs without a FASTA on disk are skipped.
    rows = read_tsv(INDEX_PATH)
    row_map = {r["proteome_id"]: r for r in rows if r.get("proteome_id")}

    files_to_zip: list[Path] = []
    for pid in request.ids:
        row = row_map.get(pid)
        if not row:
            continue
        local_path_str = (row.get("local_fasta_path") or row.get("local_fasta_paths") or "").strip()
        local_path_str = local_path_str.split("|")[0].strip()
        if not local_path_str:
            continue
        fasta_path = FASTA_DIR / local_path_str
        if fasta_path.exists():
            files_to_zip.append(fasta_path)

    if not files_to_zip:
        raise HTTPException(status_code=404, detail="No FASTA files found for the requested IDs.")

    # Write zip to a temp file so we can stream it without buffering everything in RAM
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".zip")
    os.close(tmp_fd)

    with ZipFile(tmp_path, "w", ZIP_DEFLATED) as zf:
        for path in files_to_zip:
            zf.write(path, path.name)

    background_tasks.add_task(os.unlink, tmp_path)

    return FileResponse(
        path=tmp_path,
        media_type="application/zip",
        filename="extremophile_fastas.zip",
    )


if APP_DIR.exists():
    app.mount("/", StaticFiles(directory=str(APP_DIR), html=True), name="frontend")
