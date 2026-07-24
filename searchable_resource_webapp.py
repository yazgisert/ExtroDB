"""
searchable_resource_webapp.py

FastAPI backend for the Extremophile Proteome searchable resource.

Endpoints
---------
GET /api/index       — full master index as JSON
GET /api/summary     — branch summary as JSON
GET /api/download/{proteome_id}
                     — stream the protein FASTA file for one proteome
GET /               — serves index.html (development only; in production nginx does this)

Run locally
-----------
    pip install fastapi uvicorn
    uvicorn searchable_resource_webapp:app --host 127.0.0.1 --port 8123 --reload
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

# ── Paths ────────────────────────────────────────────────────────────────────

BASE_DIR     = Path(__file__).resolve().parent
APP_DIR      = BASE_DIR / "searchable_resource_app"
INDEX_PATH   = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"
SUMMARY_PATH = BASE_DIR / "searchable_resource" / "extremophile_master_summary.tsv"

# ── App ───────────────────────────────────────────────────────────────────────

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

# ── Helpers ───────────────────────────────────────────────────────────────────

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def iso_timestamp(path: Path) -> str:
    modified = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)
    return modified.isoformat().replace("+00:00", "Z")


def fasta_streamer(path: Path) -> Iterator[bytes]:
    """Yield the FASTA file in 64 KB chunks so large files never fully load into memory."""
    with path.open("rb") as fh:
        while chunk := fh.read(65536):
            yield chunk


# ── API routes ────────────────────────────────────────────────────────────────

@app.get("/api/index")
def get_index():
    """Return the full master index as JSON."""
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
    """Return the branch summary as JSON."""
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
    """
    Stream the protein FASTA file for the requested proteome.

    The proteome_id is looked up in the master index to find the local FASTA
    path. The file is streamed in chunks so large proteomes (tens of thousands
    of proteins) never have to be fully loaded into memory.
    """
    # Find the row for this proteome_id
    rows = read_tsv(INDEX_PATH)
    row  = next((r for r in rows if r["proteome_id"] == proteome_id), None)

    if row is None:
        raise HTTPException(status_code=404, detail=f"Proteome '{proteome_id}' not found.")

    local_path_str = row.get("local_fasta_path", "").strip()
    if not local_path_str:
        raise HTTPException(status_code=404, detail="No local FASTA file recorded for this proteome.")

    fasta_path = BASE_DIR / local_path_str
    if not fasta_path.exists():
        raise HTTPException(status_code=404, detail="FASTA file not found on disk.")

    filename = fasta_path.name
    return StreamingResponse(
        fasta_streamer(fasta_path),
        media_type="application/octet-stream",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


# ── Static frontend (development only) ───────────────────────────────────────
# In production Docker deployment, nginx serves the frontend directly.
# This mount lets you run the full app locally with a single uvicorn command.

if APP_DIR.exists():
    app.mount("/", StaticFiles(directory=str(APP_DIR), html=True), name="frontend")
