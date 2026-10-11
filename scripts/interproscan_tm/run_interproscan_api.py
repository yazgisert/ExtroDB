# Sends the sampled sequences to the EBI InterProScan REST API in batches of 25 (API limit 30), polls for the
# results and saves them as TSV. Resumable via a progress file.

import json
import time
import urllib.parse
import urllib.request
import urllib.error
from pathlib import Path

EMAIL      = "yazgi.sert@tum.de"
BATCH_SIZE = 25
POLL_INTERVAL = 30   # seconds between status checks
MAX_POLL      = 120  # max poll attempts per job (~1 hour)
MAX_RETRIES   = 5

BASE_URL   = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
OUT_DIR    = Path("results/interproscan")
PROGRESS_F = OUT_DIR / "progress.json"

INPUTS = {
    "prokaryote": OUT_DIR / "prokaryote_sample.faa",
    "eukaryote":  OUT_DIR / "eukaryote_sample.faa",
}
OUTPUTS = {
    "prokaryote": OUT_DIR / "prokaryote_results.tsv",
    "eukaryote":  OUT_DIR / "eukaryote_results.tsv",
}


def read_fasta(path):
    records = []
    current_header = None
    current_seq    = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    records.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq    = []
            else:
                current_seq.append(line)
    if current_header:
        records.append((current_header, "".join(current_seq)))
    return records


def format_fasta(records):
    return "\n".join(f">{h}\n{s}" for h, s in records)


def load_progress():
    if PROGRESS_F.exists():
        with open(PROGRESS_F) as f:
            return json.load(f)
    return {"completed_batches": {}, "results": {}}


def save_progress(progress):
    with open(PROGRESS_F, "w") as f:
        json.dump(progress, f, indent=2)


def api_post(endpoint, data):
    url  = f"{BASE_URL}/{endpoint}"
    body = urllib.parse.urlencode(data).encode()
    req  = urllib.request.Request(url, data=body, method="POST")
    req.add_header("Content-Type", "application/x-www-form-urlencoded")
    req.add_header("Accept", "text/plain")
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode().strip()


def api_get(endpoint, accept="text/plain"):
    url = f"{BASE_URL}/{endpoint}"
    req = urllib.request.Request(url)
    req.add_header("Accept", accept)
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode().strip()


def submit_job(fasta_str, retries=MAX_RETRIES):
    for attempt in range(1, retries + 1):
        try:
            job_id = api_post("run", {
                "email":     EMAIL,
                "sequence":  fasta_str,
                "stype":     "p",
                "database":  "pfam",
                "goterms":   "true",
                "pathways":  "false",
            })
            return job_id
        except Exception as e:
            wait = 2 ** attempt
            print(f"    Submit attempt {attempt}/{retries} failed: {e}. Retrying in {wait}s ...")
            time.sleep(wait)
    raise RuntimeError("Failed to submit job after all retries.")


def poll_job(job_id, retries=MAX_RETRIES):
    for poll_num in range(1, MAX_POLL + 1):
        try:
            status = api_get(f"status/{job_id}")
        except Exception as e:
            print(f"    Poll error (attempt {poll_num}): {e}")
            time.sleep(POLL_INTERVAL)
            continue

        if status == "FINISHED":
            return True
        elif status in ("FAILURE", "ERROR", "NOT_FOUND"):
            print(f"    Job {job_id} ended with status: {status}")
            return False
        else:
            if poll_num % 4 == 0:
                print(f"    [{job_id}] status={status} (poll {poll_num})")
            time.sleep(POLL_INTERVAL)

    print(f"    Job {job_id} timed out after {MAX_POLL} polls.")
    return False


def fetch_result(job_id, retries=MAX_RETRIES):
    for attempt in range(1, retries + 1):
        try:
            return api_get(f"result/{job_id}/tsv")
        except Exception as e:
            wait = 2 ** attempt
            print(f"    Fetch attempt {attempt}/{retries} failed: {e}. Retrying in {wait}s ...")
            time.sleep(wait)
    raise RuntimeError(f"Failed to fetch result for {job_id}.")


def run_group(group_name, faa_path, out_tsv, progress):
    print(f"\nProcessing {group_name}")

    records = read_fasta(faa_path)
    print(f"  {len(records):,} sequences loaded from {faa_path}")

    batches = [records[i:i+BATCH_SIZE] for i in range(0, len(records), BATCH_SIZE)]
    print(f"  {len(batches)} batches of up to {BATCH_SIZE}")

    completed = progress["completed_batches"].get(group_name, {})
    all_results = progress["results"].get(group_name, [])

    for batch_idx, batch in enumerate(batches):
        batch_key = str(batch_idx)
        if batch_key in completed:
            print(f"  Batch {batch_idx+1}/{len(batches)} already done — skipping.")
            continue

        fasta_str = format_fasta(batch)
        print(f"\n  Batch {batch_idx+1}/{len(batches)} ({len(batch)} sequences) ...")

        print(f"    Submitting ...")
        job_id = submit_job(fasta_str)
        print(f"    Job ID: {job_id}")

        success = poll_job(job_id)
        if not success:
            print(f"    Batch {batch_idx+1} failed — skipping.")
            continue

        tsv_result = fetch_result(job_id)
        all_results.append(tsv_result)

        completed[batch_key] = job_id
        progress["completed_batches"][group_name] = completed
        progress["results"][group_name] = all_results
        save_progress(progress)

        print(f"    Done. {len(tsv_result.splitlines())} result lines.")
        time.sleep(5)   # brief pause between submissions

    header = (
        "Sequence_ID\tMD5\tLength\tAnalysis\tSignature_Accession\t"
        "Signature_Description\tStart\tStop\tScore\tStatus\tDate\t"
        "InterPro_Accession\tInterPro_Description\tGO_Annotations\tPathways\n"
    )
    with open(out_tsv, "w") as f:
        f.write(header)
        for chunk in all_results:
            if chunk.strip():
                f.write(chunk.strip() + "\n")

    print(f"\n  Results written to {out_tsv}")


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    progress = load_progress()

    for group in ["prokaryote", "eukaryote"]:
        faa_path = INPUTS[group]
        out_tsv  = OUTPUTS[group]

        if not faa_path.exists():
            print(f"Missing input: {faa_path} — run sample_for_interproscan.py first.")
            continue

        run_group(group, faa_path, out_tsv, progress)

    print("\nDone.")
