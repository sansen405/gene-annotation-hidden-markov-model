#!/usr/bin/env python3
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
DATA_ROOT = REPO / "genome_data" / "fungi_diverse"
DL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"


def download_zip(accession: str, dest_zip: Path) -> None:
    url = (
        f"{DL}/{accession}/download"
        "?include_annotation_type=GENOME_FASTA"
        "&include_annotation_type=GENOME_GFF"
    )
    subprocess.run(
        ["curl", "-sS", "-L", "--retry", "4", "--retry-delay", "3",
         "-H", "Accept: application/zip", "-o", str(dest_zip), url],
        check=True,
    )


def extract(zip_path: Path, accession: str, folder: str, full_dir: Path) -> None:
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        fna = next((n for n in names if n.endswith(".fna")), None)
        gff = next((n for n in names if n.endswith(".gff")), None)
        if not fna or not gff:
            raise RuntimeError(f"{accession}: zip missing fna/gff -> {names}")
        full_dir.mkdir(parents=True, exist_ok=True)
        with zf.open(fna) as src, (full_dir / f"{folder}_full.fna").open("wb") as dst:
            shutil.copyfileobj(src, dst)
        with zf.open(gff) as src, (full_dir / f"{folder}_full.gff").open("wb") as dst:
            shutil.copyfileobj(src, dst)


def main() -> None:
    manifest = json.loads((HERE / "manifest.json").read_text())
    todo = [m for m in manifest if m.get("accession")]
    print(f"Downloading {len(todo)} genomes into {DATA_ROOT}\n")
    for i, m in enumerate(todo, 1):
        folder, acc = m["folder"], m["accession"]
        full_dir = DATA_ROOT / folder / "full"
        fna = full_dir / f"{folder}_full.fna"
        gff = full_dir / f"{folder}_full.gff"
        if fna.exists() and gff.exists() and fna.stat().st_size > 0 and gff.stat().st_size > 0:
            print(f"  [{i:2d}/{len(todo)}] skip {folder} (already present)")
            continue
        print(f"  [{i:2d}/{len(todo)}] {folder} <- {acc} ...", flush=True)
        with tempfile.TemporaryDirectory() as tmp:
            zip_path = Path(tmp) / f"{acc}.zip"
            try:
                download_zip(acc, zip_path)
                extract(zip_path, acc, folder, full_dir)
            except Exception as exc:  # noqa: BLE001
                print(f"      FAILED {folder} ({acc}): {exc}", file=sys.stderr)
                continue
        print(
            f"      ok: {fna.stat().st_size/1e6:.1f} MB fna, "
            f"{gff.stat().st_size/1e6:.1f} MB gff"
        )
    print("\nDownload pass complete.")


if __name__ == "__main__":
    main()
