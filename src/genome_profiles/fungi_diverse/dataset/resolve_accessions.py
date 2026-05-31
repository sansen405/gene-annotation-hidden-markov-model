#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon"


def read_species(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        folder, taxon, group = (c.strip() for c in line.split("\t"))
        rows.append({"folder": folder, "taxon": taxon, "group": group})
    return rows


def fetch_reports(taxon: str, source: str) -> list[dict]:
    params = {
        "filters.assembly_source": source,
        "filters.has_annotation": "true",
        "page_size": "50",
    }
    url = f"{API}/{urllib.parse.quote(taxon)}/dataset_report?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.load(resp)
            return data.get("reports", []) or []
        except Exception as exc:  # noqa: BLE001
            if attempt == 3:
                print(f"  ! API error for {taxon} ({source}): {exc}", file=sys.stderr)
                return []
            time.sleep(2 * (attempt + 1))
    return []


def rank(report: dict) -> tuple:
    info = report.get("assembly_info", {})
    category = (info.get("refseq_category") or "").lower()
    level = (info.get("assembly_level") or "").lower()
    cat_rank = {"reference genome": 0, "representative genome": 1}.get(category, 2)
    level_rank = {
        "complete genome": 0,
        "chromosome": 1,
        "scaffold": 2,
        "contig": 3,
    }.get(level, 4)
    return (cat_rank, level_rank)


def pick(taxon: str) -> dict | None:
    for source in ("refseq", "genbank"):
        reports = fetch_reports(taxon, source)
        annotated = [
            r for r in reports
            if r.get("annotation_info", {}).get("name")
        ]
        if not annotated:
            continue
        best = sorted(annotated, key=rank)[0]
        info = best.get("assembly_info", {})
        stats = best.get("assembly_stats", {})
        return {
            "accession": best.get("accession"),
            "organism": best.get("organism", {}).get("organism_name"),
            "source": source,
            "assembly_level": info.get("assembly_level"),
            "refseq_category": info.get("refseq_category"),
            "annotation": best.get("annotation_info", {}).get("name"),
            "total_sequence_length": stats.get("total_sequence_length"),
            "number_of_contigs": stats.get("number_of_contigs"),
        }
    return None


def main() -> None:
    species = read_species(HERE / "species.tsv")
    manifest = []
    total_bp = 0
    print(f"Resolving {len(species)} species via NCBI datasets API...\n")
    for row in species:
        chosen = pick(row["taxon"])
        if chosen is None:
            print(f"  MISSING  {row['folder']:16s} {row['taxon']} (no annotated assembly)")
            manifest.append({**row, "accession": None})
            continue
        bp = int(chosen["total_sequence_length"] or 0)
        total_bp += bp
        print(
            f"  {chosen['accession']:20s} {row['folder']:16s} "
            f"{bp/1e6:6.1f} Mb  {chosen['assembly_level']:16s} "
            f"[{chosen['source']}] {row['group']}"
        )
        manifest.append({**row, **chosen})
        time.sleep(0.34)

    resolved = [m for m in manifest if m.get("accession")]
    out = HERE / "manifest.json"
    out.write_text(json.dumps(manifest, indent=2))
    print(
        f"\nResolved {len(resolved)}/{len(species)} species, "
        f"total ~{total_bp/1e6:.0f} Mb of sequence."
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
