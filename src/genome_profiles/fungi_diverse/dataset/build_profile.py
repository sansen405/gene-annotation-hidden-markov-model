#!/usr/bin/env python3
from __future__ import annotations

import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
DATA = "genome_data/fungi_diverse"
PROFILE = REPO / "src/genome_profiles/fungi_diverse/fungi_diverse.json"

ORGANELLE_KEYWORDS = (
    "mitochond", "chloroplast", "plastid", "apicoplast", "kinetoplast",
)
TEST_TARGET_FRACTION = 0.10
TEST_MAX_SINGLE = 0.25
SEED = 1729


def fasta_headers(path: Path) -> list[tuple[str, str]]:
    out = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                rest = line[1:].rstrip("\n")
                seqid = rest.split()[0]
                out.append((seqid, rest))
    return out


def cds_counts(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    with path.open() as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.split("\t")
            if len(parts) > 4 and parts[2] == "CDS":
                counts[parts[0]] = counts.get(parts[0], 0) + 1
    return counts


def is_organelle(description: str) -> bool:
    low = description.lower()
    return any(k in low for k in ORGANELLE_KEYWORDS)


def choose_test(seqids: list[str], cds: dict[str, int]) -> list[str]:
    total = sum(cds.get(s, 0) for s in seqids)
    if total == 0:
        return []
    rng = random.Random(SEED)
    order = list(seqids)
    rng.shuffle(order)
    target = TEST_TARGET_FRACTION * total
    cap = TEST_MAX_SINGLE * total
    picked, acc = [], 0
    for s in order:
        c = cds.get(s, 0)
        if c == 0 or c > cap:
            continue
        picked.append(s)
        acc += c
        if acc >= target:
            break
    if not picked:
        candidates = sorted((s for s in seqids if cds.get(s, 0) > 0), key=lambda s: cds[s])
        if candidates:
            picked = [candidates[0]]
    return picked


def species_entry(folder: str) -> dict:
    base = REPO / DATA / folder
    full_fna = base / "full" / f"{folder}_full.fna"
    full_gff = base / "full" / f"{folder}_full.gff"
    headers = fasta_headers(full_fna)
    cds = cds_counts(full_gff)

    organelles = [sid for sid, desc in headers if is_organelle(desc)]
    nuclear = [sid for sid, _ in headers if sid not in set(organelles)]
    test = choose_test(nuclear, cds)

    total_cds = sum(cds.values())
    test_cds = sum(cds.get(s, 0) for s in test)
    frac = (test_cds / total_cds) if total_cds else 0.0
    print(
        f"  {folder:16s} seqs={len(headers):4d} organelle={len(organelles)} "
        f"CDS={total_cds:6d} test_scaffolds={len(test):3d} "
        f"test_CDS={test_cds:5d} ({frac*100:4.1f}%)"
    )

    d = f"{DATA}/{folder}"
    return {
        "name": folder,
        "source_fasta": f"{d}/full/{folder}_full.fna",
        "source_gff":   f"{d}/full/{folder}_full.gff",
        "train_fasta":  f"{d}/train/{folder}_train.fna",
        "train_gff":    f"{d}/train/{folder}_train.gff",
        "test_fasta":   f"{d}/test/{folder}_test.fna",
        "test_gff":     f"{d}/test/{folder}_test.gff",
        "test_chromosomes": test,
        "excluded_chromosomes": organelles,
    }


def main() -> None:
    manifest = json.loads((HERE / "manifest.json").read_text())
    folders = [m["folder"] for m in manifest if m.get("accession")]
    print(f"Scanning {len(folders)} genomes to build holdouts...\n")
    species = [species_entry(f) for f in folders]

    def score_paths(kind: str, signal: str) -> list[str]:
        return [
            f"{DATA}/{s['name']}/{kind}/{s['name']}_{signal}_cnn_scores.tsv"
            for s in species
        ]

    def score_paths_minus(kind: str, signal: str) -> list[str]:
        return [
            f"{DATA}/{s['name']}/{kind}/{s['name']}_{signal}_cnn_scores_minus.tsv"
            for s in species
        ]

    profile = {
        "name": "fungi_diverse",
        "description": (
            "Diverse pan-fungal benchmark: 24 RefSeq-annotated genomes spanning "
            "Saccharomycotina, Pezizomycotina, Basidiomycota, Mucoromycota and "
            "Chytridiomycota, chosen to span the intron-density gradient that "
            "stresses the splice/intron model."
        ),
        "dataset": {"species": species},
        "splice_cnn": {
            "model": "src/model/cnn/splice/trained_models/fungi_diverse_splice_cnn.pt",
            "donor_scale": 1.0,
            "donor_bias": 0.0,
            "acceptor_scale": 1.0,
            "acceptor_bias": 0.0,
            "train_scores": score_paths("train", "splice"),
            "test_scores": score_paths("test", "splice"),
            "train_scores_minus": score_paths_minus("train", "splice"),
            "test_scores_minus": score_paths_minus("test", "splice"),
        },
        "start_cnn": {
            "model": "src/model/cnn/start/trained_models/fungi_diverse_start_cnn.pt",
            "start_scale": 1.0,
            "start_bias": 0.0,
            "train_scores": score_paths("train", "start"),
            "test_scores": score_paths("test", "start"),
            "train_scores_minus": score_paths_minus("train", "start"),
            "test_scores_minus": score_paths_minus("test", "start"),
        },
        "filters": {
            "min_first_cds_bp": 3,
            "min_last_cds_bp": 3,
            "min_intron_bp": 20,
            "require_3n_cds": True,
            "include_minus_strand": True,
        },
        "emissions": {
            "INTERGENIC":  {"type": "markov", "order": 1},
            "INTRON":      {"type": "markov", "order": 1},
            "EXON_FRAME":  {"type": "markov", "order": 5, "frame_tied": True},
            "DONOR":       {"type": "pssm", "window_left": 3, "window_right": 6},
            "ACCEPTOR":    {"type": "pssm", "window_left": 15, "window_right": 3},
            "START_CODON": {"type": "deterministic"},
            "STOP_CODON":  {"type": "deterministic"},
        },
        "smoothing": {"transition_alpha": 0.02, "emission_alpha": 0.1},
    }

    PROFILE.write_text(json.dumps(profile, indent=4) + "\n")
    n_test = sum(len(s["test_chromosomes"]) for s in species)
    n_excl = sum(len(s["excluded_chromosomes"]) for s in species)
    print(
        f"\nWrote {PROFILE}\n"
        f"  species={len(species)}  held-out test scaffolds={n_test}  "
        f"excluded organelles={n_excl}"
    )


if __name__ == "__main__":
    main()
