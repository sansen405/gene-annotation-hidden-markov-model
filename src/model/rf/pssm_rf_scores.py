"""PSSM and PSSM+RF splice scoring, as a CNN-free alternative emission source.

Experiment: revert donor/acceptor emissions to a position-specific scoring matrix
(PSSM), then train a RandomForest over the per-column PSSM log-odds, and compare
both against the existing CNN scores. Everything is written in the same sparse
`position\tdonor\tacceptor` TSV format the C++ HMM already consumes, so swapping
emission sources needs no C++ change -- only a different score file per run.

Outputs per test species:
  *_splice_pssm_scores.tsv      summed PSSM log-odds (pure PSSM emission)
  *_splice_pssm_rf_scores.tsv   RF-calibrated log-odds (RF over PSSM features)

Also writes sibling profiles fission_yeasts_pssm.json / _pssm_rf.json with the
test_scores paths swapped, so validation can be run directly against each.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_fscore_support, roc_auc_score

CNN_DIR = Path(__file__).resolve().parent.parent / "cnn"
sys.path.insert(0, str(CNN_DIR))
from train_splice_cnn_scores import (  # noqa: E402
    FastaDataset,
    candidate_splice_positions,
    read_fasta,
    splice_sites_from_gff,
)

BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
LOGIT_CLAMP = 15.0


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def is_donor_candidate(seq: str, pos: int) -> bool:
    return seq[pos : pos + 2] == "GT"


def is_acceptor_candidate(seq: str, pos: int) -> bool:
    return pos >= 1 and seq[pos - 1 : pos + 1] == "AG"


def donor_window(seq: str, pos: int, left: int, right: int) -> list[int]:
    return [base_at(seq, c) for c in range(pos - left, pos + right)]


def acceptor_window(seq: str, pos: int, left: int, right: int) -> list[int]:
    return [base_at(seq, c) for c in range(pos - (left - 1), pos + right + 1)]


def base_at(seq: str, idx: int) -> int:
    if idx < 0 or idx >= len(seq):
        return -1
    return BASE_INDEX.get(seq[idx], -1)


class Pssm:
    """Position-specific log-odds: site frequencies vs background, per column."""

    def __init__(self, width: int, alpha: float) -> None:
        self.width = width
        self.alpha = alpha
        self.site = np.zeros((width, 4), dtype=np.float64)
        self.background = np.zeros((width, 4), dtype=np.float64)
        self.log_odds: np.ndarray | None = None

    def add(self, columns: list[int], is_site: bool) -> None:
        target = self.site if is_site else self.background
        for col, base in enumerate(columns):
            if base >= 0:
                target[col, base] += 1.0

    def finalize(self) -> None:
        site = self.site + self.alpha
        background = self.background + self.alpha
        site /= site.sum(axis=1, keepdims=True)
        background /= background.sum(axis=1, keepdims=True)
        self.log_odds = np.log(site) - np.log(background)

    def column_scores(self, columns: list[int]) -> np.ndarray:
        assert self.log_odds is not None
        out = np.zeros(self.width, dtype=np.float64)
        for col, base in enumerate(columns):
            if base >= 0:
                out[col] = self.log_odds[col, base]
        return out


def profile_species(profile: dict) -> list[dict]:
    dataset = profile.get("dataset")
    if not isinstance(dataset, dict):
        raise ValueError("Profile must define a dataset object.")
    species = dataset.get("species")
    if not isinstance(species, list):
        raise ValueError("Profile dataset.species must be a list.")
    return species


def build_pssms(
    species: list[dict],
    donor_win: tuple[int, int],
    acceptor_win: tuple[int, int],
    min_intron_bp: int,
    require_3n_cds: bool,
    alpha: float,
) -> tuple[Pssm, Pssm]:
    donor_pssm = Pssm(sum(donor_win), alpha)
    acceptor_pssm = Pssm(sum(acceptor_win), alpha)

    for entry in species:
        dataset = read_fasta(Path(entry["train_fasta"]))
        donors, acceptors = splice_sites_from_gff(
            Path(entry["train_gff"]), dataset.offsets, min_intron_bp, require_3n_cds
        )
        for record in dataset.records:
            seq = record.sequence
            for local_pos in candidate_splice_positions(seq):
                global_pos = record.offset + local_pos
                if is_donor_candidate(seq, local_pos):
                    donor_pssm.add(
                        donor_window(seq, local_pos, *donor_win), global_pos in donors
                    )
                if is_acceptor_candidate(seq, local_pos):
                    acceptor_pssm.add(
                        acceptor_window(seq, local_pos, *acceptor_win), global_pos in acceptors
                    )

    donor_pssm.finalize()
    acceptor_pssm.finalize()
    return donor_pssm, acceptor_pssm


def collect_features(
    species: list[dict],
    fasta_key: str,
    gff_key: str,
    donor_pssm: Pssm,
    acceptor_pssm: Pssm,
    donor_win: tuple[int, int],
    acceptor_win: tuple[int, int],
    min_intron_bp: int,
    require_3n_cds: bool,
) -> dict[str, np.ndarray]:
    d_feat: list[np.ndarray] = []
    d_label: list[int] = []
    a_feat: list[np.ndarray] = []
    a_label: list[int] = []

    for entry in species:
        dataset = read_fasta(Path(entry[fasta_key]))
        donors, acceptors = splice_sites_from_gff(
            Path(entry[gff_key]), dataset.offsets, min_intron_bp, require_3n_cds
        )
        for record in dataset.records:
            seq = record.sequence
            for local_pos in candidate_splice_positions(seq):
                global_pos = record.offset + local_pos
                if is_donor_candidate(seq, local_pos):
                    cols = donor_pssm.column_scores(donor_window(seq, local_pos, *donor_win))
                    d_feat.append(np.append(cols, cols.sum()))
                    d_label.append(1 if global_pos in donors else 0)
                if is_acceptor_candidate(seq, local_pos):
                    cols = acceptor_pssm.column_scores(acceptor_window(seq, local_pos, *acceptor_win))
                    a_feat.append(np.append(cols, cols.sum()))
                    a_label.append(1 if global_pos in acceptors else 0)

    return {
        "donor_X": np.array(d_feat, dtype=np.float64),
        "donor_y": np.array(d_label, dtype=np.int8),
        "acceptor_X": np.array(a_feat, dtype=np.float64),
        "acceptor_y": np.array(a_label, dtype=np.int8),
    }


def calibrated_log_odds(model: RandomForestClassifier, features: np.ndarray) -> np.ndarray:
    if len(features) == 0:
        return np.zeros(0, dtype=np.float64)
    proba = model.predict_proba(features)[:, 1]
    proba = np.clip(proba, 1e-7, 1.0 - 1e-7)
    return np.clip(np.log(proba / (1.0 - proba)), -LOGIT_CLAMP, LOGIT_CLAMP)


def load_cnn_scores(path: Path) -> dict[int, tuple[float, float]]:
    scores: dict[int, tuple[float, float]] = {}
    if not path.exists():
        return scores
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) == 3:
                scores[int(fields[0])] = (float(fields[1]), float(fields[2]))
    return scores


def metric_row(name: str, labels: np.ndarray, score: np.ndarray, threshold: float) -> str:
    pred = (score > threshold).astype(np.int8)
    precision, recall, f1, _ = precision_recall_fscore_support(
        labels, pred, average="binary", zero_division=0
    )
    auc = roc_auc_score(labels, score) if labels.min() != labels.max() else float("nan")
    return f"P={precision:.4f}  R={recall:.4f}  F1={f1:.4f}  AUC={auc:.4f}"


def write_variant_profile(base_path: Path, profile: dict, cnn_suffix: str, new_suffix: str, tag: str) -> Path:
    variant = json.loads(json.dumps(profile))
    splice = variant.get("splice_cnn", {})
    for key in ("train_scores", "test_scores"):
        if key in splice:
            splice[key] = [p.replace(cnn_suffix, new_suffix) for p in splice[key]]
    out_path = base_path.with_name(base_path.stem + f"_{tag}" + base_path.suffix)
    with out_path.open("w") as handle:
        json.dump(variant, handle, indent=4)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="PSSM and PSSM+RF splice scoring vs CNN.")
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--n-estimators", type=int, default=300)
    parser.add_argument("--max-depth", type=int, default=8)
    parser.add_argument("--min-samples-leaf", type=int, default=20)
    parser.add_argument("--cnn-suffix", default="_splice_cnn_scores.tsv")
    parser.add_argument("--pssm-suffix", default="_splice_pssm_scores.tsv")
    parser.add_argument("--rf-suffix", default="_splice_pssm_rf_scores.tsv")
    args = parser.parse_args()

    with args.profile.open() as handle:
        profile = json.load(handle)

    species = profile_species(profile)
    emissions = profile.get("emissions", {})
    donor_cfg = emissions.get("DONOR", {})
    acceptor_cfg = emissions.get("ACCEPTOR", {})
    donor_win = (int(donor_cfg.get("window_left", 3)), int(donor_cfg.get("window_right", 6)))
    acceptor_win = (int(acceptor_cfg.get("window_left", 15)), int(acceptor_cfg.get("window_right", 3)))
    filters = profile.get("filters", {})
    min_intron_bp = int(filters.get("min_intron_bp", 20))
    require_3n_cds = bool(filters.get("require_3n_cds", True))
    alpha = float(profile.get("smoothing", {}).get("emission_alpha", 0.1))

    test_scores = [Path(p) for p in profile.get("splice_cnn", {}).get("test_scores", [])]

    log(f"building PSSMs (donor window {donor_win}, acceptor window {acceptor_win}, alpha={alpha})")
    donor_pssm, acceptor_pssm = build_pssms(
        species, donor_win, acceptor_win, min_intron_bp, require_3n_cds, alpha
    )

    log("collecting training features for RF")
    train = collect_features(
        species, "train_fasta", "train_gff", donor_pssm, acceptor_pssm,
        donor_win, acceptor_win, min_intron_bp, require_3n_cds,
    )
    forest_kwargs = dict(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        min_samples_leaf=args.min_samples_leaf,
        class_weight="balanced",
        random_state=7,
        n_jobs=-1,
    )
    log(f"training donor RF on {len(train['donor_y'])} candidates ({int(train['donor_y'].sum())} positive)")
    donor_rf = RandomForestClassifier(**forest_kwargs).fit(train["donor_X"], train["donor_y"])
    log(f"training acceptor RF on {len(train['acceptor_y'])} candidates ({int(train['acceptor_y'].sum())} positive)")
    acceptor_rf = RandomForestClassifier(**forest_kwargs).fit(train["acceptor_X"], train["acceptor_y"])

    log("collecting held-out test features")
    test = collect_features(
        species, "test_fasta", "test_gff", donor_pssm, acceptor_pssm,
        donor_win, acceptor_win, min_intron_bp, require_3n_cds,
    )

    donor_pssm_score = test["donor_X"][:, -1]
    acceptor_pssm_score = test["acceptor_X"][:, -1]
    donor_rf_score = calibrated_log_odds(donor_rf, test["donor_X"])
    acceptor_rf_score = calibrated_log_odds(acceptor_rf, test["acceptor_X"])

    print("\n=== Intrinsic held-out splice-candidate metrics ===")
    print("    (PSSM/RF threshold log-odds>0; AUC is threshold-free ranking quality)")
    print(f"  donor    PSSM   {metric_row('', test['donor_y'], donor_pssm_score, 0.0)}")
    print(f"  donor    RF     {metric_row('', test['donor_y'], donor_rf_score, 0.0)}")
    print(f"  acceptor PSSM   {metric_row('', test['acceptor_y'], acceptor_pssm_score, 0.0)}")
    print(f"  acceptor RF     {metric_row('', test['acceptor_y'], acceptor_rf_score, 0.0)}")
    print()

    log("writing PSSM and PSSM+RF score TSVs per test species")
    for entry, cnn_path in zip(species, test_scores):
        dataset = read_fasta(Path(entry["test_fasta"]))
        positions: list[int] = []
        donor_rows: list[np.ndarray] = []
        donor_idx: list[int] = []
        acceptor_rows: list[np.ndarray] = []
        acceptor_idx: list[int] = []
        for record in dataset.records:
            seq = record.sequence
            for local_pos in candidate_splice_positions(seq):
                row = len(positions)
                positions.append(record.offset + local_pos)
                if is_donor_candidate(seq, local_pos):
                    cols = donor_pssm.column_scores(donor_window(seq, local_pos, *donor_win))
                    donor_rows.append(np.append(cols, cols.sum()))
                    donor_idx.append(row)
                if is_acceptor_candidate(seq, local_pos):
                    cols = acceptor_pssm.column_scores(acceptor_window(seq, local_pos, *acceptor_win))
                    acceptor_rows.append(np.append(cols, cols.sum()))
                    acceptor_idx.append(row)

        n = len(positions)
        donor_pssm_out = np.zeros(n)
        donor_rf_out = np.zeros(n)
        acceptor_pssm_out = np.zeros(n)
        acceptor_rf_out = np.zeros(n)
        if donor_rows:
            donor_matrix = np.array(donor_rows)
            donor_pssm_out[donor_idx] = donor_matrix[:, -1]
            donor_rf_out[donor_idx] = calibrated_log_odds(donor_rf, donor_matrix)
        if acceptor_rows:
            acceptor_matrix = np.array(acceptor_rows)
            acceptor_pssm_out[acceptor_idx] = acceptor_matrix[:, -1]
            acceptor_rf_out[acceptor_idx] = calibrated_log_odds(acceptor_rf, acceptor_matrix)

        pssm_path = Path(str(cnn_path).replace(args.cnn_suffix, args.pssm_suffix))
        rf_path = Path(str(cnn_path).replace(args.cnn_suffix, args.rf_suffix))
        with pssm_path.open("w") as pssm_out, rf_path.open("w") as rf_out:
            for i, global_pos in enumerate(positions):
                pssm_out.write(f"{global_pos}\t{donor_pssm_out[i]:.8f}\t{acceptor_pssm_out[i]:.8f}\n")
                rf_out.write(f"{global_pos}\t{donor_rf_out[i]:.8f}\t{acceptor_rf_out[i]:.8f}\n")
        log(f"  wrote {pssm_path.name} and {rf_path.name}")

    pssm_profile = write_variant_profile(args.profile, profile, args.cnn_suffix, args.pssm_suffix, "pssm")
    rf_profile = write_variant_profile(args.profile, profile, args.cnn_suffix, args.rf_suffix, "pssm_rf")
    log(f"wrote variant profiles: {pssm_profile.name}, {rf_profile.name}")
    print("\nValidate each emission source with the SAME binary, swapping only the profile:")
    print(f"  /tmp/full_genome_validation --profile {args.profile}      # CNN baseline")
    print(f"  /tmp/full_genome_validation --profile {pssm_profile}   # pure PSSM")
    print(f"  /tmp/full_genome_validation --profile {rf_profile}   # PSSM + RF")


if __name__ == "__main__":
    main()
