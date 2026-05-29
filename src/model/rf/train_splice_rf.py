"""Option 2 splice-score recalibrator.

Trains a RandomForest on the single CNN-logit feature to map raw donor/acceptor
logits to calibrated log-odds, replacing the hand-tuned scale/bias sweep. Reads
the sparse CNN score TSVs produced by train_splice_cnn_scores.py and writes
recalibrated TSVs in the identical `position\tdonor\tacceptor` format so the C++
HMM consumes them with no change.
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
    read_fasta,
    splice_sites_from_gff,
)

LOGIT_CLAMP = 15.0  # log-odds bound, keeps log(p/(1-p)) finite for the HMM


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def load_scores(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    positions: list[int] = []
    donor: list[float] = []
    acceptor: list[float] = []
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 3:
                continue
            positions.append(int(fields[0]))
            donor.append(float(fields[1]))
            acceptor.append(float(fields[2]))
    return (
        np.array(positions, dtype=np.int64),
        np.array(donor, dtype=np.float64),
        np.array(acceptor, dtype=np.float64),
    )


def labels_for_positions(positions: np.ndarray, true_sites: set[int]) -> np.ndarray:
    return np.fromiter((pos in true_sites for pos in positions), dtype=np.int8, count=len(positions))


def calibrated_log_odds(model: RandomForestClassifier, logits: np.ndarray) -> np.ndarray:
    proba = model.predict_proba(logits.reshape(-1, 1))[:, 1]
    proba = np.clip(proba, 1e-7, 1.0 - 1e-7)
    return np.clip(np.log(proba / (1.0 - proba)), -LOGIT_CLAMP, LOGIT_CLAMP)


def profile_species(profile: dict) -> list[dict]:
    dataset = profile.get("dataset")
    if not isinstance(dataset, dict):
        raise ValueError("Profile must define a dataset object.")
    species = dataset.get("species")
    if not isinstance(species, list):
        raise ValueError("Profile dataset.species must be a list.")
    return species


def build_dataset(
    species: list[dict],
    score_paths: list[Path],
    fasta_key: str,
    gff_key: str,
    min_intron_bp: int,
    require_3n_cds: bool,
) -> dict[str, np.ndarray]:
    donor_logits: list[np.ndarray] = []
    donor_labels: list[np.ndarray] = []
    acceptor_logits: list[np.ndarray] = []
    acceptor_labels: list[np.ndarray] = []

    for entry, score_path in zip(species, score_paths):
        dataset = read_fasta(Path(entry[fasta_key]))
        donors, acceptors = splice_sites_from_gff(
            Path(entry[gff_key]), dataset.offsets, min_intron_bp, require_3n_cds
        )
        positions, donor_score, acceptor_score = load_scores(score_path)
        donor_logits.append(donor_score)
        donor_labels.append(labels_for_positions(positions, donors))
        acceptor_logits.append(acceptor_score)
        acceptor_labels.append(labels_for_positions(positions, acceptors))

    return {
        "donor_logits": np.concatenate(donor_logits),
        "donor_labels": np.concatenate(donor_labels),
        "acceptor_logits": np.concatenate(acceptor_logits),
        "acceptor_labels": np.concatenate(acceptor_labels),
    }


def report_row(name: str, labels: np.ndarray, cnn_logits: np.ndarray, rf_logits: np.ndarray) -> None:
    cnn_pred = (cnn_logits > 0.0).astype(np.int8)
    rf_pred = (rf_logits > 0.0).astype(np.int8)
    for tag, pred, score in (("CNN", cnn_pred, cnn_logits), ("RF ", rf_pred, rf_logits)):
        precision, recall, f1, _ = precision_recall_fscore_support(
            labels, pred, average="binary", zero_division=0
        )
        auc = roc_auc_score(labels, score) if labels.min() != labels.max() else float("nan")
        print(
            f"  {name:<9} {tag}  P={precision:.4f}  R={recall:.4f}  "
            f"F1={f1:.4f}  AUC={auc:.4f}  (pos={int(labels.sum())}/{len(labels)})"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Recalibrate CNN splice logits with a RandomForest (Option 2).")
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--n-estimators", type=int, default=200)
    parser.add_argument("--max-depth", type=int, default=6)
    parser.add_argument("--min-samples-leaf", type=int, default=50)
    parser.add_argument("--output-suffix", default="_splice_rf_scores.tsv")
    parser.add_argument("--cnn-suffix", default="_splice_cnn_scores.tsv")
    args = parser.parse_args()

    with args.profile.open() as handle:
        profile = json.load(handle)

    species = profile_species(profile)
    splice_cnn = profile.get("splice_cnn", {})
    train_scores = [Path(p) for p in splice_cnn.get("train_scores", [])]
    test_scores = [Path(p) for p in splice_cnn.get("test_scores", [])]
    filters = profile.get("filters", {})
    min_intron_bp = int(filters.get("min_intron_bp", 20))
    require_3n_cds = bool(filters.get("require_3n_cds", True))

    if len(train_scores) != len(species) or len(test_scores) != len(species):
        raise ValueError("splice_cnn train_scores/test_scores must align with dataset.species.")

    log("building training set from CNN train scores")
    train = build_dataset(species, train_scores, "train_fasta", "train_gff", min_intron_bp, require_3n_cds)

    forest_kwargs = dict(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        min_samples_leaf=args.min_samples_leaf,
        class_weight="balanced",
        random_state=7,
        n_jobs=-1,
    )
    log(f"training donor RF on {len(train['donor_labels'])} candidates")
    donor_model = RandomForestClassifier(**forest_kwargs)
    donor_model.fit(train["donor_logits"].reshape(-1, 1), train["donor_labels"])
    log(f"training acceptor RF on {len(train['acceptor_labels'])} candidates")
    acceptor_model = RandomForestClassifier(**forest_kwargs)
    acceptor_model.fit(train["acceptor_logits"].reshape(-1, 1), train["acceptor_labels"])

    log("building held-out test set from CNN test scores")
    test = build_dataset(species, test_scores, "test_fasta", "test_gff", min_intron_bp, require_3n_cds)
    test_donor_rf = calibrated_log_odds(donor_model, test["donor_logits"])
    test_acceptor_rf = calibrated_log_odds(acceptor_model, test["acceptor_logits"])

    print("\n=== Intrinsic splice-candidate metrics (held-out test) ===")
    print("  threshold: log-odds > 0  (CNN logit raw; RF calibrated)")
    report_row("donor", test["donor_labels"], test["donor_logits"], test_donor_rf)
    report_row("acceptor", test["acceptor_labels"], test["acceptor_logits"], test_acceptor_rf)
    print()

    for entry, cnn_path in zip(species, test_scores):
        out_path = Path(str(cnn_path).replace(args.cnn_suffix, args.output_suffix))
        positions, donor_logit, acceptor_logit = load_scores(cnn_path)
        donor_out = calibrated_log_odds(donor_model, donor_logit)
        acceptor_out = calibrated_log_odds(acceptor_model, acceptor_logit)
        with out_path.open("w") as handle:
            for pos, d, a in zip(positions, donor_out, acceptor_out):
                handle.write(f"{pos}\t{d:.8f}\t{a:.8f}\n")
        log(f"wrote recalibrated scores: {out_path} ({len(positions)} candidates)")


if __name__ == "__main__":
    main()
