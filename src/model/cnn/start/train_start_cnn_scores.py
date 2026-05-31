from __future__ import annotations

import argparse
import json
import math
import random
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import torch
    from start_cnn_network import StartCNN


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str
    offset: int


@dataclass(frozen=True)
class FastaDataset:
    path: Path
    records: list[FastaRecord]
    offsets: dict[str, int]

    @property
    def sequence(self) -> str:
        return "".join(record.sequence for record in self.records)

    @property
    def base_count(self) -> int:
        return sum(len(record.sequence) for record in self.records)


@dataclass(frozen=True)
class Calibration:
    """Score-time transform that turns raw start logits into genomic-prior log-odds.

    Mirrors the splice calibration: `temperature` fixes the scale and subtracting
    the training-prior logit fixes the offset, so the written scores are
    likelihood-ratio-style log-odds and the HMM can decode the start emission at
    scale=1, bias=0.
    """

    temperature: float
    start_prior_logit: float

    def apply(self, start_logit: float) -> float:
        return start_logit / self.temperature - self.start_prior_logit


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def format_count(value: int) -> str:
    return f"{value:,}"


_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def reverse_complement(sequence: str) -> str:
    return "".join(_COMPLEMENT.get(base, "N") for base in reversed(sequence))


def read_fasta(path: Path) -> FastaDataset:
    records: list[FastaRecord] = []
    offsets: dict[str, int] = {}
    current_name: str | None = None
    current_parts: list[str] = []
    current_offset = 0

    def flush_record() -> None:
        nonlocal current_name, current_parts, current_offset
        if current_name is None:
            return
        sequence = "".join(current_parts)
        records.append(FastaRecord(current_name, sequence, current_offset))
        current_offset += len(sequence)
        current_parts = []

    log(f"loading FASTA: {path}")
    start_time = time.perf_counter()
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_record()
                current_name = line[1:].split()[0]
                offsets[current_name] = current_offset
                continue
            current_parts.append(line.upper())

    flush_record()
    dataset = FastaDataset(path=path, records=records, offsets=offsets)
    log(
        f"loaded FASTA: {path} "
        f"records={format_count(len(dataset.records))} bases={format_count(dataset.base_count)} "
        f"elapsed={time.perf_counter() - start_time:.1f}s"
    )
    return dataset


def parent_id(attributes: str) -> str:
    for field in attributes.split(";"):
        if field.startswith("Parent="):
            return field.removeprefix("Parent=")
    return ""


def start_sites_from_gff(
    gff_path: Path,
    offsets: dict[str, int],
    require_3n_cds: bool,
) -> dict[int, str]:
    """Return translation-start positions mapped to their transcript strand.

    Both strands are kept. The start is the first coding base of the transcript:
    the lowest-coordinate CDS base on the plus strand and the highest-coordinate
    CDS base on the minus strand. Minus-strand sites are reported at their plus-
    genome coordinate but tagged "-" so the caller can reverse-complement the
    window and present every true start to the model in canonical 5'->3' (ATG)
    orientation.
    """
    cds_by_parent: dict[str, list[tuple[int, int]]] = defaultdict(list)
    strand_by_parent: dict[str, str] = {}

    log(f"loading start sites from GFF: {gff_path}")
    start_time = time.perf_counter()
    with gff_path.open() as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            chrom, _, _, start, end, _, strand, _, attributes = fields
            if chrom not in offsets:
                continue

            parent = parent_id(attributes)
            if not parent:
                continue

            key = f"{chrom}|{parent}"
            start_0 = offsets[chrom] + int(start) - 1
            end_0 = offsets[chrom] + int(end) - 1
            cds_by_parent[key].append((start_0, end_0))
            strand_by_parent[key] = strand

    starts: dict[int, str] = {}
    for key, fragments in cds_by_parent.items():
        fragments.sort()
        strand = strand_by_parent.get(key, "+")
        if strand not in ("+", "-"):
            continue
        total_cds = sum(end - start + 1 for start, end in fragments)
        if require_3n_cds and total_cds % 3 != 0:
            continue

        if strand == "+":
            # Plus strand: first coding base is the A of the ATG (low coord).
            starts[fragments[0][0]] = "+"
        else:
            # Minus strand: first coding base is the high-coordinate end; the
            # window is reverse-complemented later so it reads ATG canonically.
            starts[fragments[-1][1]] = "-"

    log(
        f"loaded start sites: starts={format_count(len(starts))} "
        f"elapsed={time.perf_counter() - start_time:.1f}s"
    )
    return starts


def window_at(sequence: str, center: int, radius: int) -> str:
    chars: list[str] = []
    for pos in range(center - radius, center + radius + 1):
        if pos < 0 or pos >= len(sequence):
            chars.append("N")
        else:
            chars.append(sequence[pos])
    return "".join(chars)


def candidate_start_positions(sequence: str) -> list[int]:
    positions: list[int] = []
    for pos in range(0, len(sequence) - 2):
        if sequence[pos : pos + 3] == "ATG":
            positions.append(pos)
    return positions


def sample_training_examples(
    dataset: FastaDataset,
    starts: dict[int, str],
    radius: int,
    negatives_per_positive: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    log(f"sampling examples from {dataset.path}")
    start_time = time.perf_counter()
    positives = sorted(set(starts))
    positive_set = set(positives)
    candidate_negatives: list[tuple[FastaRecord, int, int]] = []
    positive_sites_by_record: dict[str, list[int]] = defaultdict(list)
    total_bases = dataset.base_count
    scanned_bases = 0
    next_progress = 10

    for record in dataset.records:
        record_start = record.offset
        record_end = record.offset + len(record.sequence)
        positive_sites_by_record[record.name] = [
            pos for pos in positives if record_start <= pos < record_end
        ]
        for local_pos in candidate_start_positions(record.sequence):
            global_pos = record.offset + local_pos
            if global_pos in positive_set:
                continue
            candidate_negatives.append((record, local_pos, global_pos))
        scanned_bases += len(record.sequence)
        progress = int((scanned_bases / max(1, total_bases)) * 100)
        if progress >= next_progress:
            log(
                f"scanned negative candidates: {min(progress, 100)}% "
                f"candidates={format_count(len(candidate_negatives))}"
            )
            next_progress += 10

    rng = random.Random(7)
    negative_count = min(len(candidate_negatives), len(positives) * negatives_per_positive)
    negatives = rng.sample(candidate_negatives, negative_count)
    log(
        f"selected examples: positives={format_count(len(positives))} "
        f"negatives={format_count(len(negatives))}"
    )

    windows: list[str] = []
    labels: list[list[float]] = []
    for record in dataset.records:
        for pos in positive_sites_by_record[record.name]:
            local_pos = pos - record.offset
            window = window_at(record.sequence, local_pos, radius)
            # Minus-strand starts are presented reverse-complemented so the model
            # always sees a canonical ATG window, regardless of source strand.
            if starts.get(pos) == "-":
                window = reverse_complement(window)
            windows.append(window)
            labels.append([1.0])
    for record, local_pos, _ in negatives:
        windows.append(window_at(record.sequence, local_pos, radius))
        labels.append([0.0])

    if not windows:
        window_size = radius * 2 + 1
        return torch.empty((0, 4, window_size), dtype=torch.float32), torch.empty((0, 1), dtype=torch.float32)

    log(f"encoding windows: count={format_count(len(windows))} width={radius * 2 + 1}")
    inputs = one_hot_encode_windows(windows)
    label_tensor = torch.tensor(labels, dtype=torch.float32)
    log(f"finished examples from {dataset.path} elapsed={time.perf_counter() - start_time:.1f}s")
    return inputs, label_tensor


def bce_loss(logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
    """Plain binary cross-entropy: a proper scoring rule, so the trained logits are
    calibrated log-odds (the same choice the splice trainer uses)."""
    return nn.functional.binary_cross_entropy_with_logits(logits, targets)


def fit_temperature(logits: torch.Tensor, labels: torch.Tensor) -> float:
    """Temperature scaling: fit one scalar on held-out logits so predicted
    probabilities match observed frequencies (automatic replacement for a
    hand-tuned emission scale)."""
    log_temp = torch.zeros(1, requires_grad=True)
    optimizer = torch.optim.LBFGS([log_temp], lr=0.1, max_iter=50)

    def closure() -> torch.Tensor:
        optimizer.zero_grad()
        loss = nn.functional.binary_cross_entropy_with_logits(logits / log_temp.exp(), labels)
        loss.backward()
        return loss

    optimizer.step(closure)
    return float(log_temp.exp().item())


def average_precision(scores: torch.Tensor, targets: torch.Tensor) -> float:
    """Area under the precision-recall curve for the start column. Threshold-free,
    so it reflects start discrimination on the imbalanced ATG candidate set.
    Returns NaN if no positives."""
    n_pos = float(targets.sum().item())
    if n_pos == 0:
        return float("nan")
    order = torch.argsort(scores, descending=True)
    sorted_targets = targets[order]
    tp = torch.cumsum(sorted_targets, dim=0)
    fp = torch.cumsum(1.0 - sorted_targets, dim=0)
    precision = tp / (tp + fp)
    return float((precision * sorted_targets).sum().item() / n_pos)


# MPS conv1d rejects very large single-forward batches ("Output channels > 65536"),
# so evaluation/calibration forwards run in mini-batches. The training loop already
# batches via the DataLoader.
EVAL_BATCH_SIZE = 4096


def forward_in_batches(model, inputs, device, batch_size: int = EVAL_BATCH_SIZE):
    outputs = []
    with torch.no_grad():
        for start in range(0, len(inputs), batch_size):
            chunk = inputs[start : start + batch_size].to(device)
            outputs.append(model(chunk).cpu())
    if not outputs:
        return torch.empty((0, 1), dtype=torch.float32)
    return torch.cat(outputs)


def detect_device() -> "torch.device":
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def train_model(
    model: StartCNN,
    inputs: torch.Tensor,
    labels: torch.Tensor,
    epochs: int,
    batch_size: int,
    device: "torch.device",
) -> Calibration:
    # Seeded shuffle before splitting so the held-out 10 % is representative
    # (examples are concatenated positives-first, then per species).
    generator = torch.Generator().manual_seed(13)
    perm = torch.randperm(len(inputs), generator=generator)
    inputs = inputs[perm]
    labels = labels[perm]

    n_val = max(1, len(inputs) // 10)
    n_train = len(inputs) - n_val
    train_inputs, val_inputs = inputs[:n_train], inputs[n_train:]
    train_labels, val_labels = labels[:n_train], labels[n_train:]

    loader = DataLoader(TensorDataset(train_inputs, train_labels), batch_size=batch_size, shuffle=True)
    model = model.to(device)
    # Val tensors stay on CPU; evaluation runs in mini-batches (forward_in_batches)
    # because MPS conv1d rejects a single forward over the full held-out set.

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)
    total_batches = len(loader)

    log(
        f"starting start CNN training: train={format_count(n_train)} val={format_count(n_val)} "
        f"epochs={epochs} batch_size={batch_size} batches_per_epoch={format_count(total_batches)} "
        f"device={device}"
    )
    model.train()
    for epoch in range(epochs):
        epoch_start = time.perf_counter()
        total_loss = 0.0
        progress_interval = max(1, total_batches // 10)
        for batch_index, (batch_inputs, batch_labels) in enumerate(loader, start=1):
            batch_inputs = batch_inputs.to(device)
            batch_labels = batch_labels.to(device)
            optimizer.zero_grad()
            loss = bce_loss(model(batch_inputs), batch_labels)
            loss.backward()
            optimizer.step()
            total_loss += float(loss.item()) * batch_inputs.size(0)
            if batch_index == total_batches or batch_index % progress_interval == 0:
                log(f"epoch {epoch + 1}/{epochs}: batch {batch_index}/{total_batches}")
        scheduler.step()
        model.eval()
        val_logits = forward_in_batches(model, val_inputs, device)
        val_loss = float(bce_loss(val_logits, val_labels).item())
        val_probs = torch.sigmoid(val_logits)
        val_targets = val_labels
        start_ap = average_precision(val_probs[:, 0], val_targets[:, 0])
        model.train()
        log(
            f"finished epoch {epoch + 1}/{epochs}: "
            f"train_loss={total_loss / n_train:.4f} val_loss={val_loss:.4f} "
            f"val_start_ap={start_ap:.4f} "
            f"lr={scheduler.get_last_lr()[0]:.2e} elapsed={time.perf_counter() - epoch_start:.1f}s"
        )

    # Fit the score-time calibration so the HMM can decode at scale=1, bias=0.
    model.eval()
    final_val_logits = forward_in_batches(model, val_inputs, device)

    temperature = fit_temperature(final_val_logits, val_labels)

    def to_logit(p: float) -> float:
        p = min(max(p, 1e-6), 1.0 - 1e-6)
        return math.log(p / (1.0 - p))

    calibration = Calibration(
        temperature=temperature,
        start_prior_logit=to_logit(float(train_labels[:, 0].mean())),
    )
    log(
        f"fitted calibration: temperature={calibration.temperature:.4f} "
        f"start_prior_logit={calibration.start_prior_logit:.4f}"
    )
    return calibration


def write_scores(
    model: StartCNN,
    dataset: FastaDataset,
    radius: int,
    output_path: Path,
    batch_size: int,
    device: "torch.device",
    calibration: Calibration,
    reverse: bool = False,
) -> None:
    model.eval()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Minus-strand scores are produced by scoring each record's reverse complement;
    # positions stay record.offset + local (revcomp-local) to match the decoder.
    log(
        f"writing scores: input={dataset.path} output={output_path} "
        f"bases={format_count(dataset.base_count)} strand={'-' if reverse else '+'}"
    )
    start_time = time.perf_counter()
    written_bases = 0
    next_progress = 10
    with output_path.open("w") as output, torch.no_grad():
        for record in dataset.records:
            sequence = reverse_complement(record.sequence) if reverse else record.sequence
            for start in range(0, len(sequence), batch_size):
                end = min(start + batch_size, len(sequence))
                windows = [window_at(sequence, pos, radius) for pos in range(start, end)]
                logits = model(one_hot_encode_windows(windows).to(device)).cpu()
                for offset, row in enumerate(logits):
                    position = record.offset + start + offset
                    score = calibration.apply(float(row[0]))
                    output.write(f"{position}\t{score:.8f}\n")
                written_bases += end - start
                progress = int((written_bases / max(1, dataset.base_count)) * 100)
                if progress >= next_progress:
                    log(f"score writing progress for {output_path}: {min(progress, 100)}%")
                    next_progress += 10
    log(f"finished scores: {output_path} elapsed={time.perf_counter() - start_time:.1f}s")


def write_sparse_scores(
    model: StartCNN,
    dataset: FastaDataset,
    radius: int,
    output_path: Path,
    batch_size: int,
    device: "torch.device",
    calibration: Calibration,
    reverse: bool = False,
) -> None:
    model.eval()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # When reverse is set, ATG candidates are found on the reverse complement so the
    # sparse minus-strand TSV covers start candidates in revcomp coordinates.
    sequences = [
        reverse_complement(record.sequence) if reverse else record.sequence
        for record in dataset.records
    ]
    candidate_counts = [len(candidate_start_positions(sequence)) for sequence in sequences]
    total_candidates = sum(candidate_counts)
    log(
        f"writing sparse scores: input={dataset.path} output={output_path} "
        f"bases={format_count(dataset.base_count)} candidates={format_count(total_candidates)} "
        f"strand={'-' if reverse else '+'}"
    )

    start_time = time.perf_counter()
    written_candidates = 0
    next_progress = 10
    with output_path.open("w") as output, torch.no_grad():
        for record, sequence, candidate_count in zip(dataset.records, sequences, candidate_counts):
            positions = candidate_start_positions(sequence)
            if len(positions) != candidate_count:
                raise RuntimeError("Candidate position count changed during sparse scoring.")
            for start in range(0, len(positions), batch_size):
                batch_positions = positions[start : start + batch_size]
                windows = [window_at(sequence, pos, radius) for pos in batch_positions]
                logits = model(one_hot_encode_windows(windows).to(device)).cpu()
                for local_pos, row in zip(batch_positions, logits):
                    position = record.offset + local_pos
                    score = calibration.apply(float(row[0]))
                    output.write(f"{position}\t{score:.8f}\n")
                written_candidates += len(batch_positions)
                progress = int((written_candidates / max(1, total_candidates)) * 100)
                if progress >= next_progress:
                    log(f"sparse score writing progress for {output_path}: {min(progress, 100)}%")
                    next_progress += 10
    log(f"finished sparse scores: {output_path} elapsed={time.perf_counter() - start_time:.1f}s")


def validate_matching_counts(label: str, left: list[Path], right: list[Path]) -> None:
    if len(left) != len(right):
        raise ValueError(f"Expected the same number of {label} FASTA and GFF paths.")


def validate_output_counts(label: str, inputs: list[Path], outputs: list[Path]) -> None:
    if len(inputs) != len(outputs):
        raise ValueError(f"Expected one {label} score output path per {label} FASTA.")


def list_from_json(value: object, label: str) -> list[Path]:
    if isinstance(value, str):
        return [Path(value)]
    if isinstance(value, list) and all(isinstance(item, str) for item in value):
        return [Path(item) for item in value]
    raise ValueError(f"Profile start_cnn.{label} must be a string or list of strings.")


def profile_species(profile: dict[str, object]) -> list[dict[str, object]]:
    dataset = profile.get("dataset")
    if not isinstance(dataset, dict):
        raise ValueError("Profile must define a dataset object.")

    species = dataset.get("species")
    if species is None:
        return [dataset]
    if isinstance(species, list) and all(isinstance(item, dict) for item in species):
        return species
    raise ValueError("Profile dataset.species must be a list of dataset objects.")


def path_list_from_species(species: list[dict[str, object]], key: str) -> list[Path]:
    paths: list[Path] = []
    for item in species:
        value = item.get(key)
        if not isinstance(value, str):
            raise ValueError(f"Profile species entry is missing string field: {key}")
        paths.append(Path(value))
    return paths


def apply_profile_defaults(args: argparse.Namespace) -> None:
    if args.profile is None:
        return

    with args.profile.open() as handle:
        profile = json.load(handle)
    if not isinstance(profile, dict):
        raise ValueError("Profile JSON must be an object.")

    species = profile_species(profile)
    start_cnn = profile.get("start_cnn")
    if not isinstance(start_cnn, dict):
        raise ValueError("Profile must define start_cnn paths for start CNN score generation.")

    args.train_fasta = args.train_fasta or path_list_from_species(species, "train_fasta")
    args.train_gff = args.train_gff or path_list_from_species(species, "train_gff")
    args.test_fasta = args.test_fasta or path_list_from_species(species, "test_fasta")
    if args.model_out is None:
        model_path = start_cnn.get("model")
        if not isinstance(model_path, str) or not model_path:
            raise ValueError("Profile start_cnn.model must be a non-empty string.")
        args.model_out = Path(model_path)
    args.train_scores_out = args.train_scores_out or list_from_json(start_cnn.get("train_scores"), "train_scores")
    args.test_scores_out = args.test_scores_out or list_from_json(start_cnn.get("test_scores"), "test_scores")
    # Minus-strand score paths are optional; only generated when configured.
    if args.train_scores_minus_out is None and start_cnn.get("train_scores_minus") is not None:
        args.train_scores_minus_out = list_from_json(start_cnn.get("train_scores_minus"), "train_scores_minus")
    if args.test_scores_minus_out is None and start_cnn.get("test_scores_minus") is not None:
        args.test_scores_minus_out = list_from_json(start_cnn.get("test_scores_minus"), "test_scores_minus")

    filters = profile.get("filters")
    if isinstance(filters, dict):
        if args.require_3n_cds is None:
            args.require_3n_cds = bool(filters.get("require_3n_cds", True))


def validate_required_args(args: argparse.Namespace) -> None:
    missing: list[str] = []
    for attr, option in [
        ("train_fasta", "--train-fasta"),
        ("train_gff", "--train-gff"),
        ("test_fasta", "--test-fasta"),
        ("model_out", "--model-out"),
        ("train_scores_out", "--train-scores-out"),
        ("test_scores_out", "--test-scores-out"),
    ]:
        if not getattr(args, attr):
            missing.append(option)
    if missing:
        raise ValueError("Provide --profile or explicit values for: " + ", ".join(missing))


def main() -> None:
    parser = argparse.ArgumentParser(description="Train a translation-start CNN and write HMM emission scores.")
    parser.add_argument("--profile", type=Path, help="Genome profile JSON that defines dataset and start_cnn paths.")
    parser.add_argument("--train-fasta", type=Path, nargs="+")
    parser.add_argument("--train-gff", type=Path, nargs="+")
    parser.add_argument("--test-fasta", type=Path, nargs="+")
    parser.add_argument("--model-out", type=Path)
    parser.add_argument("--train-scores-out", type=Path, nargs="+")
    parser.add_argument("--test-scores-out", type=Path, nargs="+")
    parser.add_argument("--train-scores-minus-out", type=Path, nargs="+")
    parser.add_argument("--test-scores-minus-out", type=Path, nargs="+")
    parser.add_argument("--radius", default=60, type=int)
    parser.add_argument("--epochs", default=10, type=int)
    parser.add_argument("--batch-size", default=256, type=int)
    parser.add_argument("--score-batch-size", default=8192, type=int)
    parser.add_argument("--sparse-scores", action="store_true")
    parser.add_argument("--negatives-per-positive", default=5, type=int)
    parser.add_argument("--require-3n-cds", dest="require_3n_cds", action="store_true")
    parser.add_argument("--allow-non-3n-cds", dest="require_3n_cds", action="store_false")
    parser.set_defaults(require_3n_cds=None)
    args = parser.parse_args()

    apply_profile_defaults(args)
    if args.require_3n_cds is None:
        args.require_3n_cds = True
    validate_required_args(args)

    global torch, nn, DataLoader, TensorDataset, StartCNN, one_hot_encode_windows
    import torch
    from start_cnn_network import StartCNN, one_hot_encode_windows
    from torch import nn
    from torch.utils.data import DataLoader, TensorDataset

    validate_matching_counts("training", args.train_fasta, args.train_gff)
    validate_output_counts("training", args.train_fasta, args.train_scores_out)
    validate_output_counts("test", args.test_fasta, args.test_scores_out)

    device = detect_device()
    log(f"device: {device}")

    log("starting start CNN score pipeline")
    log(f"model checkpoint: {args.model_out}")
    train_datasets = [read_fasta(path) for path in args.train_fasta]
    test_datasets = [read_fasta(path) for path in args.test_fasta]
    model = StartCNN(window_size=args.radius * 2 + 1)

    if args.model_out.exists():
        log(f"loading existing CNN checkpoint: {args.model_out}")
        checkpoint = torch.load(args.model_out, map_location="cpu")
        model.load_state_dict(checkpoint["state_dict"])
        model = model.to(device)
        calibration = Calibration(**checkpoint["calibration"])
        log(
            f"loaded calibration: temperature={calibration.temperature:.4f} "
            f"start_prior_logit={calibration.start_prior_logit:.4f}"
        )
    else:
        all_starts: dict[int, str] = {}
        input_batches: list[torch.Tensor] = []
        label_batches: list[torch.Tensor] = []
        for dataset, gff_path in zip(train_datasets, args.train_gff):
            starts = start_sites_from_gff(
                gff_path, dataset.offsets,
                require_3n_cds=args.require_3n_cds,
            )
            all_starts.update(starts)
            inputs, labels = sample_training_examples(
                dataset, starts, args.radius, args.negatives_per_positive,
            )
            if len(labels) > 0:
                input_batches.append(inputs)
                label_batches.append(labels)
        if not input_batches:
            raise ValueError("No start-site examples were found in the training inputs.")
        log(f"training totals: starts={format_count(len(all_starts))}")
        log("combining training tensors")
        inputs = torch.cat(input_batches)
        labels = torch.cat(label_batches)
        calibration = train_model(model, inputs, labels, args.epochs, args.batch_size, device)
        args.model_out.parent.mkdir(parents=True, exist_ok=True)
        log(f"saving CNN checkpoint: {args.model_out}")
        # Save on CPU so the checkpoint is device-agnostic; calibration travels
        # with the weights so reloading reproduces identical scores.
        torch.save(
            {
                "state_dict": {k: v.cpu() for k, v in model.state_dict().items()},
                "calibration": {
                    "temperature": calibration.temperature,
                    "start_prior_logit": calibration.start_prior_logit,
                },
            },
            args.model_out,
        )

    def emit(dataset: FastaDataset, score_path: Path, reverse: bool) -> None:
        if args.sparse_scores:
            write_sparse_scores(model, dataset, args.radius, score_path, args.score_batch_size, device, calibration, reverse)
        else:
            write_scores(model, dataset, args.radius, score_path, args.score_batch_size, device, calibration, reverse)

    # Scores are written as calibrated log-odds, so the HMM can decode the start
    # emission at scale=1, bias=0 with no grid tuning.
    for dataset, score_path in zip(train_datasets, args.train_scores_out):
        emit(dataset, score_path, reverse=False)
    for dataset, score_path in zip(test_datasets, args.test_scores_out):
        emit(dataset, score_path, reverse=False)

    # Minus-strand (reverse-complement) scores, when configured for dual-strand decode.
    if args.train_scores_minus_out:
        validate_output_counts("training minus", args.train_fasta, args.train_scores_minus_out)
        for dataset, score_path in zip(train_datasets, args.train_scores_minus_out):
            emit(dataset, score_path, reverse=True)
    if args.test_scores_minus_out:
        validate_output_counts("test minus", args.test_fasta, args.test_scores_minus_out)
        for dataset, score_path in zip(test_datasets, args.test_scores_minus_out):
            emit(dataset, score_path, reverse=True)
    log("finished start CNN score pipeline")


if __name__ == "__main__":
    main()
