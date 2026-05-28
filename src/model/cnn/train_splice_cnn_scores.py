from __future__ import annotations

import argparse
import json
import random
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import torch
    from splice_cnn_network import SpliceCNN


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


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def format_count(value: int) -> str:
    return f"{value:,}"


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


def splice_sites_from_gff(
    gff_path: Path,
    offsets: dict[str, int],
    min_intron_bp: int,
    require_3n_cds: bool,
) -> tuple[set[int], set[int]]:
    cds_by_parent: dict[str, list[tuple[int, int]]] = defaultdict(list)
    supported_parent: dict[str, bool] = defaultdict(lambda: True)

    log(f"loading splice sites from GFF: {gff_path}")
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
            if strand != "+":
                supported_parent[key] = False

    donors: set[int] = set()
    acceptors: set[int] = set()
    for key, fragments in cds_by_parent.items():
        fragments.sort()
        total_cds = sum(end - start + 1 for start, end in fragments)
        if not supported_parent[key]:
            continue
        if require_3n_cds and total_cds % 3 != 0:
            continue

        valid = True
        for left, right in zip(fragments, fragments[1:]):
            intron_bp = right[0] - left[1] - 1
            if intron_bp < min_intron_bp:
                valid = False
                break
        if not valid:
            continue

        for left, right in zip(fragments, fragments[1:]):
            donors.add(left[1] + 1)
            acceptors.add(right[0] - 1)

    log(
        f"loaded splice sites: donors={format_count(len(donors))} "
        f"acceptors={format_count(len(acceptors))} elapsed={time.perf_counter() - start_time:.1f}s"
    )
    return donors, acceptors


def window_at(sequence: str, center: int, radius: int) -> str:
    chars: list[str] = []
    for pos in range(center - radius, center + radius + 1):
        if pos < 0 or pos >= len(sequence):
            chars.append("N")
        else:
            chars.append(sequence[pos])
    return "".join(chars)


def sample_training_examples(
    dataset: FastaDataset,
    donors: set[int],
    acceptors: set[int],
    radius: int,
    negatives_per_positive: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    log(f"sampling examples from {dataset.path}")
    start_time = time.perf_counter()
    positives = sorted(donors | acceptors)
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
        for local_pos in range(1, len(record.sequence) - 1):
            global_pos = record.offset + local_pos
            if global_pos in positive_set:
                continue
            if (
                record.sequence[local_pos : local_pos + 2] == "GT"
                or record.sequence[local_pos - 1 : local_pos + 1] == "AG"
            ):
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
            windows.append(window_at(record.sequence, local_pos, radius))
            labels.append([1.0 if pos in donors else 0.0, 1.0 if pos in acceptors else 0.0])
    for record, local_pos, _ in negatives:
        windows.append(window_at(record.sequence, local_pos, radius))
        labels.append([0.0, 0.0])

    if not windows:
        window_size = radius * 2 + 1
        return torch.empty((0, 4, window_size), dtype=torch.float32), torch.empty((0, 2), dtype=torch.float32)

    log(f"encoding windows: count={format_count(len(windows))} width={radius * 2 + 1}")
    inputs = one_hot_encode_windows(windows)
    label_tensor = torch.tensor(labels, dtype=torch.float32)
    log(f"finished examples from {dataset.path} elapsed={time.perf_counter() - start_time:.1f}s")
    return inputs, label_tensor


def train_model(
    model: SpliceCNN,
    inputs: torch.Tensor,
    labels: torch.Tensor,
    epochs: int,
    batch_size: int,
) -> None:
    loader = DataLoader(TensorDataset(inputs, labels), batch_size=batch_size, shuffle=True)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.BCEWithLogitsLoss()
    total_batches = len(loader)

    log(
        f"starting CNN training: examples={format_count(len(inputs))} "
        f"epochs={epochs} batch_size={batch_size} batches_per_epoch={format_count(total_batches)}"
    )
    model.train()
    for epoch in range(epochs):
        epoch_start = time.perf_counter()
        total_loss = 0.0
        progress_interval = max(1, total_batches // 10)
        for batch_index, (batch_inputs, batch_labels) in enumerate(loader, start=1):
            optimizer.zero_grad()
            loss = loss_fn(model(batch_inputs), batch_labels)
            loss.backward()
            optimizer.step()
            total_loss += float(loss.item()) * batch_inputs.size(0)
            if batch_index == total_batches or batch_index % progress_interval == 0:
                log(f"epoch {epoch + 1}/{epochs}: batch {batch_index}/{total_batches}")
        log(
            f"finished epoch {epoch + 1}/{epochs}: "
            f"loss={total_loss / len(inputs):.4f} elapsed={time.perf_counter() - epoch_start:.1f}s"
        )


def write_scores(model: SpliceCNN, dataset: FastaDataset, radius: int, output_path: Path, batch_size: int) -> None:
    model.eval()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    log(
        f"writing scores: input={dataset.path} output={output_path} "
        f"bases={format_count(dataset.base_count)}"
    )
    start_time = time.perf_counter()
    written_bases = 0
    next_progress = 10
    with output_path.open("w") as output, torch.no_grad():
        for record in dataset.records:
            for start in range(0, len(record.sequence), batch_size):
                end = min(start + batch_size, len(record.sequence))
                windows = [window_at(record.sequence, pos, radius) for pos in range(start, end)]
                logits = model(one_hot_encode_windows(windows))
                for offset, row in enumerate(logits):
                    position = record.offset + start + offset
                    output.write(f"{position}\t{float(row[0]):.8f}\t{float(row[1]):.8f}\n")
                written_bases += end - start
                progress = int((written_bases / max(1, dataset.base_count)) * 100)
                if progress >= next_progress:
                    log(f"score writing progress for {output_path}: {min(progress, 100)}%")
                    next_progress += 10
    log(f"finished scores: {output_path} elapsed={time.perf_counter() - start_time:.1f}s")


def candidate_splice_positions(sequence: str) -> list[int]:
    positions: set[int] = set()
    for pos in range(1, len(sequence) - 1):
        if sequence[pos : pos + 2] == "GT" or sequence[pos - 1 : pos + 1] == "AG":
            positions.add(pos)
    return sorted(positions)


def write_sparse_scores(model: SpliceCNN, dataset: FastaDataset, radius: int, output_path: Path, batch_size: int) -> None:
    model.eval()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    candidate_counts = [len(candidate_splice_positions(record.sequence)) for record in dataset.records]
    total_candidates = sum(candidate_counts)
    log(
        f"writing sparse scores: input={dataset.path} output={output_path} "
        f"bases={format_count(dataset.base_count)} candidates={format_count(total_candidates)}"
    )

    start_time = time.perf_counter()
    written_candidates = 0
    next_progress = 10
    with output_path.open("w") as output, torch.no_grad():
        for record, candidate_count in zip(dataset.records, candidate_counts):
            positions = candidate_splice_positions(record.sequence)
            if len(positions) != candidate_count:
                raise RuntimeError("Candidate position count changed during sparse scoring.")
            for start in range(0, len(positions), batch_size):
                batch_positions = positions[start : start + batch_size]
                windows = [window_at(record.sequence, pos, radius) for pos in batch_positions]
                logits = model(one_hot_encode_windows(windows))
                for local_pos, row in zip(batch_positions, logits):
                    position = record.offset + local_pos
                    output.write(f"{position}\t{float(row[0]):.8f}\t{float(row[1]):.8f}\n")
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
    raise ValueError(f"Profile splice_cnn.{label} must be a string or list of strings.")


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
    splice_cnn = profile.get("splice_cnn")
    if not isinstance(splice_cnn, dict):
        raise ValueError("Profile must define splice_cnn paths for CNN score generation.")

    args.train_fasta = args.train_fasta or path_list_from_species(species, "train_fasta")
    args.train_gff = args.train_gff or path_list_from_species(species, "train_gff")
    args.test_fasta = args.test_fasta or path_list_from_species(species, "test_fasta")
    if args.model_out is None:
        model_path = splice_cnn.get("model")
        if not isinstance(model_path, str) or not model_path:
            raise ValueError("Profile splice_cnn.model must be a non-empty string.")
        args.model_out = Path(model_path)
    args.train_scores_out = args.train_scores_out or list_from_json(splice_cnn.get("train_scores"), "train_scores")
    args.test_scores_out = args.test_scores_out or list_from_json(splice_cnn.get("test_scores"), "test_scores")

    filters = profile.get("filters")
    if isinstance(filters, dict):
        if args.min_intron_bp is None:
            args.min_intron_bp = int(filters.get("min_intron_bp", 20))
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
    parser = argparse.ArgumentParser(description="Train a splice-site CNN and write HMM emission scores.")
    parser.add_argument("--profile", type=Path, help="Genome profile JSON that defines dataset and splice_cnn paths.")
    parser.add_argument("--train-fasta", type=Path, nargs="+")
    parser.add_argument("--train-gff", type=Path, nargs="+")
    parser.add_argument("--test-fasta", type=Path, nargs="+")
    parser.add_argument("--model-out", type=Path)
    parser.add_argument("--train-scores-out", type=Path, nargs="+")
    parser.add_argument("--test-scores-out", type=Path, nargs="+")
    parser.add_argument("--radius", default=40, type=int)
    parser.add_argument("--epochs", default=5, type=int)
    parser.add_argument("--batch-size", default=256, type=int)
    parser.add_argument("--score-batch-size", default=8192, type=int)
    parser.add_argument("--sparse-scores", action="store_true")
    parser.add_argument("--negatives-per-positive", default=5, type=int)
    parser.add_argument("--min-intron-bp", type=int)
    parser.add_argument("--require-3n-cds", dest="require_3n_cds", action="store_true")
    parser.add_argument("--allow-non-3n-cds", dest="require_3n_cds", action="store_false")
    parser.set_defaults(require_3n_cds=None)
    args = parser.parse_args()

    apply_profile_defaults(args)
    if args.min_intron_bp is None:
        args.min_intron_bp = 20
    if args.require_3n_cds is None:
        args.require_3n_cds = True
    validate_required_args(args)

    global torch, nn, DataLoader, TensorDataset, SpliceCNN, one_hot_encode_windows
    import torch
    from splice_cnn_network import SpliceCNN, one_hot_encode_windows
    from torch import nn
    from torch.utils.data import DataLoader, TensorDataset

    validate_matching_counts("training", args.train_fasta, args.train_gff)
    validate_output_counts("training", args.train_fasta, args.train_scores_out)
    validate_output_counts("test", args.test_fasta, args.test_scores_out)

    log("starting splice CNN score pipeline")
    log(f"model checkpoint: {args.model_out}")
    train_datasets = [read_fasta(path) for path in args.train_fasta]
    test_datasets = [read_fasta(path) for path in args.test_fasta]
    model = SpliceCNN(window_size=args.radius * 2 + 1)

    if args.model_out.exists():
        log(f"loading existing CNN checkpoint: {args.model_out}")
        model.load_state_dict(torch.load(args.model_out, map_location="cpu"))
    else:
        input_batches: list[torch.Tensor] = []
        label_batches: list[torch.Tensor] = []
        donor_total = 0
        acceptor_total = 0

        for dataset, gff_path in zip(train_datasets, args.train_gff):
            donors, acceptors = splice_sites_from_gff(
                gff_path,
                dataset.offsets,
                min_intron_bp=args.min_intron_bp,
                require_3n_cds=args.require_3n_cds,
            )
            donor_total += len(donors)
            acceptor_total += len(acceptors)
            inputs, labels = sample_training_examples(
                dataset,
                donors,
                acceptors,
                args.radius,
                args.negatives_per_positive,
            )
            if len(labels) > 0:
                input_batches.append(inputs)
                label_batches.append(labels)

        if not input_batches:
            raise ValueError("No splice-site examples were found in the training inputs.")

        log(f"training totals: donors={format_count(donor_total)} acceptors={format_count(acceptor_total)}")
        log("combining training tensors")
        inputs = torch.cat(input_batches)
        labels = torch.cat(label_batches)

        train_model(model, inputs, labels, args.epochs, args.batch_size)
        args.model_out.parent.mkdir(parents=True, exist_ok=True)
        log(f"saving CNN checkpoint: {args.model_out}")
        torch.save(model.state_dict(), args.model_out)

    for dataset, score_path in zip(train_datasets, args.train_scores_out):
        if args.sparse_scores:
            write_sparse_scores(model, dataset, args.radius, score_path, args.score_batch_size)
        else:
            write_scores(model, dataset, args.radius, score_path, args.score_batch_size)
    for dataset, score_path in zip(test_datasets, args.test_scores_out):
        if args.sparse_scores:
            write_sparse_scores(model, dataset, args.radius, score_path, args.score_batch_size)
        else:
            write_scores(model, dataset, args.radius, score_path, args.score_batch_size)
    log("finished splice CNN score pipeline")


if __name__ == "__main__":
    main()
