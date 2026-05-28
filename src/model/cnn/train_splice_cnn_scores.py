from __future__ import annotations

import argparse
import random
from collections import defaultdict
from pathlib import Path

import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

from splice_cnn_network import SpliceCNN, one_hot_encode_windows


def read_fasta(path: Path) -> tuple[str, dict[str, int]]:
    sequence_parts: list[str] = []
    offsets: dict[str, int] = {}
    current_name: str | None = None
    current_length = 0

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    current_length = len("".join(sequence_parts))
                current_name = line[1:].split()[0]
                offsets[current_name] = current_length
                continue
            sequence_parts.append(line.upper())

    return "".join(sequence_parts), offsets


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
    sequence: str,
    donors: set[int],
    acceptors: set[int],
    radius: int,
    negatives_per_positive: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    positives = sorted(donors | acceptors)
    positive_set = set(positives)
    candidate_negatives = [
        pos
        for pos in range(1, len(sequence) - 1)
        if pos not in positive_set
        and (
            sequence[pos : pos + 2] == "GT"
            or sequence[pos - 1 : pos + 1] == "AG"
        )
    ]

    rng = random.Random(7)
    negative_count = min(len(candidate_negatives), len(positives) * negatives_per_positive)
    negatives = rng.sample(candidate_negatives, negative_count)

    windows: list[str] = []
    labels: list[list[float]] = []
    for pos in positives:
        windows.append(window_at(sequence, pos, radius))
        labels.append([1.0 if pos in donors else 0.0, 1.0 if pos in acceptors else 0.0])
    for pos in negatives:
        windows.append(window_at(sequence, pos, radius))
        labels.append([0.0, 0.0])

    return one_hot_encode_windows(windows), torch.tensor(labels, dtype=torch.float32)


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

    model.train()
    for epoch in range(epochs):
        total_loss = 0.0
        for batch_inputs, batch_labels in loader:
            optimizer.zero_grad()
            loss = loss_fn(model(batch_inputs), batch_labels)
            loss.backward()
            optimizer.step()
            total_loss += float(loss.item()) * batch_inputs.size(0)
        print(f"epoch={epoch + 1} loss={total_loss / len(inputs):.4f}")


def write_scores(model: SpliceCNN, sequence: str, radius: int, output_path: Path, batch_size: int) -> None:
    model.eval()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as output, torch.no_grad():
        for start in range(0, len(sequence), batch_size):
            end = min(start + batch_size, len(sequence))
            windows = [window_at(sequence, pos, radius) for pos in range(start, end)]
            logits = model(one_hot_encode_windows(windows))
            for offset, row in enumerate(logits):
                position = start + offset
                output.write(f"{position}\t{float(row[0]):.8f}\t{float(row[1]):.8f}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Train a splice-site CNN and write HMM emission scores.")
    parser.add_argument("--train-fasta", required=True, type=Path)
    parser.add_argument("--train-gff", required=True, type=Path)
    parser.add_argument("--test-fasta", required=True, type=Path)
    parser.add_argument("--model-out", required=True, type=Path)
    parser.add_argument("--train-scores-out", required=True, type=Path)
    parser.add_argument("--test-scores-out", required=True, type=Path)
    parser.add_argument("--radius", default=40, type=int)
    parser.add_argument("--epochs", default=5, type=int)
    parser.add_argument("--batch-size", default=256, type=int)
    parser.add_argument("--negatives-per-positive", default=5, type=int)
    args = parser.parse_args()

    train_sequence, train_offsets = read_fasta(args.train_fasta)
    test_sequence, _ = read_fasta(args.test_fasta)
    model = SpliceCNN(window_size=args.radius * 2 + 1)

    if args.model_out.exists():
        print(f"loading existing model: {args.model_out}")
        model.load_state_dict(torch.load(args.model_out, map_location="cpu"))
    else:
        donors, acceptors = splice_sites_from_gff(
            args.train_gff,
            train_offsets,
            min_intron_bp=20,
            require_3n_cds=True,
        )

        print(f"training donors={len(donors)} acceptors={len(acceptors)}")
        inputs, labels = sample_training_examples(
            train_sequence,
            donors,
            acceptors,
            args.radius,
            args.negatives_per_positive,
        )

        train_model(model, inputs, labels, args.epochs, args.batch_size)
        args.model_out.parent.mkdir(parents=True, exist_ok=True)
        torch.save(model.state_dict(), args.model_out)

    write_scores(model, train_sequence, args.radius, args.train_scores_out, args.batch_size)
    write_scores(model, test_sequence, args.radius, args.test_scores_out, args.batch_size)


if __name__ == "__main__":
    main()
