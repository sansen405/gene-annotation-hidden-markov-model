from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PROFILE = REPO_ROOT / "src/genome_profiles/fission_yeasts.json"
DEFAULT_BINARY = Path("/tmp/train_hmm_matrices")


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def run(command: list[str], description: str) -> None:
    start_time = time.perf_counter()
    log(description)
    log("+ " + " ".join(command))
    subprocess.run(command, cwd=REPO_ROOT, check=True)
    log(f"finished: {description} elapsed={time.perf_counter() - start_time:.1f}s")


def compile_hmm_trainer(binary_path: Path, cxx: str) -> None:
    run(
        [
            cxx,
            "-std=c++17",
            "-Isrc",
            "-I/opt/homebrew/include",
            "src/model/training_pipeline/train_hmm_matrices.cpp",
            "src/model/transition/Transition_Model.cpp",
            "src/model/emission/Emission_Model.cpp",
            "src/model/cnn/Splice_CNN_Scores.cpp",
            "src/genome_profiles/Genome_Profile.cpp",
            "src/parsers/FNA_Parser.cpp",
            "src/parsers/GFF_Parser.cpp",
            "-o",
            str(binary_path),
        ],
        "compiling HMM matrix trainer",
    )


def refresh_cnn_scores(profile_path: Path, python: str) -> None:
    with profile_path.open() as handle:
        profile = json.load(handle)

    dataset = profile["dataset"]
    splice_cnn = profile["splice_cnn"]

    run(
        [
            python,
            "src/model/cnn/train_splice_cnn_scores.py",
            "--train-fasta",
            dataset["train_fasta"],
            "--train-gff",
            dataset["train_gff"],
            "--test-fasta",
            dataset["test_fasta"],
            "--model-out",
            splice_cnn["model"],
            "--train-scores-out",
            splice_cnn["train_scores"],
            "--test-scores-out",
            splice_cnn["test_scores"],
        ],
        "training/loading CNN and refreshing splice score files",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train/cache the CNN-backed HMM transition and emission artifacts."
    )
    parser.add_argument("--profile", default=DEFAULT_PROFILE, type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--binary", default=DEFAULT_BINARY, type=Path)
    parser.add_argument("--python", default="python3")
    parser.add_argument("--cxx", default="clang++")
    parser.add_argument("--skip-cnn", action="store_true")
    parser.add_argument("--skip-compile", action="store_true")
    args = parser.parse_args()

    profile_path = args.profile
    if not profile_path.is_absolute():
        profile_path = REPO_ROOT / profile_path

    log(f"starting cached model pipeline profile={profile_path}")
    if not args.skip_cnn:
        refresh_cnn_scores(profile_path, args.python)
    else:
        log("skipping CNN score refresh")

    if not args.skip_compile:
        compile_hmm_trainer(args.binary, args.cxx)
    else:
        log("skipping HMM trainer compile")

    command = [
        str(args.binary),
        "--profile",
        str(profile_path.relative_to(REPO_ROOT) if profile_path.is_relative_to(REPO_ROOT) else profile_path),
    ]
    if args.out_dir is not None:
        command.extend(["--out-dir", str(args.out_dir)])
    run(command, "training and saving HMM transition/emission artifacts")
    log("finished cached model pipeline")


if __name__ == "__main__":
    main()
