#!/usr/bin/env bash
# Train HMM models for every genome profile in src/genome_profiles/ and save them to models/.
# Run this locally before building the Docker image.
# Requires: genome training data in the paths referenced by each profile JSON.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODELS_DIR="$ROOT_DIR/models"
TRAIN_BIN="$ROOT_DIR/train_model"
SRC_DIRS=(
    "$ROOT_DIR/src/model/Transition_Model.cpp"
    "$ROOT_DIR/src/model/Emission_Model.cpp"
    "$ROOT_DIR/src/model/Model_IO.cpp"
    "$ROOT_DIR/src/genome_profiles/Genome_Profile.cpp"
    "$ROOT_DIR/src/parsers/FNA_Parser.cpp"
    "$ROOT_DIR/src/parsers/GFF_Parser.cpp"
)

# ── Build train_model if needed ───────────────────────────────────────────────
if [ ! -f "$TRAIN_BIN" ]; then
    echo "Building train_model..."
    clang++ -std=c++17 -O2 \
        "-I$ROOT_DIR/src" \
        "-I/opt/homebrew/include" \
        "$ROOT_DIR/train_model.cpp" \
        "${SRC_DIRS[@]}" \
        -o "$TRAIN_BIN"
    echo "Built: $TRAIN_BIN"
fi

mkdir -p "$MODELS_DIR"

# ── Train one model per genome profile ───────────────────────────────────────
shopt -s nullglob
profiles=("$ROOT_DIR/src/genome_profiles/"*.json)

if [ ${#profiles[@]} -eq 0 ]; then
    echo "No genome profiles found in src/genome_profiles/." >&2
    exit 1
fi

for profile in "${profiles[@]}"; do
    name="$(basename "$profile" .json)"
    output="$MODELS_DIR/${name}_model.json"
    echo "─── Training $name ───"
    "$TRAIN_BIN" --profile "$profile" --output "$output"
    echo "Saved: $output"
done

echo ""
echo "All models trained. Ready to build Docker image:"
echo "  docker compose build"
echo "  docker compose up"
