# Gene_Analysis_HMM

A C++ Hidden Markov Model for ab-initio gene prediction in eukaryotic genomes.
The model uses a 21-state topology (intergenic, start/stop codons, frame-tracked
exons, donor/acceptor splice sites, and frame-tracked introns) and is trained
from a genome FASTA plus a GFF annotation. Per-genome behavior (training data,
filters, emission models, smoothing) is driven by a JSON profile in
`src/genome_profiles/`.

## Repository Layout

```text
Gene_Analysis_HMM/
├── README.md
├── .vscode/
│   └── settings.json
├── src/
│   ├── topology/
│   │   └── Topology.hpp           # State enum, Nucleotide enum, allowed transitions
│   ├── parsers/
│   │   ├── FNA_Parser.hpp         # FASTA -> nucleotide sequence + chromosome ranges
│   │   ├── FNA_Parser.cpp
│   │   ├── GFF_Parser.hpp         # GFF  -> per-base region/state labels
│   │   └── GFF_Parser.cpp
│   ├── model/
│   │   ├── transition/
│   │   │   ├── Transition_Model.hpp
│   │   │   └── Transition_Model.cpp
│   │   ├── emission/
│   │   │   ├── Emission_Model.hpp
│   │   │   └── Emission_Model.cpp
│   │   ├── training_pipeline/
│   │   │   ├── train_cached_model.py
│   │   │   ├── train_hmm_matrices.cpp
│   │   │   └── trained_models/     # local cached HMM artifacts, gitignored
│   │   └── cnn/
│   │       ├── Splice_CNN_Scores.hpp/.cpp
│   │       ├── splice_cnn_network.py
│   │       ├── train_splice_cnn_scores.py
│   │       └── trained_models/     # local .pt checkpoints, gitignored
│   ├── decoding/
│   │   ├── Viterbi.hpp
│   │   ├── Viterbi.cpp
│   │   ├── Forward_Backward.hpp
│   │   └── Forward_Backward.cpp
│   ├── tools/
│   │   ├── predict_fna.cpp        # JSON prediction CLI for uploaded .fna input
│   │   └── split_genome_data.cpp  # writes train/test FASTA+GFF from a profile
│   └── genome_profiles/
│       └── fission_yeasts.json
├── validation/
│   ├── full_genome_validation.cpp # train/test holdout validation runner
│   └── diagnostics/               # extra boundary and intron diagnostics
└── genome_data/
    └── fission_yeasts/
        ├── full/                  # unsplit source assembly + annotation
        ├── train/                 # training chromosomes only
        └── test/                  # held-out test chromosomes only
```

## Components

### `src/topology/`
Defines the shared vocabulary used across the project:

- `Nucleotide` and `State` enums.
- `NUM_STATES`, `NUM_NUCLEOTIDES`, `MEMORY_WINDOW` constants.
- `Transitions`: the allowed state-to-state edges in the HMM topology.

### `src/parsers/`
- **`FNA_Parser`** — reads a FASTA file into a `vector<Nucleotide>` and reports
  per-chromosome offsets and ranges so downstream code can treat the genome as
  one concatenated sequence without crossing chromosome boundaries.
- **`GFF_Parser`** — converts a GFF annotation into a per-base label vector
  (regions, then HMM states) aligned to the FASTA. CDS groups that do not pass
  the active profile filters are marked as ignored annotation before state
  labeling, so unsupported genes do not contribute training or evaluation
  signal.

### `src/model/transition/`
- **`Transition_Model`** — given the per-base state vector and chromosome
  ranges, counts bigram transitions, computes row sums, and produces a
  log-probability transition matrix with additive (`alpha`) smoothing.

### `src/model/emission/`
- **`Emission_Model`** — given the per-base state and nucleotide vectors,
  counts HMM emissions and combines them with CNN splice-site scores:
  - **Markov order-1** (`INTERGENIC`, `INTRON`) — 4×4 context→emission table.
  - **Frame-specific Markov order-5** (`EXON_FRAME_1/2/3`) — one 1024×4 table
    per coding frame, so the model can use codon periodicity.
  - **CNN splice scores** (`DONOR`, `ACCEPTOR`) — donor/acceptor states require
    per-position log-odds produced by the PyTorch splice CNN.
  - **Signal log-odds table** (`START_CODON`) — the start window is scored
    against matched non-start `ATG` background.
  - **Deterministic** (`START_CODON`, `STOP_CODON`) — returns 0.0 or −∞ based
    on the expected nucleotide context. Stop codons are restricted to `TAA`,
    `TAG`, and `TGA`.

### `src/model/cnn/`
- **`splice_cnn_network.py`** — PyTorch 1D CNN that scores donor and acceptor
  splice-site potential from one-hot DNA windows.
- **`train_splice_cnn_scores.py`** — trains the CNN when no cached checkpoint is
  present, saves the checkpoint under `trained_models/`, and writes train/test
  per-position score TSV files.
- **`Splice_CNN_Scores`** — C++ bridge that loads those score TSV files and
  exposes donor/acceptor log-odds to `Emission_Model`.

### `src/model/training_pipeline/`
- **`train_cached_model.py`** — one-command local cache refresh. It trains or
  reloads the CNN checkpoint, regenerates CNN score TSVs, compiles the HMM
  matrix trainer, and writes cached transition/emission artifacts.
- **`train_hmm_matrices.cpp`** — trains the HMM transition matrix and reusable
  emission tables from the profile train FASTA/GFF and saves them as JSON under
  `trained_models/<profile>/`.

Emission summary:

| State family | Emission type | Context/window | What is scored |
| --- | --- | --- | --- |
| `INTERGENIC` | Markov order-1 | Previous base | `log P(current base \| previous base)` learned from intergenic sequence. |
| `INTRON_1/2/3` | Markov order-1 | Previous base | `log P(current base \| previous base)` learned from intron sequence. |
| `EXON_FRAME_1/2/3` | Frame-specific Markov order-5 | Previous 5 bases | `log P(current base \| 5-base context)` with a separate table for each coding frame. |
| `DONOR_1/2/3` | CNN log-odds | Per-base score TSV | Requires canonical `GT`, then uses the loaded CNN donor score for that position. |
| `ACCEPTOR_1/2/3` | CNN log-odds | Per-base score TSV | Requires canonical `AG`, then uses the loaded CNN acceptor score for that position. |
| `START_CODON_1` | Deterministic motif + log-odds PSSM | Default `window_left=6`, `window_right=9` | Requires full `ATG`, then scores the surrounding start window against non-start `ATG` background. |
| `START_CODON_2/3` | Deterministic | Current base | Requires the second and third start-codon bases, `T` then `G`. |
| `STOP_CODON_1/2/3` | Deterministic motif | Three-base stop context | Requires one of `TAA`, `TAG`, or `TGA`. |

### `src/decoding/`
`Viterbi` decodes the best state path. The validation runner uses the extended
decoder overload with an intron body length cap and a gene-start penalty.
`Forward_Backward` computes posterior state probabilities and returns per-base
confidence for a predicted path.

### `src/genome_profiles/`
JSON files that fully describe a training/decoding run: input paths, held-out
test chromosomes, gene-quality filters, per-state emission model choices, and
smoothing hyperparameters.

Example: [`src/genome_profiles/fission_yeasts.json`](src/genome_profiles/fission_yeasts.json).

### `genome_data/`
The checked-in fission yeast profile keeps each species in its own folder under
`genome_data/fission_yeasts/`:

```text
genome_data/fission_yeasts/
├── s_pombe/
│   ├── full/
│   ├── train/
│   └── test/
├── s_japonicus/
├── s_octosporus/
└── s_cryophilus/
```

The profile trains one shared model across four Schizosaccharomyces assemblies,
but the FASTA/GFF files are not merged on disk:

- `S. pombe` 972h- (`GCF_000002945.1`)
- `S. japonicus` yFS275 (`GCA_000149845.2`)
- `S. octosporus` yFS286 (`GCA_000150505.2`)
- `S. cryophilus` OY26 (`GCA_000004155.2`)

The current split has 66 training chromosomes and 4 held-out evaluation
chromosomes, one held-out chromosome/scaffold per species. The training pipeline
and validation runner load all four species' `train/` files and combine them in
memory for parameter fitting.

Regenerate train/test files after changing source data or holdout chromosomes:

```sh
clang++ -std=c++17 -Isrc -I/opt/homebrew/include \
  src/tools/split_genome_data.cpp \
  src/genome_profiles/Genome_Profile.cpp \
  -o /tmp/split_genome_data

/tmp/split_genome_data src/genome_profiles/fission_yeasts.json
```

## Genome Profile Schema

Each profile in `src/genome_profiles/` follows this shape:

```json
{
  "name": "fission_yeasts",
  "description": "...",

  "dataset": {
    "species": [
      {
        "name": "s_pombe",
        "source_fasta": "genome_data/fission_yeasts/s_pombe/full/s_pombe_full.fna",
        "source_gff":   "genome_data/fission_yeasts/s_pombe/full/s_pombe_full.gff",
        "train_fasta":  "genome_data/fission_yeasts/s_pombe/train/s_pombe_train.fna",
        "train_gff":    "genome_data/fission_yeasts/s_pombe/train/s_pombe_train.gff",
        "test_fasta":   "genome_data/fission_yeasts/s_pombe/test/s_pombe_test.fna",
        "test_gff":     "genome_data/fission_yeasts/s_pombe/test/s_pombe_test.gff",
        "test_chromosomes": ["NC_003424.3"],
        "excluded_chromosomes": ["NC_001326.1"]
      }
    ]
  },

  "splice_cnn": {
    "model": "src/model/cnn/trained_models/fission_yeasts_splice_cnn.pt",
    "train_scores": [
      "genome_data/fission_yeasts/s_pombe/train/s_pombe_splice_cnn_scores.tsv"
    ],
    "test_scores": [
      "genome_data/fission_yeasts/s_pombe/test/s_pombe_splice_cnn_scores.tsv"
    ]
  },

  "filters": {
    "min_first_cds_bp":     3,
    "min_last_cds_bp":      3,
    "min_intron_bp":        20,
    "require_3n_cds":       true,
    "include_minus_strand": false
  },

  "emissions": {
    "INTERGENIC":  { "type": "markov", "order": 1 },
    "INTRON":      { "type": "markov", "order": 1 },
    "EXON_FRAME":  { "type": "markov", "order": 5, "frame_tied": true },
    "DONOR":       { "type": "pssm", "window_left": 3,  "window_right": 6 },
    "ACCEPTOR":    { "type": "pssm", "window_left": 15, "window_right": 3 },
    "START_CODON": { "type": "deterministic" },
    "STOP_CODON":  { "type": "deterministic" }
  },

  "smoothing": {
    "transition_alpha": 0.02,
    "emission_alpha":   0.1
  }
}
```

## Build

There is no Makefile or CMake project yet. Build directly with `clang++`.
The project expects `nlohmann/json.hpp` to be available from Homebrew at
`/opt/homebrew/include`.

Build and run unit tests:

```sh
clang++ -std=c++17 -Isrc -I/opt/homebrew/include \
  src/main.cpp \
  src/decoding/Viterbi.cpp \
  src/decoding/Forward_Backward.cpp \
  src/model/transition/Transition_Model.cpp \
  src/genome_profiles/Genome_Profile.cpp \
  src/parsers/FNA_Parser.cpp \
  src/parsers/GFF_Parser.cpp \
  src/model/emission/Emission_Model.cpp \
  src/model/cnn/Splice_CNN_Scores.cpp \
  -o /tmp/gene_hmm_tests

/tmp/gene_hmm_tests
```

Build and run the default fission yeast holdout validation:

```sh
clang++ -std=c++17 -Isrc -I/opt/homebrew/include \
  validation/full_genome_validation.cpp \
  src/decoding/Viterbi.cpp \
  src/model/transition/Transition_Model.cpp \
  src/genome_profiles/Genome_Profile.cpp \
  src/parsers/FNA_Parser.cpp \
  src/parsers/GFF_Parser.cpp \
  src/model/emission/Emission_Model.cpp \
  src/model/cnn/Splice_CNN_Scores.cpp \
  -o /tmp/full_genome_validation

/tmp/full_genome_validation
```

CNN splice logits can be calibrated before Viterbi with profile fields
(`donor_scale`, `donor_bias`, `acceptor_scale`, `acceptor_bias` under
`splice_cnn`) or with validation CLI flags:

```sh
/tmp/full_genome_validation --profile src/genome_profiles/fission_yeasts.json \
  --tune-cnn-calibration
```

The tuner tests a compact set of donor/acceptor scale and bias candidates and
selects the run that maximizes the average of intron F1, donor boundary F1, and
acceptor boundary F1.

The validation runner expects CNN splice-score files at the paths configured in
`src/genome_profiles/fission_yeasts.json`. Generate or refresh them with:

```sh
python3 src/model/cnn/train_splice_cnn_scores.py \
  --train-fasta \
    genome_data/fission_yeasts/s_pombe/train/s_pombe_train.fna \
    genome_data/fission_yeasts/s_japonicus/train/s_japonicus_train.fna \
    genome_data/fission_yeasts/s_octosporus/train/s_octosporus_train.fna \
    genome_data/fission_yeasts/s_cryophilus/train/s_cryophilus_train.fna \
  --train-gff \
    genome_data/fission_yeasts/s_pombe/train/s_pombe_train.gff \
    genome_data/fission_yeasts/s_japonicus/train/s_japonicus_train.gff \
    genome_data/fission_yeasts/s_octosporus/train/s_octosporus_train.gff \
    genome_data/fission_yeasts/s_cryophilus/train/s_cryophilus_train.gff \
  --test-fasta \
    genome_data/fission_yeasts/s_pombe/test/s_pombe_test.fna \
    genome_data/fission_yeasts/s_japonicus/test/s_japonicus_test.fna \
    genome_data/fission_yeasts/s_octosporus/test/s_octosporus_test.fna \
    genome_data/fission_yeasts/s_cryophilus/test/s_cryophilus_test.fna \
  --train-scores-out \
    genome_data/fission_yeasts/s_pombe/train/s_pombe_splice_cnn_scores.tsv \
    genome_data/fission_yeasts/s_japonicus/train/s_japonicus_splice_cnn_scores.tsv \
    genome_data/fission_yeasts/s_octosporus/train/s_octosporus_splice_cnn_scores.tsv \
    genome_data/fission_yeasts/s_cryophilus/train/s_cryophilus_splice_cnn_scores.tsv \
  --test-scores-out \
    genome_data/fission_yeasts/s_pombe/test/s_pombe_splice_cnn_scores.tsv \
    genome_data/fission_yeasts/s_japonicus/test/s_japonicus_splice_cnn_scores.tsv \
    genome_data/fission_yeasts/s_octosporus/test/s_octosporus_splice_cnn_scores.tsv \
    genome_data/fission_yeasts/s_cryophilus/test/s_cryophilus_splice_cnn_scores.tsv
```

If the `.pt` checkpoint already exists, the script loads it and only regenerates
the score TSVs. The checkpoint and generated score files are local artifacts and
are ignored by git.

To refresh the whole cached model pipeline in one step:

```sh
python3 src/model/training_pipeline/train_cached_model.py \
  --profile src/genome_profiles/fission_yeasts.json
```

This writes local artifacts under
`src/model/training_pipeline/trained_models/fission_yeasts/`:
`transition_matrix.json`, `emission_matrix.json`, and `metadata.json`. The
emission artifact stores the trained HMM emission tables and points donor/
acceptor emissions at the CNN checkpoint and per-position score files from the
profile.

By default, CNN score generation writes sparse score TSVs for canonical splice
candidate positions only (`GT` donor candidates and `AG` acceptor candidates).
This is much faster than scoring every base and is compatible with the current
decoder because donor/acceptor emissions already reject noncanonical motifs. To
force dense every-base score TSVs, add `--dense-cnn-scores`.

Build the JSON decoder:

```sh
clang++ -std=c++17 -Isrc -I/opt/homebrew/include \
  src/tools/predict_fna.cpp \
  src/decoding/Viterbi.cpp \
  src/decoding/Forward_Backward.cpp \
  src/model/transition/Transition_Model.cpp \
  src/genome_profiles/Genome_Profile.cpp \
  src/parsers/FNA_Parser.cpp \
  src/parsers/GFF_Parser.cpp \
  src/model/emission/Emission_Model.cpp \
  src/model/cnn/Splice_CNN_Scores.cpp \
  -o /tmp/hmm_predict_fna
```

To decode one specific chromosome, first write that chromosome to its own FASTA
and generate matching CNN scores from the cached CNN checkpoint:

```sh
python3 - <<'PY'
from pathlib import Path

chrom = "NC_003424.3"
source = Path("genome_data/fission_yeasts/s_pombe/test/s_pombe_test.fna")
out = Path(f"/tmp/{chrom}.fna")
writing = False

with source.open() as input_file, out.open("w") as output_file:
    for line in input_file:
        if line.startswith(">"):
            writing = line[1:].split()[0] == chrom
        if writing:
            output_file.write(line)

print(out)
PY

python3 src/model/cnn/train_splice_cnn_scores.py \
  --train-fasta \
    genome_data/fission_yeasts/s_pombe/train/s_pombe_train.fna \
    genome_data/fission_yeasts/s_japonicus/train/s_japonicus_train.fna \
    genome_data/fission_yeasts/s_octosporus/train/s_octosporus_train.fna \
    genome_data/fission_yeasts/s_cryophilus/train/s_cryophilus_train.fna \
  --train-gff \
    genome_data/fission_yeasts/s_pombe/train/s_pombe_train.gff \
    genome_data/fission_yeasts/s_japonicus/train/s_japonicus_train.gff \
    genome_data/fission_yeasts/s_octosporus/train/s_octosporus_train.gff \
    genome_data/fission_yeasts/s_cryophilus/train/s_cryophilus_train.gff \
  --test-fasta /tmp/NC_003424.3.fna \
  --train-scores-out \
    /tmp/s_pombe_train_splice_cnn_scores.tsv \
    /tmp/s_japonicus_train_splice_cnn_scores.tsv \
    /tmp/s_octosporus_train_splice_cnn_scores.tsv \
    /tmp/s_cryophilus_train_splice_cnn_scores.tsv \
  --test-scores-out /tmp/NC_003424.3_splice_cnn_scores.tsv
```

Then run the decoder on that chromosome:

```sh
/tmp/hmm_predict_fna \
  --profile src/genome_profiles/fission_yeasts.json \
  --fna /tmp/NC_003424.3.fna \
  --splice-cnn-scores /tmp/NC_003424.3_splice_cnn_scores.tsv \
  > /tmp/NC_003424.3_predictions.json
```

Run validation with an explicit profile:

```sh
/tmp/full_genome_validation --profile src/genome_profiles/fission_yeasts.json
```

## Local Frontend

The `frontend/` directory contains a local React UI plus a local Node API. The
API accepts FASTA/FNA uploads, runs the C++ HMM predictor on the same machine,
and keeps prediction runs in memory for the current browser session only. There
is no cloud storage or external service integration.

Frontend features:

- **Multi-file FASTA/FNA upload** with per-file selection, scaffold counts, base
  counts, and local run timing/status.
- **Saved local runs** in the sidebar for switching between completed analyses
  during the current browser session.
- **Prediction summary cards** for total bases, scaffolds, predicted genes,
  coding segments, and introns.
- **Genome track view** with scaffold selection, coordinate range controls, and
  draggable region selection for zooming into predicted genes.
- **Per-base confidence display** from the Forward-Backward posterior output,
  with limits for large regions.
- **Searchable and paginated gene table** with scaffold/range filters and row
  selection that updates the genome view.
- **Selected prediction details** including exon/intron counts and copyable
  predicted nucleotide sequence.
- **Export endpoints and buttons** for GFF3, CSV, BED, and FASTA output.

```sh
cd frontend
npm install
npm run dev:all
```

Open the Vite URL printed by the terminal, usually
`http://localhost:5173/`. The local API listens on `http://localhost:5174/`.

The API builds and runs `frontend/local_data/bin/hmm_predict_fna` from
`src/tools/predict_fna.cpp`. The predictor requires `--splice-cnn-scores` when
decoding an input FASTA because donor/acceptor emissions no longer fall back to
PSSM scores.

## Model Changes (3.0): Calibrated CNN Splice Model (No Tuning Required)

Version 3 is the current model. The splice CNN now **surpasses** the PSSM splice
model it replaced, and the donor/acceptor scores are self-calibrating, so the HMM
decodes at `scale=1, bias=0` with no `--tune-cnn-calibration` grid search.

Held-out fission-yeast validation at defaults (all four test chromosomes,
~11.7M bases):

| Metric | V2.4 (tuned) | V3 (defaults) |
| --- | --- | --- |
| Intron F1 | ~0.61–0.73 | **0.80** (P 0.79 / R 0.82) |
| Coding F1 | ~0.958 | **0.968** |
| Donor recall | ~0.43 | **~0.66** |
| Acceptor recall | ~0.48 | **~0.67** |

What changed, and how the model is trained differently from earlier versions:

- **Position-aware CNN, not a global-pooled one.** The backbone now ends in
  `AdaptiveMaxPool1d(8)` (eight positional bins) instead of `AdaptiveMaxPool1d(1)`,
  so each head keeps coarse information about *where* a motif sits relative to the
  splice dinucleotide. The earlier global max-pool discarded position, which is the
  single most informative cue for splice sites and the main reason the V2 CNN
  trailed the positionally-anchored PSSM.
- **Per-head asymmetric windows.** Inside `forward()`, the donor head crops a window
  biased downstream of the GT (5' consensus) and the acceptor head crops a wide
  upstream window (branch point + polypyrimidine tract before the AG). The input
  windows are still 121 bp, so score-file generation and the TSV format are
  unchanged.
- **Both strands are trained on.** `splice_sites_from_gff` now keeps minus-strand
  transcripts, tags each site with its strand, and `sample_training_examples`
  reverse-complements minus-strand windows so every true site is presented in
  canonical 5'->3' (GT/AG) orientation. This roughly **doubles** the positive set
  (~6.6k -> ~13k donors/acceptors) without changing the plus-only decoder. Negatives
  stay plus-orientation, so there is no strand leakage.
- **Honest validation split + a real metric.** Training holds out a *seeded,
  shuffled* 10% split (the previous tail split put only easy negatives from the last
  species in validation) and logs per-epoch donor/acceptor **average precision
  (AUPRC)** — a threshold-free measure of splice discrimination that is independent
  of HMM calibration.
- **Calibrated log-odds replace the tuning step.** Training uses plain
  `BCEWithLogits` instead of focal loss (a proper scoring rule, so logits are
  calibrated), then fits a single **temperature** on the held-out split and subtracts
  the measured **training-prior logit** per head. The resulting transform
  (`logit / temperature - prior_logit`) is applied when writing the score TSVs and is
  saved inside the checkpoint. Because the training prior is *measured and removed*,
  the negative:positive sampling ratio (1:5) no longer needs to match the genome's
  (~1:450); the HMM transitions supply the true genomic prior. Typical fitted values:
  `temperature ~ 1.37`, `donor/acceptor prior_logit ~ -2.40`.
- **`--tune-cnn-calibration` is no longer needed.** It still exists for diagnostics,
  but V3 validation is run without it; the calibrated scores already decode at
  `scale=1, bias=0`. This supersedes the V2.4 conclusion that calibration tuning was
  required and that the CNN could not beat the PSSM.
- **Checkpoint format changed.** The `.pt` now stores
  `{"state_dict": ..., "calibration": {temperature, donor_prior_logit,
  acceptor_prior_logit}}` instead of a bare state dict. Old V2 checkpoints are
  incompatible (both the architecture and the format changed) and must be retrained;
  delete the cached `.pt` before refreshing the pipeline.

Train/refresh exactly as before — the one-command pipeline is unchanged:

```sh
rm -f src/model/cnn/trained_models/fission_yeasts_splice_cnn.pt
python3 src/model/training_pipeline/train_cached_model.py \
  --profile src/genome_profiles/fission_yeasts.json
```

Then run validation **without** the calibration flag:

```sh
/tmp/full_genome_validation --profile src/genome_profiles/fission_yeasts.json
```

## Model Changes (2.3): Single CNN Calibration Layer

- **The CNN scorer writes raw logits.**
  `train_splice_cnn_scores.py` no longer applies an analytic log-odds prior
  shift to the donor/acceptor scores. The score TSV is now the untransformed
  CNN output, so there is one well-defined source of truth for the model's
  splice signal.
- **All splice calibration lives in one place.**
  Donor/acceptor emissions are calibrated solely by the HMM's
  `donor_scale`/`donor_bias`/`acceptor_scale`/`acceptor_bias` (profile fields or
  validation CLI flags). The previous setup stacked a Python prior shift on top
  of a separate C++ bias, and stale carried-over biases caused the decoder to
  over-predict introns (intron precision collapsed to ~0.35). Folding everything
  into the fitted `scale*logit + bias` removes that double correction.
- **The calibration tuner searches a wider bias range.**
  `--tune-cnn-calibration` now grid-searches donor/acceptor bias over
  `{-4 ... 4}` so the optimum is not pinned to a grid edge. Because opening an
  intron forces the decoded path through both a donor and an acceptor state, the
  effective knobs are a shared `scale` and the *sum* `donor_bias + acceptor_bias`
  (an "open an intron" prior); the tuner reflects this with ties along constant
  bias sums.
- **Calibration recovers the regression but does not by itself beat the PSSM
  splice model.** After tuning, intron F1 recovers from ~0.47 to ~0.61 and exact
  accuracy from ~0.906 to ~0.941, versus the PSSM baseline of ~0.73 / ~0.956.
  The remaining gap is CNN splice discrimination, not calibration.

## Model Changes (2.2): Sparse CNN Scoring

- **CNN score generation is sparse by default.**
  The training pipeline writes donor/acceptor score TSV rows only for canonical
  splice candidates (`GT` donor positions and `AG` acceptor positions) instead
  of every base in the genome.
- **Sparse score files use the same generalized CNN checkpoint.**
  The model is still trained once across all four fission yeast species; sparse
  scoring only changes how the trained CNN is applied to each FASTA.
- **Decoding remains compatible with sparse scores.**
  `Emission_Model` already rejects noncanonical donor/acceptor motifs before
  reading CNN scores, and `Splice_CNN_Scores` leaves missing positions at the
  neutral default.
- **Dense scoring is still available for debugging.**
  Add `--dense-cnn-scores` to `train_cached_model.py` to force the old every-base
  score TSV behavior.

## Model Changes (2.1): Cached Training Pipeline

- **The full training pipeline can be refreshed with one command.**
  `src/model/training_pipeline/train_cached_model.py` trains or reloads the CNN,
  regenerates splice-score TSVs, compiles the HMM matrix trainer, and writes the
  cached HMM artifacts.
- **Trained HMM artifacts are saved locally.**
  `src/model/training_pipeline/train_hmm_matrices.cpp` writes
  `transition_matrix.json`, `emission_matrix.json`, and `metadata.json` under
  `src/model/training_pipeline/trained_models/<profile>/`.
- **CNN splice emissions are part of the cached model metadata.**
  The emission artifact stores the reusable HMM emission tables and records the
  CNN checkpoint plus train/test splice-score paths used for donor/acceptor
  emissions.
- **Generated model artifacts stay out of git.**
  The cache folders keep `.gitkeep` placeholders, while local `.pt`, score TSV,
  and trained matrix artifacts are ignored so large retrained outputs are not
  committed accidentally.
- **Training output now reports progress.**
  The CNN trainer logs FASTA/GFF loading, negative sampling, window encoding,
  epoch/batch progress, checkpoint saves, and score-writing progress so long
  runs do not appear stalled.

## Model Changes (2.0): HMM + CNN Emissions

- **Donor and acceptor emissions now come from CNN score files.**
  `DONOR_1/2/3` and `ACCEPTOR_1/2/3` still require canonical `GT`/`AG`, but
  their emission log-odds are loaded from CNN-produced per-position scores.
- **There is no PSSM fallback for splice emissions.** If CNN scores are missing,
  `Emission_Model` prints an error and throws instead of silently decoding with
  donor/acceptor PSSMs.
- **The PyTorch CNN lives under `src/model/cnn/`.**
  `splice_cnn_network.py` defines the network, while
  `train_splice_cnn_scores.py` trains or reloads the local checkpoint and writes
  score TSVs for the HMM.
- **The trained fission yeast CNN checkpoint is local.**
  The profile points to
  `src/model/cnn/trained_models/fission_yeasts_splice_cnn.pt`; this file is
  intentionally gitignored so it can be cached locally without committing model
  weights.
- **The model folder is split by responsibility.** HMM transition logic lives in
  `src/model/transition/`, HMM emissions in `src/model/emission/`, and CNN
  training/score loading in `src/model/cnn/`.
- **The active dataset is the combined fission yeast family.** The current
  profile pools `S. pombe`, `S. japonicus`, `S. octosporus`, and
  `S. cryophilus`, with 66 training chromosomes and 4 held-out evaluation
  chromosomes.

## Recent Model Updates/Features (1.3)

- **Forward-Backward posterior confidence was added.**
  `src/decoding/Forward_Backward.*` computes posterior log probabilities for
  each HMM state at each base using the same transition and emission model as
  decoding.
- **Per-base confidence is available for predicted paths.**
  `Forward_Backward::confidence(...)` returns a `vector<double>` aligned to the
  input sequence, where each value is
  `P(predicted_state_at_position_i | full sequence)`.
- **Confidence tests were added.** The unit test runner now checks empty input,
  forced-path confidence near `1.0`, and a simple ambiguous one-base posterior
  case.

## Recent Model Updates/Features (1.2)

- **Annotation filters now affect the training labels.** `GFF_Parser` groups CDS
  fragments by parent transcript, applies the active profile's CDS-length,
  intron-length, and reading-frame filters plus the current plus-strand support
  rule, and marks failed groups as `IGNORED_REGION` instead of treating them as
  ordinary intergenic sequence.
- **Ignored annotation is excluded from model fitting and holdout scoring.**
  The validation runner splits training and evaluation chromosomes into usable
  intervals around ignored regions before computing transition/emission
  probabilities or decoding held-out sequence.
- **Region codes are named constants.** The parser now exposes
  `INTERGENIC_REGION`, `CDS_REGION`, `INTRON_REGION`, and `IGNORED_REGION`,
  and parser tests include the ignored count in their region sanity check.
- **Validation output was refreshed.** The result files now use the compact
  metric tables and report usable interval-aware scores.

## Recent Model Updates/Features (1.1)

- **Splice signals use log-odds now.** Donor and acceptor windows are scored as
  `log P(window | true site) - log P(window | matched background)`. Donor
  background is non-splice `GT`; acceptor background is non-splice `AG`.
- **Canonical splice motifs are enforced.** Donors require `GT`; acceptors
  require `AG`.
- **Coding emissions are frame-specific.** `EXON_FRAME_1`, `EXON_FRAME_2`, and
  `EXON_FRAME_3` each train their own Markov5 table.
- **Stop codons are context-checked.** Only `TAA`, `TAG`, and `TGA` are legal.
- **Start codons have context scoring.** `START_CODON_1` requires full `ATG`
  context and adds a trained start-window log-odds score against non-start `ATG`
  background.
- **Viterbi applies a simple intron length constraint.** Validation learns a
  p95 intron-body cap from the training chromosomes and passes it to the
  decoder.
- **Gene starts have a small prior penalty.** In
  `validation/full_genome_validation.cpp`, `gene_start_penalty` is currently
  set to `1.0` and is subtracted on `INTERGENIC -> START_CODON_1` transitions.
- **Updated validation logging** The validation runner prints compact summary, classification, and boundary
metric tables for the held-out chromosome set.
