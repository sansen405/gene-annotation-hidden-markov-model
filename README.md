# Gene_Analysis_HMM

A C++ Hidden Markov Model for ab-initio gene prediction in eukaryotic genomes.
A 21-state topology (intergenic, start/stop codons, frame-tracked exons,
donor/acceptor splice sites, frame-tracked introns) is trained from a genome
FASTA plus a GFF annotation. Per-genome behavior is driven by a JSON profile in
`src/genome_profiles/`. Donor/acceptor splice sites and translation starts are
scored by PyTorch CNNs; the HMM consumes their per-position log-odds.

**Current version: 4.1** — decodes both strands and merges predictions. See
[Version history](#version-history).

## Results (V4.1, fission-yeast holdout)

Four held-out chromosomes (one per species), both strands, 23,565,610 evaluated
bases. Reports live in [`validation/results/version_4/`](validation/results/version_4/).

| Metric | Value |
| --- | --- |
| Exact 21-state accuracy | 0.9701 |
| Coding F1 | 0.9740 (P 0.9687 / R 0.9793) |
| Intron F1 | 0.8169 (P 0.9267 / R 0.7304) |
| Start boundary | P 0.8428 / R 0.8409 |
| Stop boundary | P 0.9168 / R 0.9147 |
| Donor boundary | P 0.9121 / R 0.8039 |
| Acceptor boundary | P 0.9150 / R 0.8064 |
| Predicted genes | 7426 (gold 7443) |

Per-species exact 21-state accuracy: `s_pombe` 0.9755, `s_japonicus` 0.9508,
`s_octosporus` 0.9808, `s_cryophilus` 0.9765.

## Model

The HMM emits one symbol per base. Each state family draws from a different
emission model:

| State family | Emission | Detail |
| --- | --- | --- |
| `INTERGENIC` | Markov order-1 | `log P(base \| previous base)` |
| `INTRON_1/2/3` | Markov order-1 | `log P(base \| previous base)` |
| `EXON_FRAME_1/2/3` | Markov order-5, per frame | codon-periodic 1024×4 table per frame |
| `DONOR_1/2/3` | splice CNN log-odds | requires canonical `GT` |
| `ACCEPTOR_1/2/3` | splice CNN log-odds | requires canonical `AG` |
| `START_CODON_1` | start CNN log-odds | requires `ATG` |
| `START_CODON_2/3` | deterministic | second/third start bases `T`, `G` |
| `STOP_CODON_1/2/3` | deterministic | one of `TAA`, `TAG`, `TGA` |

CNN scores are **required** at decode time: donor/acceptor and start emissions
throw if their score files are not loaded (no PSSM fallback). The decoder also
applies a semi-Markov intron length model and a gene-start penalty.

## Repository layout

```text
Gene_Analysis_HMM/
├── src/
│   ├── topology/Topology.hpp            # State/Nucleotide enums, transitions
│   ├── parsers/                         # FASTA -> sequence, GFF -> per-base labels
│   ├── model/
│   │   ├── transition/                  # bigram transition matrix
│   │   ├── emission/                    # HMM emission tables + CNN bridge
│   │   ├── cnn/{splice,start}/          # PyTorch CNNs + C++ score loaders
│   │   └── training_pipeline/           # cached-model build scripts
│   ├── decoding/                        # Viterbi + Forward-Backward
│   ├── tools/                           # predict_fna, split_genome_data
│   └── genome_profiles/                 # per-profile folders (JSON + dataset scripts) + loader
├── validation/                          # holdout runner + results/
├── genome_data/                         # per-species FASTA/GFF (train/test/full)
└── frontend/                            # local React UI + Node API
```

The checked-in `fission_yeasts` profile pools four Schizosaccharomyces
assemblies (`S. pombe`, `S. japonicus`, `S. octosporus`, `S. cryophilus`):
66 training chromosomes and 4 held-out test chromosomes. FASTA/GFF files stay
per-species on disk and are combined in memory for fitting.

## Genome profile

Each profile in `src/genome_profiles/` declares input paths, held-out test
chromosomes, gene-quality filters, CNN score paths (plus `*_minus` variants for
dual-strand decode), emission choices, and smoothing. See
[`fission_yeasts.json`](src/genome_profiles/fission_yeasts/fission_yeasts.json). Set
`filters.include_minus_strand` to train and decode both strands.

## Build and run

No Makefile; build with `clang++`. `nlohmann/json.hpp` is expected at
`/opt/homebrew/include`.

Unit tests:

```sh
clang++ -std=c++17 -Isrc -I/opt/homebrew/include \
  src/main.cpp src/decoding/Viterbi.cpp src/decoding/Forward_Backward.cpp \
  src/model/transition/Transition_Model.cpp src/genome_profiles/Genome_Profile.cpp \
  src/parsers/FNA_Parser.cpp src/parsers/GFF_Parser.cpp src/model/emission/Emission_Model.cpp \
  src/model/cnn/splice/Splice_CNN_Scores.cpp src/model/cnn/start/Start_CNN_Scores.cpp \
  -o /tmp/gene_hmm_tests
/tmp/gene_hmm_tests src/genome_profiles/fission_yeasts/fission_yeasts.json
```

Holdout validation:

```sh
clang++ -std=c++17 -Isrc -I/opt/homebrew/include \
  validation/full_genome_validation.cpp src/decoding/Viterbi.cpp \
  src/model/transition/Transition_Model.cpp src/genome_profiles/Genome_Profile.cpp \
  src/parsers/FNA_Parser.cpp src/parsers/GFF_Parser.cpp src/model/emission/Emission_Model.cpp \
  src/model/cnn/splice/Splice_CNN_Scores.cpp src/model/cnn/start/Start_CNN_Scores.cpp \
  -o /tmp/full_genome_validation
/tmp/full_genome_validation --profile src/genome_profiles/fission_yeasts/fission_yeasts.json
```

The runner loads CNN scores from the profile (`splice_cnn`/`start_cnn`, with the
`*_minus` paths when `include_minus_strand` is set).

## Training pipeline

Refresh the whole cached model (CNNs + HMM matrices) in one step:

```sh
python3 src/model/training_pipeline/train_cached_model.py \
  --profile src/genome_profiles/fission_yeasts/fission_yeasts.json
```

`--cnn {splice,start}` selects which CNNs to (re)train (default both);
`--skip-cnn` skips them. Each trainer reloads an existing `.pt` checkpoint
instead of retraining, so an already-trained model is just a fast score refresh.
Scores are sparse by default (canonical `GT`/`AG`/`ATG` candidates only); add
`--dense-cnn-scores` to force every-base TSVs.

This writes `transition_matrix.json`, `emission_matrix.json`, and
`metadata.json` under `src/model/training_pipeline/trained_models/<profile>/`.
CNN checkpoints and score TSVs are local artifacts (gitignored, except committed
test scores).

## Prediction CLI

```sh
clang++ -std=c++17 -Isrc -I/opt/homebrew/include \
  src/tools/predict_fna.cpp src/decoding/Viterbi.cpp src/decoding/Forward_Backward.cpp \
  src/model/transition/Transition_Model.cpp src/genome_profiles/Genome_Profile.cpp \
  src/parsers/FNA_Parser.cpp src/parsers/GFF_Parser.cpp src/model/emission/Emission_Model.cpp \
  src/model/cnn/splice/Splice_CNN_Scores.cpp src/model/cnn/start/Start_CNN_Scores.cpp \
  -o /tmp/hmm_predict_fna

/tmp/hmm_predict_fna --profile src/genome_profiles/fission_yeasts/fission_yeasts.json \
  --fna INPUT.fna --splice-cnn-scores SPLICE.tsv --start-cnn-scores START.tsv \
  > predictions.json
```

`--splice-cnn-scores` and `--start-cnn-scores` are required (no PSSM fallback).

## Frontend

`frontend/` is a local React UI plus a Node API that accepts FASTA uploads, runs
the C++ predictor locally, and keeps runs in memory for the browser session
(no cloud storage). It supports multi-file upload, a genome track view with
per-base Forward-Backward confidence, a searchable gene table, and GFF3/CSV/BED/
FASTA export.

```sh
cd frontend
npm install
npm run dev:all
```

Vite serves the UI (usually `http://localhost:5173/`); the API listens on
`http://localhost:5174/`.

## Version history

- **4.1** — dual-strand decoding: minus-strand transcripts are scored in
  reverse-complement coordinates and merged with plus-strand genes by Viterbi
  path log-probability.
- **4.0** — translation start becomes a position-aware start CNN (replaces the
  start PSSM); a training-fit calibration bias sets the operating point.
- **3.1** — semi-Markov intron length model: charges `log P(length)` once when an
  intron closes instead of a flat geometric self-loop.
- **3.0** — calibrated splice CNN that surpasses the splice PSSM and decodes at
  `scale=1, bias=0` with no tuning step.
- **2.x** — CNN donor/acceptor emissions replace splice PSSMs; cached training
  pipeline; sparse CNN scoring; single CNN calibration layer.
- **1.x** — frame-specific exon emissions; log-odds splice/start scoring with
  canonical-motif enforcement; intron length cap and gene-start penalty;
  annotation filters; Forward-Backward posterior confidence.
