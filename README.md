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
│   │   ├── Transition_Model.hpp   # bigram counts -> log-prob transition matrix
│   │   ├── Transition_Model.cpp
│   │   ├── Emission_Model.hpp     # emission counts -> log-prob emission tables
│   │   └── Emission_Model.cpp
│   ├── decoding/
│   │   ├── Viterbi.hpp
│   │   └── Viterbi.cpp
│   └── genome_profiles/
│       └── yeast.json             # S. cerevisiae S288C (RefSeq GCF_000146045.2)
├── validation/
│   ├── full_genome_validation.cpp # train/test holdout validation runner
│   └── diagnostics/               # extra boundary and intron diagnostics
└── genome_data/
    └── yeast_data/
        ├── GCF_000146045.2_R64_genomic.fna
        └── genomic.gff
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
  (regions, then HMM states) aligned to the FASTA.

### `src/model/`
- **`Transition_Model`** — given the per-base state vector and chromosome
  ranges, counts bigram transitions, computes row sums, and produces a
  log-probability transition matrix with additive (`alpha`) smoothing.
- **`Emission_Model`** — given the per-base state and nucleotide vectors,
  counts emissions and produces smoothed log-probability tables:
  - **Markov order-1** (`INTERGENIC`, `INTRON`) — 4×4 context→emission table.
  - **Frame-specific Markov order-5** (`EXON_FRAME_1/2/3`) — one 1024×4 table
    per coding frame, so the model can use codon periodicity.
  - **Signal log-odds tables** (`START_CODON`, `DONOR`, `ACCEPTOR`) — short
    windows are scored against matched background windows instead of raw motif
    probability.
  - **Deterministic** (`START_CODON`, `STOP_CODON`) — returns 0.0 or −∞ based
    on the expected nucleotide context. Stop codons are restricted to `TAA`,
    `TAG`, and `TGA`.

Emission summary:

| State family | Emission type | Context/window | What is scored |
| --- | --- | --- | --- |
| `INTERGENIC` | Markov order-1 | Previous base | `log P(current base \| previous base)` learned from intergenic sequence. |
| `INTRON_1/2/3` | Markov order-1 | Previous base | `log P(current base \| previous base)` learned from intron sequence. |
| `EXON_FRAME_1/2/3` | Frame-specific Markov order-5 | Previous 5 bases | `log P(current base \| 5-base context)` with a separate table for each coding frame. |
| `DONOR_1/2/3` | Log-odds PSSM | `window_left=3`, `window_right=6` | Requires canonical `GT`, then scores the donor window against non-splice `GT` background. |
| `ACCEPTOR_1/2/3` | Log-odds PSSM | `window_left=15`, `window_right=3` | Requires canonical `AG`, then scores the acceptor window against non-splice `AG` background. |
| `START_CODON_1` | Deterministic motif + log-odds PSSM | Default `window_left=6`, `window_right=9` | Requires full `ATG`, then scores the surrounding start window against non-start `ATG` background. |
| `START_CODON_2/3` | Deterministic | Current base | Requires the second and third start-codon bases, `T` then `G`. |
| `STOP_CODON_1/2/3` | Deterministic motif | Three-base stop context | Requires one of `TAA`, `TAG`, or `TGA`. |

### `src/decoding/`
`Viterbi` decodes the best state path. The validation runner uses the extended
decoder overload with an intron body length cap and a gene-start penalty.

### `genome_data/`
Source FASTA and GFF for *Saccharomyces cerevisiae* S288C
(RefSeq assembly `GCF_000146045.2`, R64). Referenced from `yeast.json`.

### `src/genome_profiles/`
JSON files that fully describe a training/decoding run for one organism:
input paths, held-out test chromosomes, gene-quality filters, per-state
emission model choices, and smoothing hyperparameters.

Example: [`src/genome_profiles/yeast.json`](src/genome_profiles/yeast.json).

## Genome Profile Schema

Each profile in `src/genome_profiles/` follows this shape:

```json
{
  "name": "yeast",
  "description": "...",

  "training": {
    "fasta": "path/to/genome.fna",
    "gff":   "path/to/annotation.gff",
    "test_chromosomes":    ["NC_001148.4"],
    "exclude_chromosomes": ["NC_001224.1"]
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
  src/model/Transition_Model.cpp \
  src/genome_profiles/Genome_Profile.cpp \
  src/parsers/FNA_Parser.cpp \
  src/parsers/GFF_Parser.cpp \
  src/model/Emission_Model.cpp \
  -o /tmp/gene_hmm_tests

/tmp/gene_hmm_tests
```

Build and run the default yeast holdout validation:

```sh
clang++ -std=c++17 -Isrc -I/opt/homebrew/include \
  validation/full_genome_validation.cpp \
  src/decoding/Viterbi.cpp \
  src/model/Transition_Model.cpp \
  src/parsers/FNA_Parser.cpp \
  src/parsers/GFF_Parser.cpp \
  src/model/Emission_Model.cpp \
  -o /tmp/full_genome_validation

/tmp/full_genome_validation
```

Run validation on a custom genome:

```sh
/tmp/full_genome_validation \
  --name aspergillus \
  --fasta genome_data/aspergillus_data/GCF_000149205.2_ASM14920v2_genomic.fna \
  --gff genome_data/aspergillus_data/GCF_000149205.2_ASM14920v2_genomic.gff \
  --test-chromosomes NT_107008.1
```

## Recent Model Fixes (1.1)

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
