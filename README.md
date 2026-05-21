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
│   ├── training/                  # (planned) trainer, gene filtering, counters
│   ├── decoding/                  # (planned) Viterbi / forward-backward
│   └── genome_profiles/
│       └── yeast.json             # S. cerevisiae S288C (RefSeq GCF_000146045.2)
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
  counts emissions for three model types and produces smoothed log-probability
  tables:
  - **Markov order-1** (`INTERGENIC`, `INTRON`) — 4×4 context→emission table.
  - **Markov order-5** (`EXON_FRAME`, per frame) — 1024×4 table (4^5 contexts).
  - **PSSM** (`DONOR`, `ACCEPTOR`) — W×4 position-specific table; donor uses a
    9-position window (3 left + 6 right), acceptor uses 18 (15 left + 3 right).
  - **Deterministic** (`START_CODON`, `STOP_CODON`) — returns 0.0 or −∞ based
    on the expected nucleotide at each codon position; no training required.

### `src/training/` and `src/decoding/`
Reserved for the training orchestrator (counting, filtering, emission fitting)
and the decoder (Viterbi / posterior decoding). These folders are currently
empty placeholders.

### `src/genome_profiles/`
JSON files that fully describe a training/decoding run for one organism:
input paths, held-out test chromosomes, gene-quality filters, per-state
emission model choices, and smoothing hyperparameters.

Example: [`src/genome_profiles/yeast.json`](src/genome_profiles/yeast.json).

### `genome_data/`
Source FASTA and GFF for *Saccharomyces cerevisiae* S288C
(RefSeq assembly `GCF_000146045.2`, R64). Referenced from `yeast.json`.

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
    "emission_alpha":   0.5
  }
}
```

## Build

A build system (CMake / Makefile) and a `main` translation unit have not been
checked in yet; build instructions will be added once that lands.
