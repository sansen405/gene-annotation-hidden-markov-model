src/
├── topology/
│   └── Topology.hpp           // your current Config.hpp, renamed for clarity
├── annotation/
│   ├── Sequence_Parser.{hpp,cpp}   // (existing)
│   └── GFF_Parser.{hpp,cpp}        // (existing)
├── model/
│   ├── HMM.hpp                // bundle: topology + trans + emissions + lengths
│   ├── TransitionTable.{hpp,cpp}   // counts -> log-probs
│   ├── LengthModel.{hpp,cpp}       // duration distribution per state
│   └── emission/
│       ├── EmissionModel.hpp       // abstract interface
│       ├── MarkovEmission.{hpp,cpp}    // N-th order Markov
│       ├── PSSMEmission.{hpp,cpp}      // position-specific scoring matrix
│       ├── DeterministicEmission.{hpp,cpp}  // for start/stop codons
│       └── EmissionRegistry.{hpp,cpp}  // maps state -> model
├── training/
│   ├── Trainer.{hpp,cpp}      // orchestrator
│   ├── GeneFilter.{hpp,cpp}   // gene_is_valid + reject reasons
│   └── Counter.{hpp,cpp}      // walks state path, dispatches to models
├── decoding/
│   ├── Viterbi.{hpp,cpp}
│   └── ForwardBackward.{hpp,cpp}   // for posterior decoding later
├── eval/
│   └── Evaluator.{hpp,cpp}    // Sn/Sp/F1 at nuc/exon/gene level
├── io/
│   └── ModelSerializer.{hpp,cpp}   // save/load trained HMM (JSON or binary)
└── profiles/
    ├── yeast.json
    ├── arabidopsis.json
    └── human.json