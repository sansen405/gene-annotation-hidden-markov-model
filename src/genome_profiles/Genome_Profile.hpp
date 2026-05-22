#pragma once

#include <string>
#include <vector>
#include <unordered_map>

namespace gene_hmm {

    using namespace std;

    enum class Emission_Family { MARKOV, PSSM, DETERMINISTIC };

    struct Emission_Config {
        Emission_Family type        = Emission_Family::DETERMINISTIC;
        size_t        order       = 1;      // MARKOV only
        bool          frame_tied  = false;  // MARKOV only
        size_t        window_left = 0;      // PSSM only
        size_t        window_right= 0;      // PSSM only
    };

    struct Genome_Profile {

        // ── identity ─────────────────────────────────────────────────────────
        string name;

        // ── training paths ───────────────────────────────────────────────────
        string fasta_path;
        string gff_path;
        vector<string> test_chromosomes;
        vector<string> exclude_chromosomes;

        // ── annotation filters ───────────────────────────────────────────────
        size_t min_first_cds_bp;
        size_t min_last_cds_bp;
        size_t min_intron_bp;
        bool   require_3n_cds;
        bool   include_minus_strand;

        // ── emission configs ─────────────────────────────────────────────────
        // Keys match JSON: "INTERGENIC", "INTRON", "EXON_FRAME",
        //                  "DONOR", "ACCEPTOR", "START_CODON", "STOP_CODON"
        unordered_map<string, Emission_Config> emissions;

        // ── smoothing ────────────────────────────────────────────────────────
        double transition_alpha;
        double emission_alpha;

        // Parse a JSON genome profile file and return a populated instance
        static Genome_Profile load(const string& json_path);
    };

    // Global profile — set once at startup:
    //   gene_hmm::profile = Genome_Profile::load("src/genome_profiles/yeast.json");
    extern Genome_Profile profile;

}

