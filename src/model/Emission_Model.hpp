#pragma once

#include "../topology/Topology.hpp"
#include "../parsers/FNA_Parser.hpp"
#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace gene_hmm {

    using namespace std;

    //Which of the 3 trained tables (plus the special cases) a given State pulls from
    //Used by Emission_Model::emission_log_prob to dispatch.
    enum class Emission_Type {
        SILENT,                // START, END  (no emission, returns 0.0)
        MARKOV1_INTERGENIC,    // INTERGENIC
        MARKOV1_INTRON,        // INTRON_*
        MARKOV5_EXON,          // EXON_FRAME_*
        PSSM_DONOR,            // DONOR_*     
        PSSM_ACCEPTOR,         // ACCEPTOR_*
        DETERMINISTIC,         // START_CODON_*, STOP_CODON_*
    };

    enum class Splice_Signal {
        DONOR,
        ACCEPTOR,
        START_CODON,
    };

    class Emission_Model {
    public:
        //nuc is short for nucleotide
        //context_nuc is short for context nucleotide
        //emitted_nuc is short for emitted nucleotide

        Emission_Model();

        static const size_t MARKOV5_CONTEXTS = 1024; // 4^5 contexts

        //Markov1_Log_Prob[context_nuc][emitted_nuc]   = log P(nuc | context_nuc)
        //Markov5_Log_Prob[5mer_index][emitted_nuc]    = log P(nuc | 5mer_context)
        //PSSM_Log_Prob[window_position][emitted_nuc]  = log-odds P_site(nuc | position) / P_background(nuc | position)
        using Markov1_Log_Prob = array<array<Log_Prob, NUM_NUCLEOTIDES>, NUM_NUCLEOTIDES>;
        using Markov5_Log_Prob = array<array<Log_Prob, NUM_NUCLEOTIDES>, MARKOV5_CONTEXTS>;
        using Frame_Markov5_Log_Prob = array<Markov5_Log_Prob, 3>;
        using PSSM_Log_Prob    = vector<array<Log_Prob, NUM_NUCLEOTIDES>>;

        //Trained probability tables (one per state family)
        Markov1_Log_Prob intergenic_lp{};
        Markov1_Log_Prob intron_lp{};
        Markov5_Log_Prob exon_lp{};
        Frame_Markov5_Log_Prob exon_frame_lp{};
        PSSM_Log_Prob    start_codon_lp;
        PSSM_Log_Prob    donor_lp;
        PSSM_Log_Prob    acceptor_lp;

        //PSSM window sizes (hardcoded for now)
        size_t start_window_left;
        size_t start_window_right;
        size_t donor_window_left;
        size_t donor_window_right;
        size_t acceptor_window_left;
        size_t acceptor_window_right;

        //State -> emission family
        static const unordered_map<State, Emission_Type> State_To_Emission_Type;

        //Markov1_Count[context_nuc][emitted_nuc]    = count(context_nuc -> emitted_nuc)
        //Markov5_Count[5mer_index][emitted_nuc]     = count(5mer_context -> emitted_nuc)
        //PSSM_Count[window_position][emitted_nuc]   = count(nuc seen at window_position)
        using Markov1_Count    = array<array<uint64_t, NUM_NUCLEOTIDES>, NUM_NUCLEOTIDES>;
        using Markov5_Count    = array<array<uint64_t, NUM_NUCLEOTIDES>, MARKOV5_CONTEXTS>;
        using PSSM_Count       = vector<array<uint64_t, NUM_NUCLEOTIDES>>;

        //(1A) Count context_nuc -> emitted_nuc pairs for target states
        static Markov1_Count count_markov1_emissions(
            const vector<State>& states,
            const vector<Nucleotide>& nucleotides,
            const vector<Chromosome_Range>& chromosomes,
            const vector<State>& target_states); //ex. {INTERGENIC} or {INTRON_1, INTRON_2, INTRON_3}

        //(1B) Count 5mer_context -> emitted_nuc pairs for target states
        static Markov5_Count count_markov5_emissions(
            const vector<State>& states,
            const vector<Nucleotide>& nucleotides,
            const vector<Chromosome_Range>& chromosomes,
            const vector<State>& target_states); //ex. {EXON_FRAME_1}

        //(1C) Count nuc at each window position around target states
        static PSSM_Count count_pssm_emissions(
            const vector<State>& states,
            const vector<Nucleotide>& nucleotides,
            const vector<Chromosome_Range>& chromosomes,
            const vector<State>& target_states, //ex. {DONOR_1, DONOR_2, DONOR_3}
            size_t window_left,
            size_t window_right);

        //(1D) Count matched non-splice motif windows for PSSM background
        static PSSM_Count count_pssm_background_emissions(
            const vector<State>& states,
            const vector<Nucleotide>& nucleotides,
            const vector<Chromosome_Range>& chromosomes,
            const vector<State>& target_states,
            size_t window_left,
            size_t window_right,
            Splice_Signal signal);

        //(2) Compute log probability matrices from counts
        static Markov1_Log_Prob compute_markov1_log_probs(const Markov1_Count& counts);
        static Markov5_Log_Prob compute_markov5_log_probs(const Markov5_Count& counts);
        static PSSM_Log_Prob    compute_pssm_log_probs(const PSSM_Count& counts);
        static PSSM_Log_Prob    compute_pssm_log_odds(const PSSM_Count& site_counts, const PSSM_Count& background_counts);

        //(3) Return log P(nuc | state) for deterministic states (START/STOP codons)
        static Log_Prob get_deterministic_log_prob(State state, Nucleotide nuc);

        //(4) Returns log P(nucleotides[t] | state, context).
        Log_Prob emission_log_prob(State state, size_t t, const vector<Nucleotide>& nucleotides) const;
    };

}
