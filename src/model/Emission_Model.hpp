#pragma once

#include "../topology/Topology.hpp"
#include "../parsers/FNA_Parser.hpp"
#include <array>
#include <cstdint>
#include <vector>

namespace gene_hmm {

    using namespace std;

    class Emission_Model {
    public:
        //nuc is short for nucleotide
        //context_nuc is short for context nucleotide
        //emitted_nuc is short for emitted nucleotide

        static const size_t MARKOV5_CONTEXTS = 1024; // 4^5 contexts

        //Markov1_Count[context_nuc][emitted_nuc]    = count(context_nuc -> emitted_nuc)
        //Markov5_Count[5mer_index][emitted_nuc]     = count(5mer_context -> emitted_nuc)
        //PSSM_Count[window_position][emitted_nuc]   = count(nuc seen at window_position)
        using Markov1_Count    = array<array<uint64_t, NUM_NUCLEOTIDES>, NUM_NUCLEOTIDES>;
        using Markov5_Count    = array<array<uint64_t, NUM_NUCLEOTIDES>, MARKOV5_CONTEXTS>;
        using PSSM_Count       = vector<array<uint64_t, NUM_NUCLEOTIDES>>;

        //Markov1_Log_Prob[context_nuc][emitted_nuc]   = log P(nuc | context_nuc)
        //Markov5_Log_Prob[5mer_index][emitted_nuc]    = log P(nuc | 5mer_context)
        //PSSM_Log_Prob[window_position][emitted_nuc]  = log P(nuc | position)
        using Markov1_Log_Prob = array<array<Log_Prob, NUM_NUCLEOTIDES>, NUM_NUCLEOTIDES>;
        using Markov5_Log_Prob = array<array<Log_Prob, NUM_NUCLEOTIDES>, MARKOV5_CONTEXTS>;
        using PSSM_Log_Prob    = vector<array<Log_Prob, NUM_NUCLEOTIDES>>;

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

        //(2) Compute log probability matrices from counts
        static Markov1_Log_Prob compute_markov1_log_probs(const Markov1_Count& counts, double alpha);
        static Markov5_Log_Prob compute_markov5_log_probs(const Markov5_Count& counts, double alpha);
        static PSSM_Log_Prob    compute_pssm_log_probs(const PSSM_Count& counts, double alpha);

        //(3) Return log P(nuc | state) for deterministic states (START/STOP codons)
        static Log_Prob get_deterministic_log_prob(State state, Nucleotide nuc);
    };

}