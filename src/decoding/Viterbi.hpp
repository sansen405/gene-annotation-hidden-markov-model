#pragma once

#include "../model/transition/Transition_Model.hpp"
#include "../model/emission/Emission_Model.hpp"
#include "../topology/Topology.hpp"
#include <vector>
#include <array>
#include <limits>

namespace gene_hmm {

    using namespace std;
    class Viterbi {
        public:
            //V[t][s] = log prob of the best path ending in state s at position t
            //B[t][s] = the previous state p that achieved V[t][s]
            using Viterbi_Matrix = vector<array<Log_Prob, NUM_STATES>>;
            using Backpointer_Matrix = vector<array<State, NUM_STATES>>;


            static vector<State> decode(const vector<Nucleotide>& nucleotides,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model);

            static vector<State> decode(const vector<Nucleotide>& nucleotides,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model,
                size_t min_intron_body_length,
                size_t max_intron_body_length);

            static vector<State> decode(const vector<Nucleotide>& nucleotides,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model,
                size_t min_intron_body_length,
                size_t max_intron_body_length,
                Log_Prob gene_start_penalty);

            // Semi-Markov intron duration: when intron_length_log_prob is non-empty
            // (indexed by intron body length), the geometric self-loop cost is
            // dropped and log P(length) is charged once when the intron closes.
            static vector<State> decode(const vector<Nucleotide>& nucleotides,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model,
                size_t min_intron_body_length,
                size_t max_intron_body_length,
                Log_Prob gene_start_penalty,
                const vector<Log_Prob>& intron_length_log_prob);


    };
}

