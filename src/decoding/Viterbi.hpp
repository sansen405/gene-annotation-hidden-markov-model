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

            //when intron_length_log_prob is non-empty, intron duration is semi-markov: geometric self-loop cost dropped, log P(length) charged once at close
            static vector<State> decode(const vector<Nucleotide>& nucleotides,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model,
                size_t min_intron_body_length,
                size_t max_intron_body_length,
                Log_Prob gene_start_penalty,
                const vector<Log_Prob>& intron_length_log_prob);

            static Log_Prob path_log_prob(const vector<State>& states,
                const vector<Nucleotide>& nucleotides,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model,
                Log_Prob gene_start_penalty,
                size_t start,
                size_t end);


    };
}

