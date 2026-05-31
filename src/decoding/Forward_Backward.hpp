#pragma once

#include "../model/transition/Transition_Model.hpp"
#include "../model/emission/Emission_Model.hpp"
#include "../topology/Topology.hpp"
#include <array>
#include <vector>

namespace gene_hmm {

    using namespace std;

    class Forward_Backward {
        public:
            using Probability_Matrix = vector<array<Log_Prob, NUM_STATES>>;

            static Probability_Matrix posterior_log_probs(
                const vector<Nucleotide>& nucleotides,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model);

            static Probability_Matrix posterior_log_probs(
                const vector<Nucleotide>& nucleotides,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model,
                Log_Prob gene_start_penalty);

            static vector<double> confidence(
                const vector<Nucleotide>& nucleotides,
                const vector<State>& predicted_states,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model);

            static vector<double> confidence(
                const vector<Nucleotide>& nucleotides,
                const vector<State>& predicted_states,
                const Transition_Model::Log_Prob_Matrix& transition_log_probs,
                const Emission_Model& emission_model,
                Log_Prob gene_start_penalty);
    };
}

