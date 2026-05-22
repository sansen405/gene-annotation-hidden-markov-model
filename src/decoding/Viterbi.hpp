#pragma once

#include "../model/Transition_Model.hpp"
#include "../model/Emission_Model.hpp"
#include "../topology/Topology.hpp"
#include <vector>
#include <array>

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


    };
}

