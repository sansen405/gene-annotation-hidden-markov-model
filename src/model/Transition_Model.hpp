#pragma once

#include "../topology/Topology.hpp"
#include "../parsers/FNA_Parser.hpp"
#include <array>
#include <cstdint>
#include <vector>

namespace gene_hmm {

    using namespace std;

    class Transition_Model {
    public:
        //Count_Matrix[i][j] = count(state i -> state j)
        //Row_Sum_Vector[i] = sum_k count(state i -> state k)
        //Log_Prob_Matrix[i][j] = log P(state j | state i)
        using Count_Matrix = array<array<uint64_t, NUM_STATES>, NUM_STATES>;
        using Row_Sum_Vector = array<uint64_t, NUM_STATES>;
        using Log_Prob_Matrix = array<array<Log_Prob, NUM_STATES>, NUM_STATES>;

        //(1) Count all i -> j transitions
        //(2) Count all outgoing transitions from each state
        //(3) Compute the log probability matrix.
        static Count_Matrix count_bigrams(const vector<State>& states, const vector<Chromosome_Range>& chromosomes);
        static Row_Sum_Vector count_outgoing(const vector<State>& states, const vector<Chromosome_Range>& chromosomes);
        static Log_Prob_Matrix compute_log_probs(const vector<State>& states, const vector<Chromosome_Range>& chromosomes,
                                               double alpha);
    };

}
