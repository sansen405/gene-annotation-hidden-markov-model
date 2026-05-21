#include "../model/Transition_Model.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cctype>
#include <cmath>

namespace gene_hmm{
    using namespace std; 

    //Need to include chromosome range because naively going from state[i] to state[i+1] will also count transitions between different chromosomes
    Transition_Model::Count_Matrix Transition_Model::count_bigrams(const vector<State> & state_sequence, const vector<Chromosome_Range>& chromosome_range){
        Count_Matrix bigram_count = {};
        for (const auto& chromosome: chromosome_range){
            bigram_count[static_cast<size_t>(State::START)][static_cast<size_t>(state_sequence[chromosome.start])]++;
            for(size_t pos = chromosome.start; pos < chromosome.end - 1; ++pos){
                size_t curr_state = static_cast<size_t>(state_sequence[pos]);
                size_t next_state = static_cast<size_t>(state_sequence[pos + 1]);
                bigram_count[curr_state][next_state]++;
            }

            bigram_count[static_cast<size_t>(state_sequence[chromosome.end -1])][static_cast<size_t>(State::END)]++;
        }
        return bigram_count;
    }

    Transition_Model::Row_Sum_Vector Transition_Model::count_outgoing(const vector<State> & state_sequence, const vector<Chromosome_Range>& chromosome_range){
        Row_Sum_Vector outgoing_count = {};
        for (const auto& chromosome: chromosome_range){
            outgoing_count[static_cast<size_t>(State::START)]++;
            for(size_t pos = chromosome.start; pos < chromosome.end; ++pos){
                size_t curr_state = static_cast<size_t>(state_sequence[pos]);
                outgoing_count[curr_state]++;
            }
        }
        return outgoing_count;
    }

    Transition_Model::Log_Prob_Matrix Transition_Model::compute_log_probs(const vector<State> & state_sequence, const vector<Chromosome_Range>& chromosome_range, double alpha){
        Log_Prob_Matrix prob_matrix = {};
        Count_Matrix bigram_count = Transition_Model::count_bigrams(state_sequence, chromosome_range);
        Row_Sum_Vector outgoing_count = Transition_Model::count_outgoing(state_sequence, chromosome_range);
        for(size_t i = 0; i < bigram_count.size(); ++i){
            for(size_t j = 0; j < bigram_count[i].size(); ++j){
                prob_matrix[i][j] = log((bigram_count[i][j] + alpha)/(outgoing_count[i] + alpha*NUM_STATES));
            }
        }
        return prob_matrix;
    }
}