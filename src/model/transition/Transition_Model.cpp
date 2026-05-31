#include "Transition_Model.hpp"
#include "../../genome_profiles/Genome_Profile.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cctype>
#include <cmath>

namespace gene_hmm{
    using namespace std; 

    //chromosome ranges prevent counting transitions across chromosome boundaries
    Transition_Model::Count_Matrix Transition_Model::count_bigrams(const vector<State> & state_sequence, const vector<Chromosome_Range>& chromosome_range){
        Count_Matrix bigram_count = {};
        for (const auto& chromosome: chromosome_range){
            bigram_count[idx(State::START)][idx(state_sequence[chromosome.start])]++;
            for(size_t pos = chromosome.start; pos < chromosome.end - 1; ++pos){
                bigram_count[idx(state_sequence[pos])][idx(state_sequence[pos + 1])]++;
            }

            bigram_count[idx(state_sequence[chromosome.end - 1])][idx(State::END)]++;
        }
        return bigram_count;
    }

    Transition_Model::Row_Sum_Vector Transition_Model::count_outgoing(const vector<State> & state_sequence, const vector<Chromosome_Range>& chromosome_range){
        Row_Sum_Vector outgoing_count = {};
        for (const auto& chromosome: chromosome_range){
            outgoing_count[idx(State::START)]++;
            for(size_t pos = chromosome.start; pos < chromosome.end; ++pos){
                outgoing_count[idx(state_sequence[pos])]++;
            }
        }
        return outgoing_count;
    }

    Transition_Model::Log_Prob_Matrix Transition_Model::compute_log_probs(const vector<State> & state_sequence, const vector<Chromosome_Range>& chromosome_range){
        const double alpha = profile.transition_alpha;
        Log_Prob_Matrix prob_matrix = {};
        Count_Matrix bigram_count = Transition_Model::count_bigrams(state_sequence, chromosome_range);
        for(auto& row : prob_matrix){
            for(auto& value : row){
                value = LOG_ZERO;
            }
        }

        for(size_t i = 0; i < bigram_count.size(); ++i){
            State from = st(i);
            const auto& allowed = Transitions.at(from);
            if(allowed.empty()) continue;

            uint64_t allowed_outgoing = 0;
            for(State to : allowed){
                allowed_outgoing += bigram_count[i][idx(to)];
            }

            double denominator = allowed_outgoing + alpha * allowed.size();
            for(State to : allowed){
                prob_matrix[i][idx(to)] = log((bigram_count[i][idx(to)] + alpha) / denominator);
            }
        }
        return prob_matrix;
    }
}