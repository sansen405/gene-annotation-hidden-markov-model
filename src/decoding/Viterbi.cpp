#include "Viterbi.hpp"
#include <cmath>
#include <algorithm>

namespace gene_hmm {

    using namespace std;

    vector<State> Viterbi::decode(
        const vector<Nucleotide>& nucleotides,
        const Transition_Model::Log_Prob_Matrix& transition_log_probs,
        const Emission_Model& emission_model)
    {
        const size_t T = nucleotides.size();
        const size_t S = NUM_STATES;

        if (T == 0) return {};

        Viterbi_Matrix     V(T, array<Log_Prob, S>{});
        Backpointer_Matrix B(T, array<State, S>{});

        for(State s = State::START; s < st(S); s = st(idx(s)+1)){
            V[0][idx(s)]  = transition_log_probs[idx(State::START)][idx(s)];
            V[0][idx(s)] += emission_model.emission_log_prob(s, 0, nucleotides);
        }

        for(size_t t = 1; t < T; t++){
            for(State s = State::START; s < st(S); s = st(idx(s)+1)){
                Log_Prob max_prob = LOG_ZERO;
                State max_prob_prev = State::START;
                for(State p = State::START; p < st(S); p = st(idx(p)+1)){
                    Log_Prob curr_prob = V[t-1][idx(p)]+ transition_log_probs[idx(p)][idx(s)]+ emission_model.emission_log_prob(s, t, nucleotides);
                    if(curr_prob > max_prob){
                        max_prob = curr_prob;
                        max_prob_prev = p;
                    }
                }
                V[t][idx(s)] = max_prob;
                B[t][idx(s)] = max_prob_prev;
            }
        }

        Log_Prob max_prob_final = LOG_ZERO;
        State final_state = State::START;
        for(State s = State::START; s < st(S); s = st(idx(s)+1)){
            Log_Prob curr_prob_final = V[T-1][idx(s)];
            if (curr_prob_final > max_prob_final){
                max_prob_final = curr_prob_final;
                final_state = s;
            }
        }

        vector<State> genome_annotation (T, State::START);
        genome_annotation[T-1] = final_state;
        
        State curr_state = final_state;
        for (size_t t = T-1; t > 0; t--){
            genome_annotation[t-1] = B[t][idx(curr_state)];
            curr_state = genome_annotation[t-1];
        }

        return genome_annotation;

    }

}
