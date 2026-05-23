#include "Viterbi.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace gene_hmm {

    using namespace std;

    static bool is_intron_body_state(State state) {
        return state == State::INTRON_1 ||
               state == State::INTRON_2 ||
               state == State::INTRON_3;
    }

    static bool is_acceptor_state(State state) {
        return state == State::ACCEPTOR_1 ||
               state == State::ACCEPTOR_2 ||
               state == State::ACCEPTOR_3;
    }

    static bool duration_allows_transition(
        State previous,
        State current,
        size_t previous_intron_body_length,
        size_t min_intron_body_length,
        size_t max_intron_body_length)
    {
        if(is_intron_body_state(previous) && previous == current){
            return previous_intron_body_length < max_intron_body_length;
        }
        if(is_intron_body_state(previous) && is_acceptor_state(current)){
            return previous_intron_body_length >= min_intron_body_length;
        }
        return true;
    }

    static Log_Prob gene_entry_penalty(State previous, State current, Log_Prob penalty) {
        if(previous == State::INTERGENIC && current == State::START_CODON_1) {
            return penalty;
        }
        return 0.0;
    }

    vector<State> Viterbi::decode(
        const vector<Nucleotide>& nucleotides,
        const Transition_Model::Log_Prob_Matrix& transition_log_probs,
        const Emission_Model& emission_model)
    {
        return decode(
            nucleotides,
            transition_log_probs,
            emission_model,
            0,
            numeric_limits<size_t>::max(),
            0.0);
    }

    vector<State> Viterbi::decode(
        const vector<Nucleotide>& nucleotides,
        const Transition_Model::Log_Prob_Matrix& transition_log_probs,
        const Emission_Model& emission_model,
        size_t min_intron_body_length,
        size_t max_intron_body_length)
    {
        return decode(
            nucleotides,
            transition_log_probs,
            emission_model,
            min_intron_body_length,
            max_intron_body_length,
            0.0);
    }

    vector<State> Viterbi::decode(
        const vector<Nucleotide>& nucleotides,
        const Transition_Model::Log_Prob_Matrix& transition_log_probs,
        const Emission_Model& emission_model,
        size_t min_intron_body_length,
        size_t max_intron_body_length,
        Log_Prob gene_start_penalty)
    {
        const size_t T = nucleotides.size();
        const size_t S = NUM_STATES;

        if (T == 0) return {};

        Viterbi_Matrix     V(T, array<Log_Prob, S>{});
        Backpointer_Matrix B(T, array<State, S>{});
        vector<array<size_t, S>> intron_body_length(T, array<size_t, S>{});

        for(auto& row : V){
            for(auto& value : row){
                value = LOG_ZERO;
            }
        }

        for(auto& row : B){
            for(auto& value : row){
                value = State::START;
            }
        }

        for(State s = State::INTERGENIC; s < State::END; s = st(idx(s)+1)){
            V[0][idx(s)]  = transition_log_probs[idx(State::START)][idx(s)];
            V[0][idx(s)] += emission_model.emission_log_prob(s, 0, nucleotides);
            if(is_intron_body_state(s)){
                intron_body_length[0][idx(s)] = 1;
            }
        }

        for(size_t t = 1; t < T; t++){
            for(State s = State::INTERGENIC; s < State::END; s = st(idx(s)+1)){
                Log_Prob max_prob = LOG_ZERO;
                State max_prob_prev = State::START;
                size_t max_prob_intron_body_length = 0;
                for(State p = State::INTERGENIC; p < State::END; p = st(idx(p)+1)){
                    if(!duration_allows_transition(
                        p,
                        s,
                        intron_body_length[t - 1][idx(p)],
                        min_intron_body_length,
                        max_intron_body_length)) {
                        continue;
                    }

                    Log_Prob curr_prob = V[t-1][idx(p)] +
                                         transition_log_probs[idx(p)][idx(s)] -
                                         gene_entry_penalty(p, s, gene_start_penalty) +
                                         emission_model.emission_log_prob(s, t, nucleotides);
                    if(curr_prob > max_prob){
                        max_prob = curr_prob;
                        max_prob_prev = p;
                        if(is_intron_body_state(s)){
                            max_prob_intron_body_length = (p == s) ? intron_body_length[t - 1][idx(p)] + 1 : 1;
                        } else {
                            max_prob_intron_body_length = 0;
                        }
                    }
                }
                V[t][idx(s)] = max_prob;
                B[t][idx(s)] = max_prob_prev;
                intron_body_length[t][idx(s)] = max_prob_intron_body_length;
            }
        }

        Log_Prob max_prob_final = LOG_ZERO;
        State final_state = State::START;
        for(State s = State::INTERGENIC; s < State::END; s = st(idx(s)+1)){
            Log_Prob curr_prob_final = V[T-1][idx(s)] + transition_log_probs[idx(s)][idx(State::END)];
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
