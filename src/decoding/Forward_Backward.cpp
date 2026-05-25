#include "Forward_Backward.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace gene_hmm {

    using namespace std;

    static Log_Prob log_sum_exp(Log_Prob a, Log_Prob b) {
        if(a == LOG_ZERO) return b;
        if(b == LOG_ZERO) return a;
        Log_Prob max_value = max(a, b);
        return max_value + log(exp(a - max_value) + exp(b - max_value));
    }

    static Log_Prob gene_entry_penalty(State previous, State current, Log_Prob penalty) {
        if(previous == State::INTERGENIC && current == State::START_CODON_1) {
            return penalty;
        }
        return 0.0;
    }

    static Forward_Backward::Probability_Matrix make_probability_matrix(size_t T) {
        Forward_Backward::Probability_Matrix matrix(T, array<Log_Prob, NUM_STATES>{});
        for(auto& row : matrix){
            for(auto& value : row){
                value = LOG_ZERO;
            }
        }
        return matrix;
    }

    static Log_Prob sequence_log_prob(
        const Forward_Backward::Probability_Matrix& forward,
        const Transition_Model::Log_Prob_Matrix& transition_log_probs)
    {
        Log_Prob total = LOG_ZERO;
        if(forward.empty()) return total;

        size_t last = forward.size() - 1;
        for(State s = State::INTERGENIC; s < State::END; s = st(idx(s)+1)){
            Log_Prob curr = forward[last][idx(s)] + transition_log_probs[idx(s)][idx(State::END)];
            total = log_sum_exp(total, curr);
        }
        return total;
    }

    Forward_Backward::Probability_Matrix Forward_Backward::posterior_log_probs(
        const vector<Nucleotide>& nucleotides,
        const Transition_Model::Log_Prob_Matrix& transition_log_probs,
        const Emission_Model& emission_model)
    {
        return posterior_log_probs(
            nucleotides,
            transition_log_probs,
            emission_model,
            0.0);
    }

    Forward_Backward::Probability_Matrix Forward_Backward::posterior_log_probs(
        const vector<Nucleotide>& nucleotides,
        const Transition_Model::Log_Prob_Matrix& transition_log_probs,
        const Emission_Model& emission_model,
        Log_Prob gene_start_penalty)
    {
        const size_t T = nucleotides.size();
        if(T == 0) return {};

        Probability_Matrix forward = make_probability_matrix(T);
        Probability_Matrix backward = make_probability_matrix(T);
        Probability_Matrix posterior = make_probability_matrix(T);

        for(State s = State::INTERGENIC; s < State::END; s = st(idx(s)+1)){
            forward[0][idx(s)] = transition_log_probs[idx(State::START)][idx(s)] +
                                 emission_model.emission_log_prob(s, 0, nucleotides);
        }

        for(size_t t = 1; t < T; t++){
            for(State s = State::INTERGENIC; s < State::END; s = st(idx(s)+1)){
                Log_Prob total = LOG_ZERO;
                for(State p = State::INTERGENIC; p < State::END; p = st(idx(p)+1)){
                    Log_Prob curr = forward[t - 1][idx(p)] +
                                    transition_log_probs[idx(p)][idx(s)] -
                                    gene_entry_penalty(p, s, gene_start_penalty);
                    total = log_sum_exp(total, curr);
                }
                forward[t][idx(s)] = total + emission_model.emission_log_prob(s, t, nucleotides);
            }
        }

        for(State s = State::INTERGENIC; s < State::END; s = st(idx(s)+1)){
            backward[T - 1][idx(s)] = transition_log_probs[idx(s)][idx(State::END)];
        }

        for(size_t t = T - 1; t > 0; t--){
            for(State p = State::INTERGENIC; p < State::END; p = st(idx(p)+1)){
                Log_Prob total = LOG_ZERO;
                for(State s = State::INTERGENIC; s < State::END; s = st(idx(s)+1)){
                    Log_Prob curr = transition_log_probs[idx(p)][idx(s)] -
                                    gene_entry_penalty(p, s, gene_start_penalty) +
                                    emission_model.emission_log_prob(s, t, nucleotides) +
                                    backward[t][idx(s)];
                    total = log_sum_exp(total, curr);
                }
                backward[t - 1][idx(p)] = total;
            }
        }

        Log_Prob log_likelihood = sequence_log_prob(forward, transition_log_probs);
        for(size_t t = 0; t < T; t++){
            for(State s = State::INTERGENIC; s < State::END; s = st(idx(s)+1)){
                posterior[t][idx(s)] = forward[t][idx(s)] + backward[t][idx(s)] - log_likelihood;
            }
        }

        return posterior;
    }

    vector<double> Forward_Backward::confidence(
        const vector<Nucleotide>& nucleotides,
        const vector<State>& predicted_states,
        const Transition_Model::Log_Prob_Matrix& transition_log_probs,
        const Emission_Model& emission_model)
    {
        return confidence(
            nucleotides,
            predicted_states,
            transition_log_probs,
            emission_model,
            0.0);
    }

    vector<double> Forward_Backward::confidence(
        const vector<Nucleotide>& nucleotides,
        const vector<State>& predicted_states,
        const Transition_Model::Log_Prob_Matrix& transition_log_probs,
        const Emission_Model& emission_model,
        Log_Prob gene_start_penalty)
    {
        if(nucleotides.size() != predicted_states.size()){
            throw runtime_error("Predicted state length does not match nucleotide length.");
        }

        Probability_Matrix posterior = posterior_log_probs(
            nucleotides,
            transition_log_probs,
            emission_model,
            gene_start_penalty);

        vector<double> confidences(nucleotides.size(), 0.0);
        for(size_t t = 0; t < nucleotides.size(); t++){
            double value = exp(posterior[t][idx(predicted_states[t])]);
            confidences[t] = min(1.0, max(0.0, value));
        }
        return confidences;
    }
}

