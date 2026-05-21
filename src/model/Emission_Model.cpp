#include "../model/Emission_Model.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cctype>

namespace gene_hmm{

    using namespace std;

    static vector<bool> make_target_set(const vector<State>& target_states){
        vector<bool> is_target(NUM_STATES, false);
        for(State s: target_states)
            is_target[static_cast<size_t>(s)] = true;
        return is_target;
    }

    static size_t encode_5mer(const vector<Nucleotide>& nucs, size_t pos){
        size_t context = 0;
        // build context from positions pos-5 .. pos-1
        for(size_t i = 0; i < 5; ++i){
            context = context * NUM_NUCLEOTIDES + (static_cast<size_t>(nucs[pos - 5 + i]) - 1);
        }
        return context;
    }

    Emission_Model::Markov1_Count Emission_Model::count_markov1_emissions(const vector<State>& state_sequence, 
    const vector<Nucleotide>& nucleotides_sequence, 
    const vector<Chromosome_Range>& chromosome_range, 
    const vector<State>& target_states){
        vector<bool> is_target = make_target_set(target_states);
        Markov1_Count counts = {};
        for(const auto& chr : chromosome_range){
            for(size_t pos = chr.start + 1; pos < chr.end; ++pos){
                if(!is_target[static_cast<size_t>(state_sequence[pos])]) continue;
                size_t contex = static_cast<size_t>(nucleotides_sequence[pos - 1]) - 1;
                size_t emit = static_cast<size_t>(nucleotides_sequence[pos]) - 1;
                counts[contex][emit] ++;
            }
        }
        return counts;
    }

    Emission_Model::Markov5_Count Emission_Model::count_markov5_emissions(const vector<State>& state_sequence, 
    const vector<Nucleotide>& nucleotides_sequence, 
    const vector<Chromosome_Range>& chromosome_range, 
    const vector<State>& target_states){
        vector<bool> is_target = make_target_set(target_states);
        Markov5_Count counts = {};

        for(const auto& chr : chromosome_range){
            for(size_t pos = chr.start + 5; pos < chr.end; ++ pos){
                if(!is_target[static_cast<size_t>(state_sequence[pos])]) continue;
                size_t context = encode_5mer(nucleotides_sequence, pos);
                size_t emit = static_cast<size_t>(nucleotides_sequence[pos]) - 1;
                counts[context][emit] ++;
            }
        }
        return counts;
    }

    Emission_Model::PSSM_Count Emission_Model::count_pssm_emissions(const vector<State>& state_sequence, 
    const vector<Nucleotide>& nucleotides_sequence, 
    const vector<Chromosome_Range>& chromosome_range, 
    const vector<State>& target_states, 
    size_t window_left, 
    size_t window_right){
        const size_t window_size = window_left + window_right;
        vector<bool> is_target = make_target_set(target_states);
        PSSM_Count counts(window_size, array<uint64_t, NUM_NUCLEOTIDES>{});
        for(const auto& chr: chromosome_range){
            for(size_t pos = chr.start; pos < chr.end; ++pos){
                if (!is_target[static_cast<size_t>(state_sequence[pos])]) continue;
                if(pos < chr.start + window_left) continue;
                if(pos + window_right > chr.end) continue;
                for(size_t col = 0; col < window_size; ++col){
                    size_t nucleotide = static_cast<size_t>(nucleotides_sequence[pos - window_left + col]) - 1;
                    counts[col][nucleotide]++;
                }
            }
        }
        return counts;
    }

    Emission_Model::Markov1_Log_Prob Emission_Model::compute_markov1_log_probs(const Markov1_Count& counts, double alpha){
        Markov1_Log_Prob log_prob = {};
        for(size_t context = 0; context < NUM_NUCLEOTIDES; ++ context){
            uint64_t row_sum = 0;
            for(size_t increment = 0; increment < NUM_NUCLEOTIDES; ++increment){
                row_sum += counts[context][increment];
            }
            for(size_t increment = 0; increment < NUM_NUCLEOTIDES; ++increment){
                log_prob[context][increment] = log((counts[context][increment] + alpha)/(row_sum + alpha*NUM_NUCLEOTIDES));
            }
        }
        return log_prob;
    }

    Emission_Model::Markov5_Log_Prob Emission_Model::compute_markov5_log_probs(const Markov5_Count& counts, double alpha){
        Markov5_Log_Prob log_prob = {};
        for(size_t context = 0; context < MARKOV5_CONTEXTS; ++ context){
            uint64_t row_sum = 0;
            for(size_t increment = 0; increment < NUM_NUCLEOTIDES; ++increment){
                row_sum += counts[context][increment];
            }
            for(size_t increment = 0; increment < NUM_NUCLEOTIDES; ++increment){
                log_prob[context][increment] = log((counts[context][increment] + alpha)/(row_sum + alpha*NUM_NUCLEOTIDES));
            }
        }
        return log_prob;
    }

    Emission_Model::PSSM_Log_Prob Emission_Model::compute_pssm_log_probs(const PSSM_Count& counts, double alpha){
        size_t count_size = counts.size();
        PSSM_Log_Prob log_prob(count_size, array<Log_Prob, NUM_NUCLEOTIDES>{});
        for(size_t col = 0; col < count_size; ++col){
            uint64_t col_sum = 0;
            for(size_t n = 0; n < NUM_NUCLEOTIDES; ++n) col_sum += counts[col][n];
            for(size_t n = 0; n < NUM_NUCLEOTIDES; ++n){
                log_prob[col][n] = log((counts[col][n] + alpha)/(col_sum + alpha*NUM_NUCLEOTIDES));
            }
        }

        return log_prob;
    }

    Log_Prob Emission_Model::get_deterministic_log_prob(State state, Nucleotide nucleotide){
        switch(state ){
            case State::START_CODON_1: return (nucleotide == Nucleotide::A) ? 0.0 : LOG_ZERO;
            case State::START_CODON_2: return (nucleotide == Nucleotide::T) ? 0.0 : LOG_ZERO;
            case State::START_CODON_3: return (nucleotide == Nucleotide::G) ? 0.0 : LOG_ZERO;
            case State::STOP_CODON_1: return (nucleotide == Nucleotide::T) ? 0.0 : LOG_ZERO; 
            case State::STOP_CODON_2:
                return (nucleotide == Nucleotide::A || nucleotide == Nucleotide::G) ? log(0.5) : LOG_ZERO;
            case State::STOP_CODON_3:
                return (nucleotide == Nucleotide::A || nucleotide == Nucleotide::G) ? log(0.5) : LOG_ZERO;
            default: return 0.0;
        }
    }

}

