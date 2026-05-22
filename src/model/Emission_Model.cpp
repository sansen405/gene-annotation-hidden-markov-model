#include "../model/Emission_Model.hpp"
#include "../genome_profiles/Genome_Profile.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cctype>

namespace gene_hmm{

    using namespace std;

    Emission_Model::Emission_Model(){
        const auto& em = profile.emissions;
        donor_window_left     = em.count("DONOR")    ? em.at("DONOR").window_left     : 3;
        donor_window_right    = em.count("DONOR")    ? em.at("DONOR").window_right    : 6;
        acceptor_window_left  = em.count("ACCEPTOR") ? em.at("ACCEPTOR").window_left  : 15;
        acceptor_window_right = em.count("ACCEPTOR") ? em.at("ACCEPTOR").window_right : 3;

        const array<Log_Prob, NUM_NUCLEOTIDES> uniform_row = {log(0.25), log(0.25), log(0.25), log(0.25)};
        for(auto& row : intergenic_lp) row = uniform_row;
        for(auto& row : intron_lp) row = uniform_row;
        for(auto& row : exon_lp) row = uniform_row;
        donor_lp.assign(donor_window_left + donor_window_right, uniform_row);
        acceptor_lp.assign(acceptor_window_left + acceptor_window_right, uniform_row);
    }

    const unordered_map<State, Emission_Type> Emission_Model::State_To_Emission_Type = {
        {State::START, Emission_Type::SILENT},
        {State::END, Emission_Type::SILENT},
        {State::INTERGENIC, Emission_Type::MARKOV1_INTERGENIC},
        {State::INTRON_1, Emission_Type::MARKOV1_INTRON},
        {State::INTRON_2, Emission_Type::MARKOV1_INTRON},
        {State::INTRON_3, Emission_Type::MARKOV1_INTRON},
        {State::EXON_FRAME_1, Emission_Type::MARKOV5_EXON},
        {State::EXON_FRAME_2, Emission_Type::MARKOV5_EXON},
        {State::EXON_FRAME_3, Emission_Type::MARKOV5_EXON},
        {State::DONOR_1, Emission_Type::PSSM_DONOR},
        {State::DONOR_2, Emission_Type::PSSM_DONOR},
        {State::DONOR_3, Emission_Type::PSSM_DONOR},
        {State::ACCEPTOR_1, Emission_Type::PSSM_ACCEPTOR},
        {State::ACCEPTOR_2, Emission_Type::PSSM_ACCEPTOR},
        {State::ACCEPTOR_3, Emission_Type::PSSM_ACCEPTOR},
        {State::START_CODON_1, Emission_Type::DETERMINISTIC},
        {State::START_CODON_2, Emission_Type::DETERMINISTIC},
        {State::START_CODON_3, Emission_Type::DETERMINISTIC},
        {State::STOP_CODON_1, Emission_Type::DETERMINISTIC},
        {State::STOP_CODON_2, Emission_Type::DETERMINISTIC},
        {State::STOP_CODON_3, Emission_Type::DETERMINISTIC},
    };

    //helper functions start
    static Log_Prob pssm_window_log_prob(
        const Emission_Model::PSSM_Log_Prob& log_probs,
        const vector<Nucleotide>& nucleotides,
        size_t pos,
        size_t window_left,
        size_t window_right)
    {
        const size_t window_size = window_left + window_right;
        if(pos < window_left) return LOG_ZERO;
        if(pos + window_right > nucleotides.size()) return LOG_ZERO;
        if(log_probs.size() < window_size) return LOG_ZERO;

        Log_Prob total = 0.0;
        for(size_t col = 0; col < window_size; ++col){
            size_t nuc = idx(nucleotides[pos - window_left + col]);
            total += log_probs[col][nuc];
        }
        return total;
    }

    static vector<bool> make_target_set(const vector<State>& target_states){
        vector<bool> is_target(NUM_STATES, false);
        for(State s: target_states)
            is_target[static_cast<size_t>(s)] = true;
        return is_target;
    }

    static size_t encode_5mer(const vector<Nucleotide>& nucs, size_t pos){
        size_t context = 0;
        for(size_t i = 0; i < 5; ++i){
            context = context * NUM_NUCLEOTIDES + (static_cast<size_t>(nucs[pos - 5 + i]) - 1);
        }
        return context;
    }
    //helper functions end

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

    Emission_Model::Markov1_Log_Prob Emission_Model::compute_markov1_log_probs(const Markov1_Count& counts){
        const double alpha = profile.emission_alpha;
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

    Emission_Model::Markov5_Log_Prob Emission_Model::compute_markov5_log_probs(const Markov5_Count& counts){
        const double alpha = profile.emission_alpha;
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

    Emission_Model::PSSM_Log_Prob Emission_Model::compute_pssm_log_probs(const PSSM_Count& counts){
        const double alpha = profile.emission_alpha;
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

    Log_Prob Emission_Model::emission_log_prob(State state, size_t t, const vector<Nucleotide>& nucleotides) const{
        if(t >= nucleotides.size()) return LOG_ZERO;

        auto it = State_To_Emission_Type.find(state);
        if(it == State_To_Emission_Type.end()) return LOG_ZERO;

        switch(it->second){
            case Emission_Type::SILENT:
                return 0.0;
            case Emission_Type::MARKOV1_INTERGENIC:
                if(t == 0) return log(1.0 / NUM_NUCLEOTIDES);
                return intergenic_lp[idx(nucleotides[t - 1])][idx(nucleotides[t])];
            case Emission_Type::MARKOV1_INTRON:
                if(t == 0) return log(1.0 / NUM_NUCLEOTIDES);
                return intron_lp[idx(nucleotides[t - 1])][idx(nucleotides[t])];
            case Emission_Type::MARKOV5_EXON:
                if(t < 5) return log(1.0 / NUM_NUCLEOTIDES);
                return exon_lp[encode_5mer(nucleotides, t)][idx(nucleotides[t])];
            case Emission_Type::PSSM_DONOR:
                return pssm_window_log_prob(donor_lp, nucleotides, t, donor_window_left, donor_window_right);
            case Emission_Type::PSSM_ACCEPTOR:
                return pssm_window_log_prob(acceptor_lp, nucleotides, t, acceptor_window_left, acceptor_window_right);
            case Emission_Type::DETERMINISTIC:
                return get_deterministic_log_prob(state, nucleotides[t]);
        }

        return LOG_ZERO;
    }

}