#include "Emission_Model.hpp"
#include "../../genome_profiles/Genome_Profile.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cctype>

namespace gene_hmm{

    using namespace std;

    Emission_Model::Emission_Model(){
        const auto& em = profile.emissions;
        start_window_left     = (em.count("START_CODON") && em.at("START_CODON").window_left + em.at("START_CODON").window_right > 0) ? em.at("START_CODON").window_left     : 6;
        start_window_right    = (em.count("START_CODON") && em.at("START_CODON").window_left + em.at("START_CODON").window_right > 0) ? em.at("START_CODON").window_right    : 9;
        donor_window_left     = em.count("DONOR")    ? em.at("DONOR").window_left     : 3;
        donor_window_right    = em.count("DONOR")    ? em.at("DONOR").window_right    : 6;
        acceptor_window_left  = em.count("ACCEPTOR") ? em.at("ACCEPTOR").window_left  : 15;
        acceptor_window_right = em.count("ACCEPTOR") ? em.at("ACCEPTOR").window_right : 3;

        const array<Log_Prob, NUM_NUCLEOTIDES> uniform_row = {log(0.25), log(0.25), log(0.25), log(0.25)};
        const array<Log_Prob, NUM_NUCLEOTIDES> neutral_log_odds_row = {0.0, 0.0, 0.0, 0.0};
        for(auto& row : intergenic_lp) row = uniform_row;
        for(auto& row : intron_lp) row = uniform_row;
        for(auto& row : exon_lp) row = uniform_row;
        for(auto& frame_table : exon_frame_lp)
            for(auto& row : frame_table) row = uniform_row;
        start_codon_lp.assign(start_window_left + start_window_right, neutral_log_odds_row);
        donor_lp.assign(donor_window_left + donor_window_right, neutral_log_odds_row);
        acceptor_lp.assign(acceptor_window_left + acceptor_window_right, neutral_log_odds_row);
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

    static bool has_splice_consensus(
        const vector<Nucleotide>& nucleotides,
        size_t pos,
        Splice_Signal signal)
    {
        if(signal == Splice_Signal::DONOR){
            return pos + 1 < nucleotides.size() && nucleotides[pos] == Nucleotide::G &&
                nucleotides[pos + 1] == Nucleotide::T;
        }

        if(signal == Splice_Signal::START_CODON){
            return pos + 2 < nucleotides.size() && nucleotides[pos] == Nucleotide::A &&
                nucleotides[pos + 1] == Nucleotide::T && nucleotides[pos + 2] == Nucleotide::G;
        }

        return pos > 0 &&  nucleotides[pos - 1] == Nucleotide::A && nucleotides[pos] == Nucleotide::G;
    }

    static size_t encode_5mer(const vector<Nucleotide>& nucs, size_t pos){
        size_t context = 0;
        for(size_t i = 0; i < 5; ++i){
            context = context * NUM_NUCLEOTIDES + (static_cast<size_t>(nucs[pos - 5 + i]) - 1);
        }
        return context;
    }

    static size_t exon_frame_index(State state){
        switch(state){
            case State::EXON_FRAME_1: return 0;
            case State::EXON_FRAME_2: return 1;
            case State::EXON_FRAME_3: return 2;
            default: return 0;
        }
    }

    static Log_Prob deterministic_log_prob(
        State state,
        size_t pos,
        const vector<Nucleotide>& nucleotides)
    {
        if(state == State::STOP_CODON_3){
            if(pos < 2 || nucleotides[pos - 2] != Nucleotide::T) return LOG_ZERO;
            bool is_taa = nucleotides[pos - 1] == Nucleotide::A && nucleotides[pos] == Nucleotide::A;
            bool is_tag = nucleotides[pos - 1] == Nucleotide::A && nucleotides[pos] == Nucleotide::G;
            bool is_tga = nucleotides[pos - 1] == Nucleotide::G && nucleotides[pos] == Nucleotide::A;
            return (is_taa || is_tag || is_tga) ? 0.0 : LOG_ZERO;
        }

        return Emission_Model::get_deterministic_log_prob(state, nucleotides[pos]);
    }

    static Log_Prob start_codon_signal_log_prob(
        const Emission_Model::PSSM_Log_Prob& log_probs,
        const vector<Nucleotide>& nucleotides,
        size_t pos,
        size_t window_left,
        size_t window_right)
    {
        const size_t window_size = window_left + window_right;
        if(!has_splice_consensus(nucleotides, pos, Splice_Signal::START_CODON)) return LOG_ZERO;
        if(pos < window_left) return 0.0;
        if(pos + window_right > nucleotides.size()) return 0.0;
        if(log_probs.size() < window_size) return 0.0;

        Log_Prob total = 0.0;
        for(size_t col = 0; col < window_size; ++col){
            size_t nuc = idx(nucleotides[pos - window_left + col]);
            total += log_probs[col][nuc];
        }
        return total;
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

    Emission_Model::PSSM_Count Emission_Model::count_pssm_background_emissions(
        const vector<State>& state_sequence,
        const vector<Nucleotide>& nucleotides_sequence,
        const vector<Chromosome_Range>& chromosome_range,
        const vector<State>& target_states,
        size_t window_left,
        size_t window_right,
        Splice_Signal signal)
    {
        const size_t window_size = window_left + window_right;
        vector<bool> is_target = make_target_set(target_states);
        PSSM_Count counts(window_size, array<uint64_t, NUM_NUCLEOTIDES>{});

        for(const auto& chr: chromosome_range){
            for(size_t pos = chr.start; pos < chr.end; ++pos){
                if(is_target[static_cast<size_t>(state_sequence[pos])]) continue;
                if(pos < chr.start + window_left) continue;
                if(pos + window_right > chr.end) continue;
                if(!has_splice_consensus(nucleotides_sequence, pos, signal)) continue;

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

    Emission_Model::PSSM_Log_Prob Emission_Model::compute_pssm_log_odds(
        const PSSM_Count& site_counts,
        const PSSM_Count& background_counts)
    {
        const double alpha = profile.emission_alpha;
        size_t count_size = site_counts.size();
        PSSM_Log_Prob log_odds(count_size, array<Log_Prob, NUM_NUCLEOTIDES>{});

        for(size_t col = 0; col < count_size; ++col){
            uint64_t site_sum = 0;
            uint64_t background_sum = 0;
            for(size_t n = 0; n < NUM_NUCLEOTIDES; ++n){
                site_sum += site_counts[col][n];
                background_sum += background_counts[col][n];
            }

            for(size_t n = 0; n < NUM_NUCLEOTIDES; ++n){
                double site_prob = (site_counts[col][n] + alpha) / (site_sum + alpha * NUM_NUCLEOTIDES);
                double background_prob = (background_counts[col][n] + alpha) / (background_sum + alpha * NUM_NUCLEOTIDES);
                log_odds[col][n] = log(site_prob) - log(background_prob);
            }
        }

        return log_odds;
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
                return exon_frame_lp[exon_frame_index(state)][encode_5mer(nucleotides, t)][idx(nucleotides[t])];
            case Emission_Type::PSSM_DONOR:
                if(!has_splice_consensus(nucleotides, t, Splice_Signal::DONOR)) return LOG_ZERO;
                if(!splice_cnn_scores_loaded){
                    cerr << "CNN splice scores were not loaded before donor emissions were requested.\n";
                    throw runtime_error("Missing CNN splice scores for donor emissions.");
                }
                return splice_cnn.donor_log_odds(splice_cnn_position_offset + t);
            case Emission_Type::PSSM_ACCEPTOR:
                if(!has_splice_consensus(nucleotides, t, Splice_Signal::ACCEPTOR)) return LOG_ZERO;
                if(!splice_cnn_scores_loaded){
                    cerr << "CNN splice scores were not loaded before acceptor emissions were requested.\n";
                    throw runtime_error("Missing CNN splice scores for acceptor emissions.");
                }
                return splice_cnn.acceptor_log_odds(splice_cnn_position_offset + t);
            case Emission_Type::DETERMINISTIC:
                if(state == State::START_CODON_1){
                    return start_codon_signal_log_prob(start_codon_lp, nucleotides, t, start_window_left, start_window_right);
                }
                return deterministic_log_prob(state, t, nucleotides);
        }

        return LOG_ZERO;
    }

    void Emission_Model::load_splice_cnn_scores(const string& score_path, size_t sequence_length) {
        if(score_path.empty()){
            cerr << "CNN splice score path is empty; donor/acceptor emissions cannot use CNN predictions.\n";
            throw runtime_error("Missing CNN splice score path.");
        }
        if(!splice_cnn.load_scores(score_path, sequence_length)){
            cerr << "CNN splice scores did not load from: " << score_path << "\n";
            throw runtime_error("Missing CNN splice score file.");
        }
        splice_cnn_scores_loaded = true;
    }

    void Emission_Model::set_splice_cnn_position_offset(size_t offset) {
        splice_cnn_position_offset = offset;
    }

}