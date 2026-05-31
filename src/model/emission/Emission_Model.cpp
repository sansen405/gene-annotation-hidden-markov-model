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
        const array<Log_Prob, NUM_NUCLEOTIDES> uniform_row = {log(0.25), log(0.25), log(0.25), log(0.25)};
        for(auto& row : intergenic_lp) row = uniform_row;
        for(auto& row : intron_lp) row = uniform_row;
        for(auto& frame_table : exon_frame_lp)
            for(auto& row : frame_table) row = uniform_row;
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
        {State::DONOR_1, Emission_Type::CNN_DONOR},
        {State::DONOR_2, Emission_Type::CNN_DONOR},
        {State::DONOR_3, Emission_Type::CNN_DONOR},
        {State::ACCEPTOR_1, Emission_Type::CNN_ACCEPTOR},
        {State::ACCEPTOR_2, Emission_Type::CNN_ACCEPTOR},
        {State::ACCEPTOR_3, Emission_Type::CNN_ACCEPTOR},
        {State::START_CODON_1, Emission_Type::DETERMINISTIC},
        {State::START_CODON_2, Emission_Type::DETERMINISTIC},
        {State::START_CODON_3, Emission_Type::DETERMINISTIC},
        {State::STOP_CODON_1, Emission_Type::DETERMINISTIC},
        {State::STOP_CODON_2, Emission_Type::DETERMINISTIC},
        {State::STOP_CODON_3, Emission_Type::DETERMINISTIC},
    };

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
            case Emission_Type::CNN_DONOR:
                if(!has_splice_consensus(nucleotides, t, Splice_Signal::DONOR)) return LOG_ZERO;
                if(!splice_cnn_scores_loaded){
                    cerr << "CNN splice scores were not loaded before donor emissions were requested.\n";
                    throw runtime_error("Missing CNN splice scores for donor emissions.");
                }
                return donor_cnn_scale * splice_cnn.donor_log_odds(splice_cnn_position_offset + t) + donor_cnn_bias;
            case Emission_Type::CNN_ACCEPTOR:
                if(!has_splice_consensus(nucleotides, t, Splice_Signal::ACCEPTOR)) return LOG_ZERO;
                if(!splice_cnn_scores_loaded){
                    cerr << "CNN splice scores were not loaded before acceptor emissions were requested.\n";
                    throw runtime_error("Missing CNN splice scores for acceptor emissions.");
                }
                return acceptor_cnn_scale * splice_cnn.acceptor_log_odds(splice_cnn_position_offset + t) + acceptor_cnn_bias;
            case Emission_Type::DETERMINISTIC:
                if(state == State::START_CODON_1){
                    if(!has_splice_consensus(nucleotides, t, Splice_Signal::START_CODON)) return LOG_ZERO;
                    if(!start_cnn_scores_loaded){
                        cerr << "CNN start scores were not loaded before start emissions were requested.\n";
                        throw runtime_error("Missing CNN start scores for start emissions.");
                    }
                    return start_cnn_scale * start_cnn.start_log_odds(start_cnn_position_offset + t) + start_cnn_bias;
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

    void Emission_Model::load_splice_cnn_scores(
        const vector<string>& score_paths,
        const vector<size_t>& offsets,
        size_t sequence_length)
    {
        if(score_paths.empty()){
            cerr << "CNN splice score paths are empty; donor/acceptor emissions cannot use CNN predictions.\n";
            throw runtime_error("Missing CNN splice score paths.");
        }
        if(!splice_cnn.load_scores(score_paths, offsets, sequence_length)){
            cerr << "CNN splice scores did not load from one or more profile score paths.\n";
            throw runtime_error("Missing CNN splice score file.");
        }
        splice_cnn_scores_loaded = true;
    }

    void Emission_Model::set_splice_cnn_position_offset(size_t offset) {
        splice_cnn_position_offset = offset;
    }

    void Emission_Model::set_splice_cnn_calibration(
        Log_Prob donor_scale,
        Log_Prob donor_bias,
        Log_Prob acceptor_scale,
        Log_Prob acceptor_bias)
    {
        donor_cnn_scale = donor_scale;
        donor_cnn_bias = donor_bias;
        acceptor_cnn_scale = acceptor_scale;
        acceptor_cnn_bias = acceptor_bias;
    }

    void Emission_Model::load_start_cnn_scores(const string& score_path, size_t sequence_length) {
        if(score_path.empty()){
            cerr << "CNN start score path is empty; start emissions cannot use CNN predictions.\n";
            throw runtime_error("Missing CNN start score path.");
        }
        if(!start_cnn.load_scores(score_path, sequence_length)){
            cerr << "CNN start scores did not load from: " << score_path << "\n";
            throw runtime_error("Missing CNN start score file.");
        }
        start_cnn_scores_loaded = true;
    }

    void Emission_Model::load_start_cnn_scores(
        const vector<string>& score_paths,
        const vector<size_t>& offsets,
        size_t sequence_length)
    {
        if(score_paths.empty()){
            cerr << "CNN start score paths are empty; start emissions cannot use CNN predictions.\n";
            throw runtime_error("Missing CNN start score paths.");
        }
        if(!start_cnn.load_scores(score_paths, offsets, sequence_length)){
            cerr << "CNN start scores did not load from one or more profile score paths.\n";
            throw runtime_error("Missing CNN start score file.");
        }
        start_cnn_scores_loaded = true;
    }

    void Emission_Model::set_start_cnn_position_offset(size_t offset) {
        start_cnn_position_offset = offset;
    }

    void Emission_Model::set_start_cnn_calibration(Log_Prob start_scale, Log_Prob start_bias) {
        start_cnn_scale = start_scale;
        start_cnn_bias = start_bias;
    }

}
