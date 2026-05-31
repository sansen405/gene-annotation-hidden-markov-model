#pragma once

#include "../../topology/Topology.hpp"
#include "../../parsers/FNA_Parser.hpp"
#include "../cnn/splice/Splice_CNN_Scores.hpp"
#include "../cnn/start/Start_CNN_Scores.hpp"
#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace gene_hmm {

    using namespace std;

    enum class Emission_Type {
        SILENT,
        MARKOV1_INTERGENIC,
        MARKOV1_INTRON,
        MARKOV5_EXON,
        CNN_DONOR,
        CNN_ACCEPTOR,
        DETERMINISTIC,
    };

    enum class Splice_Signal {
        DONOR,
        ACCEPTOR,
        START_CODON,
    };

    class Emission_Model {
    public:
        Emission_Model();

        static const size_t MARKOV5_CONTEXTS = 1024;

        using Markov1_Log_Prob = array<array<Log_Prob, NUM_NUCLEOTIDES>, NUM_NUCLEOTIDES>;
        using Markov5_Log_Prob = array<array<Log_Prob, NUM_NUCLEOTIDES>, MARKOV5_CONTEXTS>;
        using Frame_Markov5_Log_Prob = array<Markov5_Log_Prob, 3>;

        Markov1_Log_Prob intergenic_lp{};
        Markov1_Log_Prob intron_lp{};
        Frame_Markov5_Log_Prob exon_frame_lp{};

        static const unordered_map<State, Emission_Type> State_To_Emission_Type;

        using Markov1_Count = array<array<uint64_t, NUM_NUCLEOTIDES>, NUM_NUCLEOTIDES>;
        using Markov5_Count = array<array<uint64_t, NUM_NUCLEOTIDES>, MARKOV5_CONTEXTS>;

        static Markov1_Count count_markov1_emissions(
            const vector<State>& states,
            const vector<Nucleotide>& nucleotides,
            const vector<Chromosome_Range>& chromosomes,
            const vector<State>& target_states);

        static Markov5_Count count_markov5_emissions(
            const vector<State>& states,
            const vector<Nucleotide>& nucleotides,
            const vector<Chromosome_Range>& chromosomes,
            const vector<State>& target_states);

        static Markov1_Log_Prob compute_markov1_log_probs(const Markov1_Count& counts);
        static Markov5_Log_Prob compute_markov5_log_probs(const Markov5_Count& counts);

        static Log_Prob get_deterministic_log_prob(State state, Nucleotide nuc);

        Log_Prob emission_log_prob(State state, size_t t, const vector<Nucleotide>& nucleotides) const;

        void load_splice_cnn_scores(const string& score_path, size_t sequence_length);
        void load_splice_cnn_scores(
            const vector<string>& score_paths,
            const vector<size_t>& offsets,
            size_t sequence_length);
        void set_splice_cnn_position_offset(size_t offset);
        void set_splice_cnn_calibration(
            Log_Prob donor_scale,
            Log_Prob donor_bias,
            Log_Prob acceptor_scale,
            Log_Prob acceptor_bias);

        void load_start_cnn_scores(const string& score_path, size_t sequence_length);
        void load_start_cnn_scores(
            const vector<string>& score_paths,
            const vector<size_t>& offsets,
            size_t sequence_length);
        void set_start_cnn_position_offset(size_t offset);
        void set_start_cnn_calibration(Log_Prob start_scale, Log_Prob start_bias);

    private:
        Splice_CNN_Scores splice_cnn;
        bool splice_cnn_scores_loaded = false;
        size_t splice_cnn_position_offset = 0;
        Log_Prob donor_cnn_scale = 1.0;
        Log_Prob donor_cnn_bias = 0.0;
        Log_Prob acceptor_cnn_scale = 1.0;
        Log_Prob acceptor_cnn_bias = 0.0;

        Start_CNN_Scores start_cnn;
        bool start_cnn_scores_loaded = false;
        size_t start_cnn_position_offset = 0;
        Log_Prob start_cnn_scale = 1.0;
        Log_Prob start_cnn_bias = 0.0;
    };

}
