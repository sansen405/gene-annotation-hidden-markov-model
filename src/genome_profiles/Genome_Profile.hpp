#pragma once

#include <string>
#include <vector>
#include <unordered_map>

namespace gene_hmm {

    using namespace std;

    enum class Emission_Family { MARKOV, PSSM, DETERMINISTIC };

    struct Emission_Config {
        Emission_Family type = Emission_Family::DETERMINISTIC;
        size_t order = 1; //MARKOV
        bool frame_tied = false; //MARKOV 
        size_t window_left = 0; //PSSM
        size_t window_right= 0; //PSSM
    };

    struct Species_Dataset {
        string name;
        string source_fasta_path;
        string source_gff_path;
        string train_fasta_path;
        string train_gff_path;
        string test_fasta_path;
        string test_gff_path;
        vector<string> test_chromosomes;
        vector<string> excluded_chromosomes;
    };

    struct Splice_CNN_Config {
        string model_path;
        vector<string> train_score_paths;
        vector<string> test_score_paths;
        double donor_scale = 1.0;
        double donor_bias = 0.0;
        double acceptor_scale = 1.0;
        double acceptor_bias = 0.0;
    };

    struct Genome_Profile {
        string name;

        vector<Species_Dataset> species;
        Splice_CNN_Config splice_cnn;

        string source_fasta_path;
        string source_gff_path;
        string train_fasta_path;
        string train_gff_path;
        string test_fasta_path;
        string test_gff_path;
        vector<string> test_chromosomes;
        vector<string> excluded_chromosomes;

        size_t min_first_cds_bp;
        size_t min_last_cds_bp;
        size_t min_intron_bp;
        bool   require_3n_cds;
        bool   include_minus_strand;

        unordered_map<string, Emission_Config> emissions;

        double transition_alpha;
        double emission_alpha;

        //Parse JSON genome profile file and return a populated instance
        static Genome_Profile load(const string& json_path);
    };

    extern Genome_Profile profile;

}

