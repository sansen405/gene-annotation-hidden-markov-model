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

    struct Genome_Profile {
        string name;

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

    //gene_hmm::profile = Genome_Profile::load("src/genome_profiles/yeast.json");
    extern Genome_Profile profile;

}

