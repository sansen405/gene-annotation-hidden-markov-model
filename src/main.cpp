#include "genome_profiles/Genome_Profile.hpp"
#include "unit_tests/GFF_Parser_Test.hpp"
#include "unit_tests/FNA_Parser_Test.hpp"
#include "unit_tests/Transition_Model_Test.hpp"
#include "unit_tests/Emission_Model_Test.hpp"
#include "unit_tests/Viterbi_Test.hpp"
#include "unit_tests/Forward_Backward_Test.hpp"
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " PROFILE_JSON\n";
        return 1;
    }
    gene_hmm::profile = gene_hmm::Genome_Profile::load(argv[1]);

    gene_hmm::run_FNA_Parser_tests();

    std::string train_fasta_path = gene_hmm::profile.train_fasta_path;
    std::string train_gff_path = gene_hmm::profile.train_gff_path;
    if (!gene_hmm::profile.species.empty()) {
        train_fasta_path = gene_hmm::profile.species.front().train_fasta_path;
        train_gff_path = gene_hmm::profile.species.front().train_gff_path;
    }

    gene_hmm::run_GFF_Parser_tests(train_fasta_path, train_gff_path);

    gene_hmm::run_Transition_Model_tests();

    gene_hmm::run_Emission_Model_tests();
    
    gene_hmm::run_Viterbi_tests();

    gene_hmm::run_Forward_Backward_tests();

    return 0;
}
