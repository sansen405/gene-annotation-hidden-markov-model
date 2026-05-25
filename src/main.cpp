#include "genome_profiles/Genome_Profile.hpp"
#include "unit_tests/GFF_Parser_Test.hpp"
#include "unit_tests/FNA_Parser_Test.hpp"
#include "unit_tests/Transition_Model_Test.hpp"
#include "unit_tests/Emission_Model_Test.hpp"
#include "unit_tests/Viterbi_Test.hpp"
#include "unit_tests/Forward_Backward_Test.hpp"

int main() {
    gene_hmm::profile = gene_hmm::Genome_Profile::load("src/genome_profiles/yeast.json");

    gene_hmm::run_FNA_Parser_tests();

    gene_hmm::run_GFF_Parser_tests(
        gene_hmm::profile.fasta_path,
        gene_hmm::profile.gff_path
    );

    gene_hmm::run_Transition_Model_tests();

    gene_hmm::run_Emission_Model_tests();
    
    gene_hmm::run_Viterbi_tests();

    gene_hmm::run_Forward_Backward_tests();

    return 0;
}
