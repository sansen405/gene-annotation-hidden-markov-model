#include "unit_tests/GFF_Parser_Test.hpp"
#include "unit_tests/FNA_Parser_Test.hpp"
#include "unit_tests/Transition_Model_Test.hpp"
#include "unit_tests/Emission_Model_Test.hpp"
#include "unit_tests/Viterbi_Test.hpp"

int main() {
    gene_hmm::run_FNA_Parser_tests();

    gene_hmm::run_GFF_Parser_tests(
        "genome_data/yeast_data/GCF_000146045.2_R64_genomic.fna",
        "genome_data/yeast_data/genomic.gff"
    );

    gene_hmm::run_Transition_Model_tests();
    gene_hmm::run_Emission_Model_tests();
    gene_hmm::run_Viterbi_tests();

    return 0;
}
