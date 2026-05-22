#include "genome_profiles/Genome_Profile.hpp"
#include "unit_tests/GFF_Parser_Test.hpp"
#include "unit_tests/Emission_Model_Test.hpp"

int main() {
    gene_hmm::profile = gene_hmm::Genome_Profile::load("src/genome_profiles/yeast.json");

    // Unit Test (1): GFF Parser (parse_regions + parse_states) on yeast genome
    gene_hmm::run_GFF_Parser_tests(
        gene_hmm::profile.fasta_path,
        gene_hmm::profile.gff_path
    );

    // Unit Test (2): Emission_Model (count + log_prob + deterministic)
    gene_hmm::run_Emission_Model_tests();

    return 0;
}
