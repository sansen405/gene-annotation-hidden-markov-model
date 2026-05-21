#include "unit_tests/GFF_Parser_Test.hpp"

int main() {
    //Unit Test (1): GFF Parser (parse_regions + parse_states) on yeast genome
    gene_hmm::run_GFF_Parser_tests(
        "genome_data/yeast_data/GCF_000146045.2_R64_genomic.fna",
        "genome_data/yeast_data/genomic.gff"
    );

    return 0;
}
