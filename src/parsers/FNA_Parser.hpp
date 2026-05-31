#pragma once

#include "../topology/Topology.hpp"
#include <string>
#include <vector>
#include <unordered_map>

namespace gene_hmm {

    using namespace std;

    struct Chromosome_Range {
        string name;
        size_t start; //inclusive
        size_t end; //exclusive
    };

    class FNA_Parser{
        public:
            static vector<Nucleotide> parse_sequence(const string& fasta_path);
            static size_t get_sequence_length(const string& fasta_path);
            static unordered_map<string, size_t> get_chromosome_offsets(const string& fasta_path);
            static vector<Chromosome_Range> get_chromosome_ranges(const string& fasta_path);

            static vector<Nucleotide> reverse_complement(const vector<Nucleotide>& nucleotides);
    };

}