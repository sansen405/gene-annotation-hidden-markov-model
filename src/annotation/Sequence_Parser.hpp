#pragma once

#include "../topology/Topology.hpp"
#include <string>
#include <vector>
#include <unordered_map>

namespace gene_hmm {

    using namespace std;

    class Sequence_Parser{
        public:
            static vector<Nucleotide> parse_sequence(string& file_path);
            static size_t get_sequence_length(const string& fasta_path);
            static unordered_map<string, size_t> get_chromosome_offsets(const string& fasta_path);
    };

}