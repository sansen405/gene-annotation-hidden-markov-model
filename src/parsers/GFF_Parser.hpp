#pragma once

#include "../topology/Topology.hpp"
#include <string>
#include <vector>

namespace gene_hmm {

    using namespace std;

    class GFF_Parser{
        public:
            static vector<int> parse_regions(string& gff_path, string& fna_path);
            static vector<State> parse_states(vector<int> regions);
    };

}