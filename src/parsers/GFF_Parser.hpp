#pragma once

#include "../topology/Topology.hpp"
#include <string>
#include <vector>

namespace gene_hmm {

    using namespace std;

    class GFF_Parser{
        public:
            static constexpr int INTERGENIC_REGION = 0;
            static constexpr int CDS_REGION = 1;
            static constexpr int INTRON_REGION = 2;
            static constexpr int IGNORED_REGION = 3;

            static vector<int> parse_regions(string& gff_path, string& fna_path);
            static vector<State> parse_states(vector<int> regions);
    };

}