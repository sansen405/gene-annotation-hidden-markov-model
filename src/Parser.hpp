#pragma once

#include "Config.hpp"
#include <string>
#include <vector>

namespace config {

    using namespace std;

    class Parser{
        public:
            static vector<Nucleotide> parse_sequence(string& file_path);
    };

}