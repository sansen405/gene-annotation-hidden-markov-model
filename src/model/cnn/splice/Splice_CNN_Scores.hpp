#pragma once

#include "../../../topology/Topology.hpp"
#include <string>
#include <vector>

namespace gene_hmm {

    using namespace std;

    class Splice_CNN_Scores {
    public:
        Splice_CNN_Scores();

        bool load_scores(const string& score_path, size_t sequence_length);
        bool load_scores(const vector<string>& score_paths, const vector<size_t>& offsets, size_t sequence_length);

        Log_Prob donor_log_odds(size_t position) const;
        Log_Prob acceptor_log_odds(size_t position) const;

    private:
        vector<Log_Prob> donor_scores;
        vector<Log_Prob> acceptor_scores;
    };

}
