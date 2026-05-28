#include "Splice_CNN_Scores.hpp"
#include <fstream>
#include <stdexcept>

namespace gene_hmm {

    using namespace std;

    Splice_CNN_Scores::Splice_CNN_Scores() = default;

    bool Splice_CNN_Scores::load_scores(const string& score_path, size_t sequence_length) {
        donor_scores.assign(sequence_length, 0.0);
        acceptor_scores.assign(sequence_length, 0.0);

        ifstream file(score_path);
        if (!file.is_open()) {
            return false;
        }

        size_t position = 0;
        Log_Prob donor = 0.0;
        Log_Prob acceptor = 0.0;

        while (file >> position >> donor >> acceptor) {
            if (position >= sequence_length) {
                throw runtime_error("Splice CNN score position is outside the sequence length.");
            }

            donor_scores[position] = donor;
            acceptor_scores[position] = acceptor;
        }

        return true;
    }

    Log_Prob Splice_CNN_Scores::donor_log_odds(size_t position) const {
        if (position >= donor_scores.size()) {
            return 0.0;
        }

        return donor_scores[position];
    }

    Log_Prob Splice_CNN_Scores::acceptor_log_odds(size_t position) const {
        if (position >= acceptor_scores.size()) {
            return 0.0;
        }

        return acceptor_scores[position];
    }

}
