#include "Splice_CNN_Scores.hpp"
#include <fstream>
#include <stdexcept>

namespace gene_hmm {

    using namespace std;

    Splice_CNN_Scores::Splice_CNN_Scores() = default;

    bool Splice_CNN_Scores::load_scores(const string& score_path, size_t sequence_length) {
        return load_scores(vector<string>{score_path}, vector<size_t>{0}, sequence_length);
    }

    bool Splice_CNN_Scores::load_scores(
        const vector<string>& score_paths,
        const vector<size_t>& offsets,
        size_t sequence_length)
    {
        if (score_paths.size() != offsets.size()) {
            throw runtime_error("Splice CNN score paths and offsets must have the same length.");
        }

        donor_scores.assign(sequence_length, 0.0);
        acceptor_scores.assign(sequence_length, 0.0);

        for (size_t file_index = 0; file_index < score_paths.size(); ++file_index) {
            ifstream file(score_paths[file_index]);
            if (!file.is_open()) {
                return false;
            }

            size_t position = 0;
            Log_Prob donor = 0.0;
            Log_Prob acceptor = 0.0;
            while (file >> position >> donor >> acceptor) {
                size_t combined_position = offsets[file_index] + position;
                if (combined_position >= sequence_length) {
                    throw runtime_error("Splice CNN score position is outside the sequence length.");
                }

                donor_scores[combined_position] = donor;
                acceptor_scores[combined_position] = acceptor;
            }
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
