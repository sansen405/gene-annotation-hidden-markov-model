#include "Start_CNN_Scores.hpp"
#include <fstream>
#include <stdexcept>

namespace gene_hmm {

    using namespace std;

    Start_CNN_Scores::Start_CNN_Scores() = default;

    bool Start_CNN_Scores::load_scores(const string& score_path, size_t sequence_length) {
        return load_scores(vector<string>{score_path}, vector<size_t>{0}, sequence_length);
    }

    bool Start_CNN_Scores::load_scores(
        const vector<string>& score_paths,
        const vector<size_t>& offsets,
        size_t sequence_length)
    {
        if (score_paths.size() != offsets.size()) {
            throw runtime_error("Start CNN score paths and offsets must have the same length.");
        }

        start_scores.assign(sequence_length, 0.0);

        for (size_t file_index = 0; file_index < score_paths.size(); ++file_index) {
            ifstream file(score_paths[file_index]);
            if (!file.is_open()) {
                return false;
            }

            size_t position = 0;
            Log_Prob start = 0.0;
            while (file >> position >> start) {
                size_t combined_position = offsets[file_index] + position;
                if (combined_position >= sequence_length) {
                    throw runtime_error("Start CNN score position is outside the sequence length.");
                }

                start_scores[combined_position] = start;
            }
        }

        return true;
    }

    Log_Prob Start_CNN_Scores::start_log_odds(size_t position) const {
        if (position >= start_scores.size()) {
            return 0.0;
        }

        return start_scores[position];
    }

}
