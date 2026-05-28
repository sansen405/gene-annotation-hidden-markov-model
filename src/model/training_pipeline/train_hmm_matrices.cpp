#include "../../genome_profiles/Genome_Profile.hpp"
#include "../../model/emission/Emission_Model.hpp"
#include "../../model/transition/Transition_Model.hpp"
#include "../../parsers/FNA_Parser.hpp"
#include "../../parsers/GFF_Parser.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

using namespace gene_hmm;
using namespace std;
namespace fs = std::filesystem;

bool is_usable_region(int region) {
    return region != GFF_Parser::IGNORED_REGION;
}

vector<Chromosome_Range> split_usable_ranges(
    const vector<Chromosome_Range>& ranges,
    const vector<int>& regions)
{
    vector<Chromosome_Range> usable_ranges;

    for (const auto& range : ranges) {
        size_t pos = range.start;
        while (pos < range.end) {
            while (pos < range.end && !is_usable_region(regions[pos])) {
                pos++;
            }
            if (pos >= range.end) {
                break;
            }

            size_t start = pos;
            while (pos < range.end && is_usable_region(regions[pos])) {
                pos++;
            }
            usable_ranges.push_back({range.name, start, pos});
        }
    }

    return usable_ranges;
}

Emission_Model train_emissions(
    const vector<State>& states,
    const vector<Nucleotide>& nucleotides,
    const vector<Chromosome_Range>& train_ranges)
{
    Emission_Model model;

    model.intergenic_lp = Emission_Model::compute_markov1_log_probs(
        Emission_Model::count_markov1_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::INTERGENIC}));

    model.intron_lp = Emission_Model::compute_markov1_log_probs(
        Emission_Model::count_markov1_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::INTRON_1, State::INTRON_2, State::INTRON_3}));

    model.exon_lp = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::EXON_FRAME_1, State::EXON_FRAME_2, State::EXON_FRAME_3}));

    model.exon_frame_lp[0] = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::EXON_FRAME_1}));
    model.exon_frame_lp[1] = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::EXON_FRAME_2}));
    model.exon_frame_lp[2] = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::EXON_FRAME_3}));

    auto start_targets = vector<State>{State::START_CODON_1};
    model.start_codon_lp = Emission_Model::compute_pssm_log_odds(
        Emission_Model::count_pssm_emissions(
            states,
            nucleotides,
            train_ranges,
            start_targets,
            model.start_window_left,
            model.start_window_right),
        Emission_Model::count_pssm_background_emissions(
            states,
            nucleotides,
            train_ranges,
            start_targets,
            model.start_window_left,
            model.start_window_right,
            Splice_Signal::START_CODON));

    return model;
}

json serialize_log_prob(Log_Prob value) {
    if (!isfinite(value)) {
        return nullptr;
    }
    return value;
}

template <typename Matrix>
json serialize_matrix(const Matrix& matrix) {
    json rows = json::array();
    for (const auto& row : matrix) {
        json out_row = json::array();
        for (Log_Prob value : row) {
            out_row.push_back(serialize_log_prob(value));
        }
        rows.push_back(out_row);
    }
    return rows;
}

json serialize_pssm(const Emission_Model::PSSM_Log_Prob& matrix) {
    json rows = json::array();
    for (const auto& row : matrix) {
        json out_row = json::array();
        for (Log_Prob value : row) {
            out_row.push_back(serialize_log_prob(value));
        }
        rows.push_back(out_row);
    }
    return rows;
}

json read_splice_cnn_config(const string& profile_path) {
    ifstream input(profile_path);
    if (!input.is_open()) {
        throw runtime_error("Cannot open genome profile: " + profile_path);
    }

    json profile_json = json::parse(input);
    if (!profile_json.contains("splice_cnn")) {
        return json::object();
    }
    return profile_json["splice_cnn"];
}

string value_after_arg(int argc, char** argv, const string& name, const string& fallback) {
    for (int i = 1; i + 1 < argc; i++) {
        if (argv[i] == name) {
            return argv[i + 1];
        }
    }
    return fallback;
}

bool has_help_arg(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            return true;
        }
    }
    return false;
}

void print_usage(const string& program_name) {
    cerr << "Usage: " << program_name
         << " [--profile PATH] [--out-dir DIR]\n\n"
         << "Defaults:\n"
         << "  --profile src/genome_profiles/fission_yeasts.json\n"
         << "  --out-dir src/model/training_pipeline/trained_models/<profile-name>\n";
}

} // namespace

namespace gene_hmm {
    extern Genome_Profile profile;
}

int main(int argc, char** argv) {
    try {
        if (has_help_arg(argc, argv)) {
            print_usage(argv[0]);
            return 0;
        }

        string profile_path = value_after_arg(argc, argv, "--profile", "src/genome_profiles/fission_yeasts.json");
        gene_hmm::profile = Genome_Profile::load(profile_path);

        string default_out_dir = "src/model/training_pipeline/trained_models/" + gene_hmm::profile.name;
        string out_dir = value_after_arg(argc, argv, "--out-dir", default_out_dir);
        fs::create_directories(out_dir);

        vector<Nucleotide> training_nucleotides =
            FNA_Parser::parse_sequence(gene_hmm::profile.train_fasta_path);
        vector<Chromosome_Range> training_chromosomes =
            FNA_Parser::get_chromosome_ranges(gene_hmm::profile.train_fasta_path);

        string train_gff_path = gene_hmm::profile.train_gff_path;
        string train_fasta_path = gene_hmm::profile.train_fasta_path;
        vector<int> training_regions = GFF_Parser::parse_regions(train_gff_path, train_fasta_path);
        vector<State> training_states = GFF_Parser::parse_states(training_regions);
        vector<Chromosome_Range> train_ranges = split_usable_ranges(training_chromosomes, training_regions);

        if (training_states.size() != training_nucleotides.size()) {
            throw runtime_error("Training state length does not match nucleotide length.");
        }
        if (train_ranges.empty()) {
            throw runtime_error("No usable training intervals after applying annotation filters.");
        }

        auto transition_matrix = Transition_Model::compute_log_probs(training_states, train_ranges);
        Emission_Model emission_model = train_emissions(training_states, training_nucleotides, train_ranges);
        json splice_cnn = read_splice_cnn_config(profile_path);

        json transition_artifact = {
            {"profile", gene_hmm::profile.name},
            {"format", "log_transition_matrix_v1"},
            {"log_zero", nullptr},
            {"matrix", serialize_matrix(transition_matrix)}
        };

        json emission_artifact = {
            {"profile", gene_hmm::profile.name},
            {"format", "hmm_emission_tables_v1"},
            {"log_zero", nullptr},
            {"cnn_splice_emissions", {
                {"type", "per_position_scores"},
                {"config", splice_cnn}
            }},
            {"markov1", {
                {"intergenic", serialize_matrix(emission_model.intergenic_lp)},
                {"intron", serialize_matrix(emission_model.intron_lp)}
            }},
            {"markov5", {
                {"exon_combined", serialize_matrix(emission_model.exon_lp)},
                {"exon_frame_1", serialize_matrix(emission_model.exon_frame_lp[0])},
                {"exon_frame_2", serialize_matrix(emission_model.exon_frame_lp[1])},
                {"exon_frame_3", serialize_matrix(emission_model.exon_frame_lp[2])}
            }},
            {"pssm", {
                {"start_codon", serialize_pssm(emission_model.start_codon_lp)}
            }},
            {"windows", {
                {"start_left", emission_model.start_window_left},
                {"start_right", emission_model.start_window_right},
                {"donor_left", emission_model.donor_window_left},
                {"donor_right", emission_model.donor_window_right},
                {"acceptor_left", emission_model.acceptor_window_left},
                {"acceptor_right", emission_model.acceptor_window_right}
            }}
        };

        json metadata = {
            {"profile", gene_hmm::profile.name},
            {"train_fasta", gene_hmm::profile.train_fasta_path},
            {"train_gff", gene_hmm::profile.train_gff_path},
            {"training_chromosomes", training_chromosomes.size()},
            {"usable_training_intervals", train_ranges.size()},
            {"transition_matrix", "transition_matrix.json"},
            {"emission_matrix", "emission_matrix.json"},
            {"splice_cnn", splice_cnn}
        };

        ofstream transition_out(fs::path(out_dir) / "transition_matrix.json");
        transition_out << setw(2) << transition_artifact << "\n";

        ofstream emission_out(fs::path(out_dir) / "emission_matrix.json");
        emission_out << setw(2) << emission_artifact << "\n";

        ofstream metadata_out(fs::path(out_dir) / "metadata.json");
        metadata_out << setw(2) << metadata << "\n";

        cout << "Saved transition matrix: " << (fs::path(out_dir) / "transition_matrix.json") << "\n";
        cout << "Saved emission matrix:   " << (fs::path(out_dir) / "emission_matrix.json") << "\n";
        cout << "Saved metadata:          " << (fs::path(out_dir) / "metadata.json") << "\n";
        return 0;
    } catch (const exception& error) {
        cerr << error.what() << "\n";
        return 1;
    }
}
