#include "src/genome_profiles/Genome_Profile.hpp"
#include "src/model/Emission_Model.hpp"
#include "src/model/Model_IO.hpp"
#include "src/model/Transition_Model.hpp"
#include "src/parsers/FNA_Parser.hpp"
#include "src/parsers/GFF_Parser.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

using namespace gene_hmm;
using namespace std;

namespace {

bool is_intron_body(State s) {
    return s == State::INTRON_1 || s == State::INTRON_2 || s == State::INTRON_3;
}

vector<size_t> collect_intron_body_lengths(
    const vector<State>& states,
    const vector<Chromosome_Range>& ranges)
{
    vector<size_t> lengths;
    for (const auto& range : ranges) {
        size_t i = range.start;
        while (i < range.end) {
            if (!is_intron_body(states[i])) { i++; continue; }
            size_t start = i;
            while (i < range.end && is_intron_body(states[i])) i++;
            lengths.push_back(i - start);
        }
    }
    return lengths;
}

size_t percentile_length(vector<size_t> lengths, double percentile) {
    if (lengths.empty()) return numeric_limits<size_t>::max();
    sort(lengths.begin(), lengths.end());
    size_t index = static_cast<size_t>(percentile * static_cast<double>(lengths.size() - 1));
    return lengths[index];
}

Emission_Model train_emissions(
    const vector<State>& states,
    const vector<Nucleotide>& nucleotides,
    const vector<Chromosome_Range>& train_ranges)
{
    Emission_Model model;

    model.intergenic_lp = Emission_Model::compute_markov1_log_probs(
        Emission_Model::count_markov1_emissions(
            states, nucleotides, train_ranges, {State::INTERGENIC}));

    model.intron_lp = Emission_Model::compute_markov1_log_probs(
        Emission_Model::count_markov1_emissions(
            states, nucleotides, train_ranges,
            {State::INTRON_1, State::INTRON_2, State::INTRON_3}));

    model.exon_lp = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states, nucleotides, train_ranges,
            {State::EXON_FRAME_1, State::EXON_FRAME_2, State::EXON_FRAME_3}));

    model.exon_frame_lp[0] = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states, nucleotides, train_ranges, {State::EXON_FRAME_1}));
    model.exon_frame_lp[1] = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states, nucleotides, train_ranges, {State::EXON_FRAME_2}));
    model.exon_frame_lp[2] = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states, nucleotides, train_ranges, {State::EXON_FRAME_3}));

    auto start_targets = vector<State>{State::START_CODON_1};
    model.start_codon_lp = Emission_Model::compute_pssm_log_odds(
        Emission_Model::count_pssm_emissions(
            states, nucleotides, train_ranges,
            start_targets, model.start_window_left, model.start_window_right),
        Emission_Model::count_pssm_background_emissions(
            states, nucleotides, train_ranges,
            start_targets, model.start_window_left, model.start_window_right,
            Splice_Signal::START_CODON));

    auto donor_targets = vector<State>{State::DONOR_1, State::DONOR_2, State::DONOR_3};
    model.donor_lp = Emission_Model::compute_pssm_log_odds(
        Emission_Model::count_pssm_emissions(
            states, nucleotides, train_ranges,
            donor_targets, model.donor_window_left, model.donor_window_right),
        Emission_Model::count_pssm_background_emissions(
            states, nucleotides, train_ranges,
            donor_targets, model.donor_window_left, model.donor_window_right,
            Splice_Signal::DONOR));

    auto acceptor_targets = vector<State>{State::ACCEPTOR_1, State::ACCEPTOR_2, State::ACCEPTOR_3};
    model.acceptor_lp = Emission_Model::compute_pssm_log_odds(
        Emission_Model::count_pssm_emissions(
            states, nucleotides, train_ranges,
            acceptor_targets, model.acceptor_window_left, model.acceptor_window_right),
        Emission_Model::count_pssm_background_emissions(
            states, nucleotides, train_ranges,
            acceptor_targets, model.acceptor_window_left, model.acceptor_window_right,
            Splice_Signal::ACCEPTOR));

    return model;
}

} // namespace

int main(int argc, char** argv) {
    string profile_path, output_path;
    for (int i = 1; i + 1 < argc; i += 2) {
        string flag = argv[i];
        if (flag == "--profile") profile_path = argv[i + 1];
        else if (flag == "--output") output_path = argv[i + 1];
    }

    if (profile_path.empty() || output_path.empty()) {
        cerr << "Usage: " << argv[0] << " --profile <genome_profile.json> --output <model.json>\n";
        return 1;
    }

    gene_hmm::profile = Genome_Profile::load(profile_path);
    cerr << "Training model for: " << gene_hmm::profile.name << "\n";

    vector<Nucleotide> nucleotides = FNA_Parser::parse_sequence(gene_hmm::profile.fasta_path);
    vector<Chromosome_Range> ranges = FNA_Parser::get_chromosome_ranges(gene_hmm::profile.fasta_path);

    vector<int> regions  = GFF_Parser::parse_regions(gene_hmm::profile.gff_path, gene_hmm::profile.fasta_path);
    vector<State> states = GFF_Parser::parse_states(regions);

    set<string> test_set(gene_hmm::profile.test_chromosomes.begin(), gene_hmm::profile.test_chromosomes.end());
    set<string> excl_set(gene_hmm::profile.exclude_chromosomes.begin(), gene_hmm::profile.exclude_chromosomes.end());

    vector<Chromosome_Range> train_ranges;
    for (const auto& r : ranges) {
        if (!excl_set.count(r.name) && !test_set.count(r.name))
            train_ranges.push_back(r);
    }

    if (train_ranges.empty()) {
        cerr << "No training chromosomes found.\n";
        return 1;
    }
    cerr << "Training on " << train_ranges.size() << " chromosomes.\n";

    Transition_Model::Log_Prob_Matrix trans = Transition_Model::compute_log_probs(states, train_ranges);
    Emission_Model emission = train_emissions(states, nucleotides, train_ranges);

    auto intron_lengths = collect_intron_body_lengths(states, train_ranges);
    size_t max_intron   = percentile_length(intron_lengths, 0.95);

    Trained_Model trained;
    trained.genome_name            = gene_hmm::profile.name;
    trained.transition_log_probs   = trans;
    trained.emission_model         = emission;
    trained.max_intron_body_length = max_intron;
    trained.gene_start_penalty     = 1.0;

    Model_IO::save(output_path, trained);
    cerr << "Model saved to: " << output_path << "\n";
    cerr << "Max intron body length (p95): " << max_intron << "\n";
    return 0;
}
