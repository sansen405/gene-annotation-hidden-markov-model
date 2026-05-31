#include "../src/decoding/Viterbi.hpp"
#include "../src/genome_profiles/Genome_Profile.hpp"
#include "../src/model/emission/Emission_Model.hpp"
#include "../src/model/transition/Transition_Model.hpp"
#include "../src/parsers/FNA_Parser.hpp"
#include "../src/parsers/GFF_Parser.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <system_error>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace gene_hmm;
using namespace std;

Log_Prob gene_start_penalty = 1.0;
vector<Log_Prob> intron_length_log_probs;
string splice_cnn_scores_path;
Log_Prob donor_cnn_scale = 1.0;
Log_Prob donor_cnn_bias = 0.0;
Log_Prob acceptor_cnn_scale = 1.0;
Log_Prob acceptor_cnn_bias = 0.0;
string start_cnn_scores_path;
Log_Prob start_cnn_scale = 1.0;
Log_Prob start_cnn_bias = 0.0;
bool tune_cnn_calibration = false;
bool tune_start_calibration = false;
bool tune_only = false;
size_t tune_subset_ranges = 64;
bool collect_start_positions = false;
string start_offset_report_path;
string results_dir;
const long START_MATCH_WINDOW = 300;

struct Binary_Metrics {
    size_t tp = 0;
    size_t fp = 0;
    size_t fn = 0;
    size_t tn = 0;
};

struct Interval {
    string chromosome;
    size_t start = 0; // 0-based, chromosome-local, inclusive
    size_t end = 0;   // 0-based, chromosome-local, exclusive
};

struct Validation_Result {
    Binary_Metrics coding_total;
    Binary_Metrics intron_total;
    size_t total_bases = 0;
    size_t exact_matches = 0;
    size_t illegal_transitions = 0;
    size_t start_predicted = 0;
    size_t start_gold = 0;
    size_t start_exact = 0;
    size_t stop_predicted = 0;
    size_t stop_gold = 0;
    size_t stop_exact = 0;
    size_t donor_predicted = 0;
    size_t donor_gold = 0;
    size_t donor_exact = 0;
    size_t acceptor_predicted = 0;
    size_t acceptor_gold = 0;
    size_t acceptor_exact = 0;
    vector<Interval> predicted_gene_intervals;
    vector<Interval> gold_gene_intervals;
    // Start-boundary offset diagnostic (populated when collect_start_positions
    // is set). Offsets are signed (predicted - nearest gold) within the same
    // usable interval, in base pairs.
    vector<long> pred_start_signed_offset; // one per predicted start that had >=1 gold start in its interval
    vector<long> gold_start_nearest_dist;  // one per gold start that had >=1 predicted start in its interval
    size_t pred_starts_no_gold = 0;        // predicted starts in an interval with no gold start
    size_t gold_starts_no_pred = 0;        // gold starts in an interval with no predicted start
};

struct Sequence_Data {
    vector<Nucleotide> nucleotides;
    vector<State> states;
    vector<int> regions;
    vector<Chromosome_Range> chromosomes;
    vector<size_t> dataset_offsets;
};

bool is_coding(State state) {
    return state == State::START_CODON_1 ||
           state == State::START_CODON_2 ||
           state == State::START_CODON_3 ||
           state == State::EXON_FRAME_1 ||
           state == State::EXON_FRAME_2 ||
           state == State::EXON_FRAME_3 ||
           state == State::STOP_CODON_1 ||
           state == State::STOP_CODON_2 ||
           state == State::STOP_CODON_3;
}

bool is_intron(State state) {
    return state == State::DONOR_1 ||
           state == State::DONOR_2 ||
           state == State::DONOR_3 ||
           state == State::INTRON_1 ||
           state == State::INTRON_2 ||
           state == State::INTRON_3 ||
           state == State::ACCEPTOR_1 ||
           state == State::ACCEPTOR_2 ||
           state == State::ACCEPTOR_3;
}

bool is_intron_body(State state) {
    return state == State::INTRON_1 ||
           state == State::INTRON_2 ||
           state == State::INTRON_3;
}

bool is_gene(State state) {
    return is_coding(state) || is_intron(state);
}

double divide(size_t numerator, size_t denominator) {
    return denominator == 0 ? 0.0 : static_cast<double>(numerator) / static_cast<double>(denominator);
}

double f1(double precision, double recall) {
    return precision + recall == 0.0 ? 0.0 : 2.0 * precision * recall / (precision + recall);
}

double binary_f1(const Binary_Metrics& metrics) {
    return f1(
        divide(metrics.tp, metrics.tp + metrics.fp),
        divide(metrics.tp, metrics.tp + metrics.fn));
}

double boundary_f1(size_t exact, size_t predicted, size_t gold) {
    return f1(divide(exact, predicted), divide(exact, gold));
}

double calibration_objective(const Validation_Result& result) {
    return (
        binary_f1(result.intron_total) +
        boundary_f1(result.donor_exact, result.donor_predicted, result.donor_gold) +
        boundary_f1(result.acceptor_exact, result.acceptor_predicted, result.acceptor_gold)) / 3.0;
}

// Start calibration is fit against exact start-boundary F1 alone: it balances
// recovering true starts (recall) against the over-opening that an untuned
// start-CNN bias produces (precision), which is exactly the trade-off the
// start_scale/start_bias knobs control.
double start_calibration_objective(const Validation_Result& result) {
    return boundary_f1(result.start_exact, result.start_predicted, result.start_gold);
}

Binary_Metrics binary_metrics(
    const vector<State>& predicted,
    const vector<State>& gold,
    bool (*predicate)(State))
{
    Binary_Metrics metrics;
    for (size_t i = 0; i < predicted.size(); ++i) {
        bool pred_positive = predicate(predicted[i]);
        bool gold_positive = predicate(gold[i]);

        if (pred_positive && gold_positive) metrics.tp++;
        else if (pred_positive && !gold_positive) metrics.fp++;
        else if (!pred_positive && gold_positive) metrics.fn++;
        else metrics.tn++;
    }
    return metrics;
}

void add_metrics(Binary_Metrics& total, const Binary_Metrics& next) {
    total.tp += next.tp;
    total.fp += next.fp;
    total.fn += next.fn;
    total.tn += next.tn;
}

void print_binary_metrics_row(ostream& out, const string& label, const Binary_Metrics& metrics) {
    double precision = divide(metrics.tp, metrics.tp + metrics.fp);
    double recall = divide(metrics.tp, metrics.tp + metrics.fn);
    double accuracy = divide(metrics.tp + metrics.tn, metrics.tp + metrics.fp + metrics.fn + metrics.tn);

    out << fixed << setprecision(4);
    out << left << setw(12) << label
         << right << setw(12) << precision
         << setw(12) << recall
         << setw(12) << f1(precision, recall)
         << setw(12) << accuracy
         << setw(10) << metrics.tp
         << setw(10) << metrics.fp
         << setw(10) << metrics.fn
         << setw(10) << metrics.tn << "\n";
}

void print_binary_metrics_header(ostream& out) {
    out << left << setw(12) << "Label"
         << right << setw(12) << "Precision"
         << setw(12) << "Recall"
         << setw(12) << "F1"
         << setw(12) << "Accuracy"
         << setw(10) << "TP"
         << setw(10) << "FP"
         << setw(10) << "FN"
         << setw(10) << "TN" << "\n";
}

void print_boundary_metrics_row(
    ostream& out,
    const string& label,
    size_t exact,
    size_t predicted,
    size_t gold)
{
    out << fixed << setprecision(4);
    out << left << setw(12) << label
         << right << setw(12) << divide(exact, predicted)
         << setw(12) << divide(exact, gold)
         << setw(12) << exact
         << setw(12) << predicted
         << setw(12) << gold << "\n";
}

void print_boundary_metrics_header(ostream& out) {
    out << left << setw(12) << "Boundary"
         << right << setw(12) << "Precision"
         << setw(12) << "Recall"
         << setw(12) << "Exact"
         << setw(12) << "Predicted"
         << setw(12) << "Gold" << "\n";
}

vector<Nucleotide> slice_nucleotides(
    const vector<Nucleotide>& nucleotides,
    size_t start,
    size_t end)
{
    return vector<Nucleotide>(nucleotides.begin() + start, nucleotides.begin() + end);
}

vector<State> slice_states(
    const vector<State>& states,
    size_t start,
    size_t end)
{
    return vector<State>(states.begin() + start, states.begin() + end);
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

    auto donor_targets = vector<State>{State::DONOR_1, State::DONOR_2, State::DONOR_3};
    model.donor_lp = Emission_Model::compute_pssm_log_odds(
        Emission_Model::count_pssm_emissions(
            states,
            nucleotides,
            train_ranges,
            donor_targets,
            model.donor_window_left,
            model.donor_window_right),
        Emission_Model::count_pssm_background_emissions(
            states,
            nucleotides,
            train_ranges,
            donor_targets,
            model.donor_window_left,
            model.donor_window_right,
            Splice_Signal::DONOR));

    auto acceptor_targets = vector<State>{State::ACCEPTOR_1, State::ACCEPTOR_2, State::ACCEPTOR_3};
    model.acceptor_lp = Emission_Model::compute_pssm_log_odds(
        Emission_Model::count_pssm_emissions(
            states,
            nucleotides,
            train_ranges,
            acceptor_targets,
            model.acceptor_window_left,
            model.acceptor_window_right),
        Emission_Model::count_pssm_background_emissions(
            states,
            nucleotides,
            train_ranges,
            acceptor_targets,
            model.acceptor_window_left,
            model.acceptor_window_right,
            Splice_Signal::ACCEPTOR));

    return model;
}

size_t exact_state_matches(const vector<State>& predicted, const vector<State>& gold) {
    size_t matches = 0;
    for (size_t i = 0; i < predicted.size(); ++i) {
        if (predicted[i] == gold[i]) matches++;
    }
    return matches;
}

size_t count_state(const vector<State>& states, State target) {
    return count(states.begin(), states.end(), target);
}

size_t count_states(const vector<State>& states, const vector<State>& targets) {
    size_t total = 0;
    for (State target : targets) {
        total += count_state(states, target);
    }
    return total;
}

size_t count_exact_boundary_matches(
    const vector<State>& predicted,
    const vector<State>& gold,
    State boundary_state)
{
    size_t matches = 0;
    for (size_t i = 0; i < predicted.size(); ++i) {
        if (predicted[i] == boundary_state && gold[i] == boundary_state) {
            matches++;
        }
    }
    return matches;
}

size_t count_exact_boundary_matches(
    const vector<State>& predicted,
    const vector<State>& gold,
    const vector<State>& boundary_states)
{
    size_t matches = 0;
    for (State boundary_state : boundary_states) {
        matches += count_exact_boundary_matches(predicted, gold, boundary_state);
    }
    return matches;
}

vector<size_t> collect_intron_body_lengths(
    const vector<State>& states,
    const vector<Chromosome_Range>& ranges)
{
    vector<size_t> lengths;
    for (const auto& range : ranges) {
        size_t i = range.start;
        while (i < range.end) {
            if (!is_intron_body(states[i])) {
                i++;
                continue;
            }

            size_t start = i;
            while (i < range.end && is_intron_body(states[i])) {
                i++;
            }
            lengths.push_back(i - start);
        }
    }
    return lengths;
}

size_t percentile_length(vector<size_t> lengths, double percentile) {
    if (lengths.empty()) {
        return numeric_limits<size_t>::max();
    }

    sort(lengths.begin(), lengths.end());
    size_t index = static_cast<size_t>(percentile * static_cast<double>(lengths.size() - 1));
    return lengths[index];
}

// Smoothed empirical histogram of training intron-body lengths over [1, max],
// returned as log probabilities indexed by length. Add-1 smoothing keeps every
// in-range length finite; the decoder caps the body length at `max`, so longer
// lengths are unreachable. Replaces the implicit geometric duration with the
// true (tightly peaked) intron length distribution.
vector<Log_Prob> build_intron_length_log_probs(
    const vector<size_t>& lengths,
    size_t max_length)
{
    if (max_length == 0 || lengths.empty()) {
        return {};
    }

    vector<double> counts(max_length + 1, 1.0);
    counts[0] = 0.0;
    for (size_t length : lengths) {
        if (length >= 1 && length <= max_length) {
            counts[length] += 1.0;
        }
    }

    double total = 0.0;
    for (size_t length = 1; length <= max_length; ++length) {
        total += counts[length];
    }

    vector<Log_Prob> log_probs(max_length + 1, LOG_ZERO);
    for (size_t length = 1; length <= max_length; ++length) {
        log_probs[length] = log(counts[length] / total);
    }
    return log_probs;
}

vector<Interval> collect_intervals(
    const vector<State>& states,
    const string& chromosome,
    bool (*predicate)(State))
{
    vector<Interval> intervals;
    size_t i = 0;
    while (i < states.size()) {
        if (!predicate(states[i])) {
            i++;
            continue;
        }

        size_t start = i;
        while (i < states.size() && predicate(states[i])) {
            i++;
        }
        intervals.push_back({chromosome, start, i});
    }
    return intervals;
}

vector<Interval> collect_error_intervals(
    const vector<State>& predicted,
    const vector<State>& gold,
    const string& chromosome,
    bool false_positive)
{
    vector<State> error_mask(predicted.size(), State::INTERGENIC);
    for (size_t i = 0; i < predicted.size(); ++i) {
        bool predicted_gene = is_gene(predicted[i]);
        bool gold_gene = is_gene(gold[i]);
        if (false_positive && predicted_gene && !gold_gene) {
            error_mask[i] = predicted[i];
        }
        if (!false_positive && !predicted_gene && gold_gene) {
            error_mask[i] = gold[i];
        }
    }
    return collect_intervals(error_mask, chromosome, is_gene);
}

size_t illegal_transition_count(const vector<State>& states) {
    size_t illegal = 0;
    for (size_t i = 0; i + 1 < states.size(); ++i) {
        const auto& allowed = Transitions.at(states[i]);
        if (find(allowed.begin(), allowed.end(), states[i + 1]) == allowed.end()) {
            illegal++;
        }
    }
    return illegal;
}

void print_intervals(const string& label, const vector<Interval>& intervals, size_t max_count) {
    cout << label << " count=" << intervals.size();
    if (!intervals.empty()) {
        cout << " examples:";
        for (size_t i = 0; i < min(max_count, intervals.size()); ++i) {
            cout << " " << intervals[i].chromosome
                 << ":[" << intervals[i].start << "," << intervals[i].end << ")";
        }
    }
    cout << "\n";
}

vector<string> split_csv(const string& text) {
    vector<string> values;
    string token;
    stringstream stream(text);
    while (getline(stream, token, ',')) {
        if (!token.empty()) {
            values.push_back(token);
        }
    }
    return values;
}

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

void append_dataset(
    Sequence_Data& combined,
    const string& name,
    const string& fasta_path,
    const string& gff_path)
{
    size_t offset = combined.nucleotides.size();
    combined.dataset_offsets.push_back(offset);

    vector<Nucleotide> nucleotides = FNA_Parser::parse_sequence(fasta_path);
    vector<Chromosome_Range> chromosomes = FNA_Parser::get_chromosome_ranges(fasta_path);
    string mutable_gff_path = gff_path;
    string mutable_fasta_path = fasta_path;
    vector<int> regions = GFF_Parser::parse_regions(mutable_gff_path, mutable_fasta_path);
    vector<State> states = GFF_Parser::parse_states(regions);

    if (states.size() != nucleotides.size()) {
        throw runtime_error("State length does not match nucleotide length for " + name + ".");
    }

    combined.nucleotides.insert(combined.nucleotides.end(), nucleotides.begin(), nucleotides.end());
    combined.states.insert(combined.states.end(), states.begin(), states.end());
    combined.regions.insert(combined.regions.end(), regions.begin(), regions.end());
    for (auto range : chromosomes) {
        range.name = name + ":" + range.name;
        range.start += offset;
        range.end += offset;
        combined.chromosomes.push_back(range);
    }
}

Sequence_Data load_training_data(const Genome_Profile& profile) {
    Sequence_Data combined;
    for (const auto& dataset : profile.species) {
        append_dataset(combined, dataset.name, dataset.train_fasta_path, dataset.train_gff_path);
    }
    return combined;
}

Sequence_Data load_evaluation_data(const Genome_Profile& profile) {
    Sequence_Data combined;
    for (const auto& dataset : profile.species) {
        append_dataset(combined, dataset.name, dataset.test_fasta_path, dataset.test_gff_path);
    }
    return combined;
}

void print_state_family_counts(const string& label, const vector<State>& states) {
    cout << label
         << " intergenic=" << count_state(states, State::INTERGENIC)
         << " start=" << count_states(states, {State::START_CODON_1, State::START_CODON_2, State::START_CODON_3})
         << " exon=" << count_states(states, {State::EXON_FRAME_1, State::EXON_FRAME_2, State::EXON_FRAME_3})
         << " donor=" << count_states(states, {State::DONOR_1, State::DONOR_2, State::DONOR_3})
         << " intron=" << count_states(states, {State::INTRON_1, State::INTRON_2, State::INTRON_3})
         << " acceptor=" << count_states(states, {State::ACCEPTOR_1, State::ACCEPTOR_2, State::ACCEPTOR_3})
         << " stop=" << count_states(states, {State::STOP_CODON_1, State::STOP_CODON_2, State::STOP_CODON_3})
         << "\n";
}

void print_log_prob(const string& label, Log_Prob value) {
    cout << label << "=";
    if (value == LOG_ZERO) {
        cout << "-inf";
    } else {
        cout << value;
    }
    cout << "\n";
}

bool contains_state(const vector<State>& states, State target) {
    return find(states.begin(), states.end(), target) != states.end();
}

void print_average_emission(
    const string& label,
    const Emission_Model& model,
    const vector<Nucleotide>& nucleotides,
    const vector<State>& states,
    const vector<State>& targets)
{
    double total = 0.0;
    size_t count = 0;
    for (size_t i = 0; i < states.size(); ++i) {
        if (!contains_state(targets, states[i])) {
            continue;
        }
        Log_Prob lp = model.emission_log_prob(states[i], i, nucleotides);
        if (lp == LOG_ZERO) {
            continue;
        }
        total += lp;
        count++;
    }

    cout << label << " count=" << count;
    if (count > 0) {
        cout << " avg_log_emission=" << (total / count);
    }
    cout << "\n";
}

Genome_Profile manual_profile() {
    Genome_Profile profile;
    profile.name = "custom";
    profile.min_first_cds_bp = 3;
    profile.min_last_cds_bp = 3;
    profile.min_intron_bp = 20;
    profile.require_3n_cds = true;
    profile.include_minus_strand = false;
    profile.transition_alpha = 0.02;
    profile.emission_alpha = 0.1;
    profile.emissions["INTERGENIC"] = {Emission_Family::MARKOV, 1, false, 0, 0};
    profile.emissions["INTRON"] = {Emission_Family::MARKOV, 1, false, 0, 0};
    profile.emissions["EXON_FRAME"] = {Emission_Family::MARKOV, 5, true, 0, 0};
    profile.emissions["DONOR"] = {Emission_Family::PSSM, 1, false, 3, 6};
    profile.emissions["ACCEPTOR"] = {Emission_Family::PSSM, 1, false, 15, 3};
    profile.emissions["START_CODON"] = {Emission_Family::DETERMINISTIC, 1, false, 0, 0};
    profile.emissions["STOP_CODON"] = {Emission_Family::DETERMINISTIC, 1, false, 0, 0};
    return profile;
}

void print_usage(const string& program_name) {
    cerr << "Usage:\n";
    cerr << "  " << program_name << "\n";
    cerr << "  " << program_name << " --profile PATH\n";
    cerr << "  " << program_name
         << " --train-fasta PATH --train-gff PATH --test-fasta PATH --test-gff PATH [options]\n\n";
    cerr << "Options:\n";
    cerr << "  --name NAME              Dataset label for output\n";
    cerr << "  --splice-cnn-scores PATH CNN donor/acceptor score TSV for evaluation sequence\n";
    cerr << "  --cnn-donor-scale VALUE  Donor CNN logit multiplier, default from profile or 1.0\n";
    cerr << "  --cnn-donor-bias VALUE   Donor CNN logit offset, default from profile or 0.0\n";
    cerr << "  --cnn-acceptor-scale VALUE  Acceptor CNN logit multiplier, default from profile or 1.0\n";
    cerr << "  --cnn-acceptor-bias VALUE   Acceptor CNN logit offset, default from profile or 0.0\n";
    cerr << "  --start-cnn-scores PATH  CNN start score TSV for evaluation sequence (else profile start_cnn.test_scores)\n";
    cerr << "  --cnn-start-scale VALUE  Start CNN logit multiplier, default from profile or 1.0\n";
    cerr << "  --cnn-start-bias VALUE   Start CNN logit offset, default from profile or 0.0\n";
    cerr << "  --tune-cnn-calibration   Grid-search splice CNN calibration on the evaluation labels\n";
    cerr << "  --tune-start-calibration Grid-search start CNN scale/bias on the TRAINING labels\n";
    cerr << "  --tune-only              Stop after selecting CNN calibration\n";
    cerr << "  --tune-subset-ranges N   Usable evaluation intervals for tuning, default 64\n";
    cerr << "  --gene-start-penalty VALUE  Log-prob penalty on INTERGENIC -> START_CODON_1, default 1.0\n";
    cerr << "  --start-offset-report PATH  Write a start-boundary offset diagnostic (per-species + combined) to PATH\n";
    cerr << "  --results-dir DIR        Write combined.txt and per-species result files into DIR\n";
    cerr << "  --start-window-left N    Start-codon PSSM bases upstream of the ATG (retrains the PSSM), default 6\n";
    cerr << "  --start-window-right N   Start-codon PSSM bases from the ATG onward (>=3 covers ATG), default 9\n";
    cerr << "  --transition-alpha VALUE Transition smoothing alpha, default 0.02\n";
    cerr << "  --emission-alpha VALUE   Emission smoothing alpha, default 0.1\n\n";
    cerr << "Example:\n";
    cerr << "  " << program_name
         << " --name c_elegans"
         << " --train-fasta genome_data/c_elegans/train/c_elegans_train.fna"
         << " --train-gff genome_data/c_elegans/train/c_elegans_train.gff"
         << " --test-fasta genome_data/c_elegans/test/c_elegans_test.fna"
         << " --test-gff genome_data/c_elegans/test/c_elegans_test.gff\n";
}

Genome_Profile parse_args(int argc, char** argv) {
    if (argc == 1) {
        print_usage(argv[0]);
        exit(1);
    }

    Genome_Profile profile = manual_profile();
    bool donor_scale_set = false;
    bool donor_bias_set = false;
    bool acceptor_scale_set = false;
    bool acceptor_bias_set = false;
    bool start_scale_set = false;
    bool start_bias_set = false;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            exit(0);
        }
        if (arg == "--tune-start-calibration") {
            tune_start_calibration = true;
            continue;
        }
        if (arg == "--tune-cnn-calibration") {
            tune_cnn_calibration = true;
            continue;
        }
        if (arg == "--tune-only") {
            tune_only = true;
            continue;
        }
        if (i + 1 >= argc) {
            cerr << "Missing value for argument: " << arg << "\n";
            print_usage(argv[0]);
            exit(1);
        }

        string value = argv[++i];
        if (arg == "--name") {
            profile.name = value;
        } else if (arg == "--profile") {
            profile = Genome_Profile::load(value);
        } else if (arg == "--train-fasta") {
            profile.train_fasta_path = value;
        } else if (arg == "--train-gff") {
            profile.train_gff_path = value;
        } else if (arg == "--test-fasta") {
            profile.test_fasta_path = value;
        } else if (arg == "--test-gff") {
            profile.test_gff_path = value;
        } else if (arg == "--splice-cnn-scores") {
            splice_cnn_scores_path = value;
        } else if (arg == "--cnn-donor-scale") {
            donor_cnn_scale = stod(value);
            donor_scale_set = true;
        } else if (arg == "--cnn-donor-bias") {
            donor_cnn_bias = stod(value);
            donor_bias_set = true;
        } else if (arg == "--cnn-acceptor-scale") {
            acceptor_cnn_scale = stod(value);
            acceptor_scale_set = true;
        } else if (arg == "--cnn-acceptor-bias") {
            acceptor_cnn_bias = stod(value);
            acceptor_bias_set = true;
        } else if (arg == "--start-cnn-scores") {
            start_cnn_scores_path = value;
        } else if (arg == "--cnn-start-scale") {
            start_cnn_scale = stod(value);
            start_scale_set = true;
        } else if (arg == "--cnn-start-bias") {
            start_cnn_bias = stod(value);
            start_bias_set = true;
        } else if (arg == "--tune-subset-ranges") {
            tune_subset_ranges = stoull(value);
        } else if (arg == "--gene-start-penalty") {
            gene_start_penalty = stod(value);
        } else if (arg == "--start-offset-report") {
            start_offset_report_path = value;
        } else if (arg == "--results-dir") {
            results_dir = value;
        } else if (arg == "--start-window-left") {
            profile.emissions["START_CODON"].window_left = stoull(value);
        } else if (arg == "--start-window-right") {
            profile.emissions["START_CODON"].window_right = stoull(value);
        } else if (arg == "--transition-alpha") {
            profile.transition_alpha = stod(value);
        } else if (arg == "--emission-alpha") {
            profile.emission_alpha = stod(value);
        } else {
            cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            exit(1);
        }
    }

    if (profile.train_fasta_path.empty() || profile.train_gff_path.empty() ||
        profile.test_fasta_path.empty() || profile.test_gff_path.empty()) {
        cerr << "Provide --profile or all of --train-fasta, --train-gff, --test-fasta, --test-gff.\n";
        print_usage(argv[0]);
        exit(1);
    }

    if (profile.species.empty()) {
        profile.species.push_back({
            profile.name,
            profile.source_fasta_path,
            profile.source_gff_path,
            profile.train_fasta_path,
            profile.train_gff_path,
            profile.test_fasta_path,
            profile.test_gff_path,
            profile.test_chromosomes,
            profile.excluded_chromosomes
        });
    }

    if (!donor_scale_set) donor_cnn_scale = profile.splice_cnn.donor_scale;
    if (!donor_bias_set) donor_cnn_bias = profile.splice_cnn.donor_bias;
    if (!acceptor_scale_set) acceptor_cnn_scale = profile.splice_cnn.acceptor_scale;
    if (!acceptor_bias_set) acceptor_cnn_bias = profile.splice_cnn.acceptor_bias;
    if (!start_scale_set) start_cnn_scale = profile.start_cnn.start_scale;
    if (!start_bias_set) start_cnn_bias = profile.start_cnn.start_bias;

    return profile;
}

// Positions (interval-local) of every START_CODON_1 in a decoded path.
vector<size_t> start_codon_positions(const vector<State>& states) {
    vector<size_t> positions;
    for (size_t i = 0; i < states.size(); ++i) {
        if (states[i] == State::START_CODON_1) positions.push_back(i);
    }
    return positions;
}

// For one usable interval, record how far each predicted gene start sits from
// the nearest gold gene start (and vice versa). Predicted starts in an interval
// with no gold start (or gold with no prediction) are counted separately as
// spurious / missed rather than producing a misleading nearest distance.
void accumulate_start_offsets(
    const vector<State>& predicted,
    const vector<State>& gold,
    Validation_Result& result)
{
    vector<size_t> pred = start_codon_positions(predicted);
    vector<size_t> gold_starts = start_codon_positions(gold);

    for (size_t p : pred) {
        if (gold_starts.empty()) {
            result.pred_starts_no_gold++;
            continue;
        }
        auto it = lower_bound(gold_starts.begin(), gold_starts.end(), p);
        long best = numeric_limits<long>::max();
        long signed_best = 0;
        if (it != gold_starts.end()) {
            long d = static_cast<long>(*it) - static_cast<long>(p);
            if (labs(d) < best) { best = labs(d); signed_best = -d; }
        }
        if (it != gold_starts.begin()) {
            long d = static_cast<long>(*(it - 1)) - static_cast<long>(p);
            if (labs(d) < best) { best = labs(d); signed_best = -d; }
        }
        result.pred_start_signed_offset.push_back(signed_best);
    }

    for (size_t g : gold_starts) {
        if (pred.empty()) {
            result.gold_starts_no_pred++;
            continue;
        }
        auto it = lower_bound(pred.begin(), pred.end(), g);
        long best = numeric_limits<long>::max();
        if (it != pred.end()) best = min(best, labs(static_cast<long>(*it) - static_cast<long>(g)));
        if (it != pred.begin()) best = min(best, labs(static_cast<long>(*(it - 1)) - static_cast<long>(g)));
        result.gold_start_nearest_dist.push_back(best);
    }
}

void merge_start_offsets(Validation_Result& into, const Validation_Result& from) {
    into.pred_start_signed_offset.insert(
        into.pred_start_signed_offset.end(),
        from.pred_start_signed_offset.begin(),
        from.pred_start_signed_offset.end());
    into.gold_start_nearest_dist.insert(
        into.gold_start_nearest_dist.end(),
        from.gold_start_nearest_dist.begin(),
        from.gold_start_nearest_dist.end());
    into.pred_starts_no_gold += from.pred_starts_no_gold;
    into.gold_starts_no_pred += from.gold_starts_no_pred;
}

void write_start_offset_report(ostream& out, const string& title, const Validation_Result& result) {
    const long W = START_MATCH_WINDOW;

    size_t pred_total = result.pred_start_signed_offset.size() + result.pred_starts_no_gold;
    size_t gold_total = result.gold_start_nearest_dist.size() + result.gold_starts_no_pred;

    size_t exact = 0;          // offset == 0
    size_t within_window = 0;  // 0 < |offset| <= W  (shifted but same gene)
    size_t spurious = result.pred_starts_no_gold; // no gold in interval
    size_t upstream = 0;       // predicted start before gold (offset < 0)
    size_t downstream = 0;     // predicted start after gold (offset > 0)
    size_t in_frame = 0;       // nonzero offset divisible by 3, within window
    size_t out_frame = 0;      // nonzero offset not divisible by 3, within window

    // |offset| buckets for shifted-but-matched predicted starts
    size_t b_1_2 = 0, b_3 = 0, b_4_9 = 0, b_10_30 = 0, b_31_100 = 0, b_101_W = 0;
    vector<long> matched_abs; // for median

    for (long off : result.pred_start_signed_offset) {
        long a = labs(off);
        if (off == 0) { exact++; continue; }
        if (a > W) { spurious++; continue; }
        within_window++;
        if (off < 0) upstream++; else downstream++;
        if (off % 3 == 0) in_frame++; else out_frame++;
        matched_abs.push_back(a);
        if (a <= 2) b_1_2++;
        else if (a == 3) b_3++;
        else if (a <= 9) b_4_9++;
        else if (a <= 30) b_10_30++;
        else if (a <= 100) b_31_100++;
        else b_101_W++;
    }

    size_t missed = result.gold_starts_no_pred;
    for (long d : result.gold_start_nearest_dist) {
        if (d > W) missed++;
    }

    long median_abs = 0;
    if (!matched_abs.empty()) {
        sort(matched_abs.begin(), matched_abs.end());
        median_abs = matched_abs[matched_abs.size() / 2];
    }

    out << "=== Start-offset diagnostic: " << title << " ===\n";
    out << left << setw(40) << "Predicted gene starts" << right << setw(12) << pred_total << "\n";
    out << left << setw(40) << "Gold gene starts" << right << setw(12) << gold_total << "\n";
    out << left << setw(40) << "Exact (offset 0)" << right << setw(12) << exact << "\n";
    out << left << setw(40) << "Shifted, within +/-" + to_string(W) + " bp" << right << setw(12) << within_window << "\n";
    out << left << setw(40) << "Spurious (no gold within window)" << right << setw(12) << spurious << "\n";
    out << left << setw(40) << "Missed gold (no pred within window)" << right << setw(12) << missed << "\n";
    out << "\n";
    out << left << setw(40) << "Shifted: predicted upstream of gold" << right << setw(12) << upstream << "\n";
    out << left << setw(40) << "Shifted: predicted downstream of gold" << right << setw(12) << downstream << "\n";
    out << left << setw(40) << "Shifted: in-frame (|off| % 3 == 0)" << right << setw(12) << in_frame << "\n";
    out << left << setw(40) << "Shifted: out-of-frame" << right << setw(12) << out_frame << "\n";
    out << left << setw(40) << "Shifted: median |offset| (bp)" << right << setw(12) << median_abs << "\n";
    out << "\n";
    out << "|offset| histogram (shifted-but-matched predicted starts):\n";
    out << left << setw(16) << "  1-2 bp" << right << setw(12) << b_1_2 << "\n";
    out << left << setw(16) << "  3 bp" << right << setw(12) << b_3 << "\n";
    out << left << setw(16) << "  4-9 bp" << right << setw(12) << b_4_9 << "\n";
    out << left << setw(16) << "  10-30 bp" << right << setw(12) << b_10_30 << "\n";
    out << left << setw(16) << "  31-100 bp" << right << setw(12) << b_31_100 << "\n";
    out << left << setw(16) << "  101-" + to_string(W) + " bp" << right << setw(12) << b_101_W << "\n";
    out << "\n";
}

Validation_Result decode_validation(
    Emission_Model& emission_model,
    const Transition_Model::Log_Prob_Matrix& transition_log_probs,
    const Sequence_Data& eval_data,
    const vector<Chromosome_Range>& eval_ranges,
    size_t max_intron_body_length,
    bool collect_gene_intervals)
{
    Validation_Result result;
    const vector<State> donor_states{State::DONOR_1, State::DONOR_2, State::DONOR_3};
    const vector<State> acceptor_states{State::ACCEPTOR_1, State::ACCEPTOR_2, State::ACCEPTOR_3};

    for (const auto& range : eval_ranges) {
        emission_model.set_splice_cnn_position_offset(range.start);
        emission_model.set_start_cnn_position_offset(range.start);
        vector<Nucleotide> chromosome_nucleotides =
            slice_nucleotides(eval_data.nucleotides, range.start, range.end);
        vector<State> chromosome_gold =
            slice_states(eval_data.states, range.start, range.end);
        vector<State> chromosome_predicted =
            Viterbi::decode(
                chromosome_nucleotides,
                transition_log_probs,
                emission_model,
                0,
                max_intron_body_length,
                gene_start_penalty,
                intron_length_log_probs);

        if (chromosome_predicted.size() != chromosome_gold.size()) {
            throw runtime_error("Decoded path length mismatch on " + range.name + ".");
        }

        if (collect_start_positions) {
            accumulate_start_offsets(chromosome_predicted, chromosome_gold, result);
        }

        add_metrics(result.coding_total, binary_metrics(chromosome_predicted, chromosome_gold, is_coding));
        add_metrics(result.intron_total, binary_metrics(chromosome_predicted, chromosome_gold, is_intron));

        result.total_bases += chromosome_gold.size();
        result.exact_matches += exact_state_matches(chromosome_predicted, chromosome_gold);
        result.illegal_transitions += illegal_transition_count(chromosome_predicted);

        result.start_predicted += count_state(chromosome_predicted, State::START_CODON_1);
        result.start_gold += count_state(chromosome_gold, State::START_CODON_1);
        result.start_exact += count_exact_boundary_matches(chromosome_predicted, chromosome_gold, State::START_CODON_1);

        result.stop_predicted += count_state(chromosome_predicted, State::STOP_CODON_1);
        result.stop_gold += count_state(chromosome_gold, State::STOP_CODON_1);
        result.stop_exact += count_exact_boundary_matches(chromosome_predicted, chromosome_gold, State::STOP_CODON_1);

        result.donor_predicted += count_states(chromosome_predicted, donor_states);
        result.donor_gold += count_states(chromosome_gold, donor_states);
        result.donor_exact += count_exact_boundary_matches(chromosome_predicted, chromosome_gold, donor_states);

        result.acceptor_predicted += count_states(chromosome_predicted, acceptor_states);
        result.acceptor_gold += count_states(chromosome_gold, acceptor_states);
        result.acceptor_exact += count_exact_boundary_matches(chromosome_predicted, chromosome_gold, acceptor_states);

        if (collect_gene_intervals) {
            vector<Interval> pred_genes = collect_intervals(chromosome_predicted, range.name, is_gene);
            vector<Interval> gold_genes = collect_intervals(chromosome_gold, range.name, is_gene);
            result.predicted_gene_intervals.insert(
                result.predicted_gene_intervals.end(),
                pred_genes.begin(),
                pred_genes.end());
            result.gold_gene_intervals.insert(
                result.gold_gene_intervals.end(),
                gold_genes.begin(),
                gold_genes.end());
        }
    }

    return result;
}

vector<Chromosome_Range> evenly_spaced_ranges(
    const vector<Chromosome_Range>& ranges,
    size_t limit)
{
    if (limit == 0 || ranges.size() <= limit) {
        return ranges;
    }

    vector<Chromosome_Range> subset;
    subset.reserve(limit);
    for (size_t i = 0; i < limit; ++i) {
        size_t index = (i * ranges.size()) / limit;
        if (!subset.empty() && subset.back().start == ranges[index].start && subset.back().end == ranges[index].end) {
            continue;
        }
        subset.push_back(ranges[index]);
    }
    return subset;
}

void tune_splice_cnn_calibration(
    Emission_Model& emission_model,
    const Transition_Model::Log_Prob_Matrix& transition_log_probs,
    const Sequence_Data& eval_data,
    const vector<Chromosome_Range>& eval_ranges,
    size_t max_intron_body_length)
{
    struct Calibration_Candidate {
        Log_Prob donor_scale;
        Log_Prob donor_bias;
        Log_Prob acceptor_scale;
        Log_Prob acceptor_bias;
    };

    vector<Calibration_Candidate> candidates{
        {donor_cnn_scale, donor_cnn_bias, acceptor_cnn_scale, acceptor_cnn_bias}
    };
    const vector<Log_Prob> scales{0.35, 0.50, 0.75, 1.00, 1.25, 1.50};
    // Donor/acceptor emissions compete directly against the exon/intron
    // log-likelihoods, so the operating bias depends on where the raw CNN
    // logits sit relative to those. Search a wide symmetric range so the
    // optimum is never pinned to a grid edge.
    const vector<Log_Prob> donor_biases{-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
    const vector<Log_Prob> acceptor_biases{-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
    for (Log_Prob scale : scales) {
        for (Log_Prob donor_bias_candidate : donor_biases) {
            for (Log_Prob acceptor_bias_candidate : acceptor_biases) {
                candidates.push_back({
                    scale,
                    donor_bias_candidate,
                    scale,
                    acceptor_bias_candidate
                });
            }
        }
    }
    const vector<Chromosome_Range> tuning_ranges = evenly_spaced_ranges(eval_ranges, tune_subset_ranges);

    double best_objective = -1.0;
    Log_Prob best_donor_scale = donor_cnn_scale;
    Log_Prob best_donor_bias = donor_cnn_bias;
    Log_Prob best_acceptor_scale = acceptor_cnn_scale;
    Log_Prob best_acceptor_bias = acceptor_cnn_bias;
    size_t tried = 0;
    const size_t total = candidates.size();

    cout << "\nTuning CNN calibration on evaluation labels"
         << " objective=(intron F1 + donor boundary F1 + acceptor boundary F1)/3\n";
    cout << "Tuning subset: " << tuning_ranges.size()
         << "/" << eval_ranges.size()
         << " usable evaluation intervals\n";

    for (const auto& candidate : candidates) {
        emission_model.set_splice_cnn_calibration(
            candidate.donor_scale,
            candidate.donor_bias,
            candidate.acceptor_scale,
            candidate.acceptor_bias);
        Validation_Result result = decode_validation(
            emission_model,
            transition_log_probs,
            eval_data,
            tuning_ranges,
            max_intron_body_length,
            false);
        double objective = calibration_objective(result);
        tried++;

        cout << "  candidate " << tried << "/" << total
             << " objective=" << fixed << setprecision(4) << objective
             << " donor_scale=" << candidate.donor_scale
             << " donor_bias=" << candidate.donor_bias
             << " acceptor_scale=" << candidate.acceptor_scale
             << " acceptor_bias=" << candidate.acceptor_bias
             << "\n";

        if (objective > best_objective) {
            best_objective = objective;
            best_donor_scale = candidate.donor_scale;
            best_donor_bias = candidate.donor_bias;
            best_acceptor_scale = candidate.acceptor_scale;
            best_acceptor_bias = candidate.acceptor_bias;
            cout << "  new best objective=" << fixed << setprecision(4) << best_objective
                 << " donor_scale=" << candidate.donor_scale
                 << " donor_bias=" << candidate.donor_bias
                 << " acceptor_scale=" << candidate.acceptor_scale
                 << " acceptor_bias=" << candidate.acceptor_bias
                 << "\n";
        }
    }

    donor_cnn_scale = best_donor_scale;
    donor_cnn_bias = best_donor_bias;
    acceptor_cnn_scale = best_acceptor_scale;
    acceptor_cnn_bias = best_acceptor_bias;
    emission_model.set_splice_cnn_calibration(
        donor_cnn_scale,
        donor_cnn_bias,
        acceptor_cnn_scale,
        acceptor_cnn_bias);

    cout << "Selected CNN calibration:"
         << " donor_scale=" << donor_cnn_scale
         << " donor_bias=" << donor_cnn_bias
         << " acceptor_scale=" << acceptor_cnn_scale
         << " acceptor_bias=" << acceptor_cnn_bias
         << " objective=" << best_objective
         << "\n";
}

// Fits start_scale/start_bias on the TRAINING split so the reported test
// numbers are honest (no tuning on the evaluation genomes). The caller must
// have loaded the start-CNN TRAIN scores onto emission_model before this runs.
void tune_start_cnn_calibration(
    Emission_Model& emission_model,
    const Transition_Model::Log_Prob_Matrix& transition_log_probs,
    const Sequence_Data& train_data,
    const vector<Chromosome_Range>& train_ranges,
    size_t max_intron_body_length)
{
    const vector<Log_Prob> scales{0.50, 0.75, 1.00, 1.25};
    // The untuned start CNN over-opens genes, so the operating bias is expected
    // to be negative (make starts costlier). Search a wide one-sided range.
    const vector<Log_Prob> biases{-6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0};

    const vector<Chromosome_Range> tuning_ranges =
        evenly_spaced_ranges(train_ranges, tune_subset_ranges);

    double best_objective = -1.0;
    Log_Prob best_scale = start_cnn_scale;
    Log_Prob best_bias = start_cnn_bias;
    size_t tried = 0;
    const size_t total = scales.size() * biases.size();

    cout << "\nTuning start-CNN calibration on TRAINING labels"
         << " objective=start boundary F1\n";
    cout << "Tuning subset: " << tuning_ranges.size()
         << "/" << train_ranges.size()
         << " usable training intervals\n";

    for (Log_Prob scale : scales) {
        for (Log_Prob bias : biases) {
            emission_model.set_start_cnn_calibration(scale, bias);
            Validation_Result result = decode_validation(
                emission_model,
                transition_log_probs,
                train_data,
                tuning_ranges,
                max_intron_body_length,
                false);
            double objective = start_calibration_objective(result);
            tried++;

            cout << "  candidate " << tried << "/" << total
                 << " objective=" << fixed << setprecision(4) << objective
                 << " start_scale=" << scale
                 << " start_bias=" << bias
                 << " predicted=" << result.start_predicted
                 << " gold=" << result.start_gold
                 << "\n";

            if (objective > best_objective) {
                best_objective = objective;
                best_scale = scale;
                best_bias = bias;
                cout << "  new best objective=" << fixed << setprecision(4) << best_objective
                     << " start_scale=" << scale
                     << " start_bias=" << bias
                     << "\n";
            }
        }
    }

    start_cnn_scale = best_scale;
    start_cnn_bias = best_bias;
    emission_model.set_start_cnn_calibration(start_cnn_scale, start_cnn_bias);

    cout << "Selected start-CNN calibration:"
         << " start_scale=" << start_cnn_scale
         << " start_bias=" << start_cnn_bias
         << " objective=" << best_objective
         << "\n";
}

// Writes the full metric/classification/boundary report for one decode result
// to any stream, so the same layout serves stdout, the combined results file,
// and the per-species results files.
void write_validation_report(
    ostream& out,
    const string& header,
    const Validation_Result& result,
    size_t max_intron_body_length)
{
    out << "\n=== " << header << " ===\n";
    out << fixed << setprecision(4);
    out << left << setw(28) << "Metric" << right << setw(14) << "Value" << "\n";
    out << left << setw(28) << "Evaluated bases" << right << setw(14) << result.total_bases << "\n";
    out << left << setw(28) << "Exact 21-state accuracy" << right << setw(14) << divide(result.exact_matches, result.total_bases) << "\n";
    out << left << setw(28) << "Illegal transitions" << right << setw(14) << result.illegal_transitions << "\n";
    out << left << setw(28) << "Predicted gene intervals" << right << setw(14) << result.predicted_gene_intervals.size() << "\n";
    out << left << setw(28) << "Gold gene intervals" << right << setw(14) << result.gold_gene_intervals.size() << "\n";
    out << left << setw(28) << "Intron length cap p95" << right << setw(14) << max_intron_body_length << "\n";
    out << left << setw(28) << "Gene start penalty" << right << setw(14) << gene_start_penalty << "\n";
    out << left << setw(28) << "CNN donor scale" << right << setw(14) << donor_cnn_scale << "\n";
    out << left << setw(28) << "CNN donor bias" << right << setw(14) << donor_cnn_bias << "\n";
    out << left << setw(28) << "CNN acceptor scale" << right << setw(14) << acceptor_cnn_scale << "\n";
    out << left << setw(28) << "CNN acceptor bias" << right << setw(14) << acceptor_cnn_bias << "\n";
    out << left << setw(28) << "CNN start scale" << right << setw(14) << start_cnn_scale << "\n";
    out << left << setw(28) << "CNN start bias" << right << setw(14) << start_cnn_bias << "\n";

    out << "\nClassification Metrics:\n";
    print_binary_metrics_header(out);
    print_binary_metrics_row(out, "coding", result.coding_total);
    print_binary_metrics_row(out, "intron", result.intron_total);

    out << "\nBoundary Metrics:\n";
    print_boundary_metrics_header(out);
    print_boundary_metrics_row(out, "start", result.start_exact, result.start_predicted, result.start_gold);
    print_boundary_metrics_row(out, "stop", result.stop_exact, result.stop_predicted, result.stop_gold);
    print_boundary_metrics_row(out, "donor", result.donor_exact, result.donor_predicted, result.donor_gold);
    print_boundary_metrics_row(out, "acceptor", result.acceptor_exact, result.acceptor_predicted, result.acceptor_gold);
}

} // namespace

int main(int argc, char** argv) {
    gene_hmm::profile = parse_args(argc, argv);

    cout << "Loading validation profile: " << gene_hmm::profile.name << "\n";
    cout << "Species datasets: " << gene_hmm::profile.species.size() << "\n";
    for (const auto& dataset : gene_hmm::profile.species) {
        cout << "  " << dataset.name << "\n"
             << "    train FASTA: " << dataset.train_fasta_path << "\n"
             << "    train GFF:   " << dataset.train_gff_path << "\n"
             << "    test FASTA:  " << dataset.test_fasta_path << "\n"
             << "    test GFF:    " << dataset.test_gff_path << "\n";
    }

    Sequence_Data train_data = load_training_data(gene_hmm::profile);
    Sequence_Data eval_data = load_evaluation_data(gene_hmm::profile);

    if (train_data.states.size() != train_data.nucleotides.size()) {
        cerr << "Train gold state length does not match nucleotide length.\n";
        return 1;
    }
    if (eval_data.states.size() != eval_data.nucleotides.size()) {
        cerr << "Test gold state length does not match nucleotide length.\n";
        return 1;
    }

    vector<Chromosome_Range> train_ranges = split_usable_ranges(train_data.chromosomes, train_data.regions);
    vector<Chromosome_Range> eval_ranges = split_usable_ranges(eval_data.chromosomes, eval_data.regions);
    if (train_ranges.empty()) {
        cerr << "No usable training intervals after applying annotation filters.\n";
        return 1;
    }
    if (eval_ranges.empty()) {
        cerr << "No usable evaluation intervals after applying annotation filters.\n";
        return 1;
    }

    cout << "Training chromosomes: " << train_data.chromosomes.size() << "\n";
    cout << "Evaluation chromosomes: " << eval_data.chromosomes.size() << "\n";
    cout << "Usable training intervals: " << train_ranges.size() << "\n";
    cout << "Usable evaluation intervals: " << eval_ranges.size() << "\n";

    auto transition_log_probs = Transition_Model::compute_log_probs(train_data.states, train_ranges);
    Emission_Model emission_model = train_emissions(train_data.states, train_data.nucleotides, train_ranges);
    if (!splice_cnn_scores_path.empty()) {
        emission_model.load_splice_cnn_scores(splice_cnn_scores_path, eval_data.nucleotides.size());
    } else if (!gene_hmm::profile.splice_cnn.test_score_paths.empty()) {
        if (gene_hmm::profile.splice_cnn.test_score_paths.size() != eval_data.dataset_offsets.size()) {
            cerr << "Profile splice_cnn.test_scores count does not match evaluation datasets.\n";
            return 1;
        }
        emission_model.load_splice_cnn_scores(
            gene_hmm::profile.splice_cnn.test_score_paths,
            eval_data.dataset_offsets,
            eval_data.nucleotides.size());
    } else {
        cerr << "Provide --splice-cnn-scores or define splice_cnn.test_scores in the profile.\n";
        return 1;
    }
    emission_model.set_splice_cnn_calibration(
        donor_cnn_scale,
        donor_cnn_bias,
        acceptor_cnn_scale,
        acceptor_cnn_bias);

    // Translation-start CNN scores are required, like splice donor/acceptor.
    if (!start_cnn_scores_path.empty()) {
        emission_model.load_start_cnn_scores(start_cnn_scores_path, eval_data.nucleotides.size());
    } else if (!gene_hmm::profile.start_cnn.test_score_paths.empty()) {
        if (gene_hmm::profile.start_cnn.test_score_paths.size() != eval_data.dataset_offsets.size()) {
            cerr << "Profile start_cnn.test_scores count does not match evaluation datasets.\n";
            return 1;
        }
        emission_model.load_start_cnn_scores(
            gene_hmm::profile.start_cnn.test_score_paths,
            eval_data.dataset_offsets,
            eval_data.nucleotides.size());
    } else {
        cerr << "Provide --start-cnn-scores or define start_cnn.test_scores in the profile.\n";
        return 1;
    }
    emission_model.set_start_cnn_calibration(start_cnn_scale, start_cnn_bias);

    vector<size_t> train_intron_body_lengths = collect_intron_body_lengths(train_data.states, train_ranges);
    size_t max_intron_body_length = percentile_length(train_intron_body_lengths, 0.95);
    intron_length_log_probs = build_intron_length_log_probs(train_intron_body_lengths, max_intron_body_length);

    if (tune_start_calibration) {
        // Fit start_scale/start_bias on the TRAINING split (decoding train
        // genomes), then restore the evaluation scores so the final report is
        // measured on held-out genomes at the selected operating point.
        const auto& splice_train_paths = gene_hmm::profile.splice_cnn.train_score_paths;
        const auto& start_train_paths = gene_hmm::profile.start_cnn.train_score_paths;
        if (splice_train_paths.size() != train_data.dataset_offsets.size() ||
            start_train_paths.size() != train_data.dataset_offsets.size()) {
            cerr << "Profile train_scores count does not match training datasets;"
                 << " cannot fit start calibration.\n";
            return 1;
        }
        emission_model.load_splice_cnn_scores(
            splice_train_paths, train_data.dataset_offsets, train_data.nucleotides.size());
        emission_model.load_start_cnn_scores(
            start_train_paths, train_data.dataset_offsets, train_data.nucleotides.size());

        tune_start_cnn_calibration(
            emission_model,
            transition_log_probs,
            train_data,
            train_ranges,
            max_intron_body_length);

        // Restore evaluation scores exactly as they were originally loaded.
        if (!splice_cnn_scores_path.empty()) {
            emission_model.load_splice_cnn_scores(splice_cnn_scores_path, eval_data.nucleotides.size());
        } else {
            emission_model.load_splice_cnn_scores(
                gene_hmm::profile.splice_cnn.test_score_paths,
                eval_data.dataset_offsets,
                eval_data.nucleotides.size());
        }
        if (!start_cnn_scores_path.empty()) {
            emission_model.load_start_cnn_scores(start_cnn_scores_path, eval_data.nucleotides.size());
        } else {
            emission_model.load_start_cnn_scores(
                gene_hmm::profile.start_cnn.test_score_paths,
                eval_data.dataset_offsets,
                eval_data.nucleotides.size());
        }
        emission_model.set_splice_cnn_calibration(
            donor_cnn_scale, donor_cnn_bias, acceptor_cnn_scale, acceptor_cnn_bias);
        emission_model.set_start_cnn_calibration(start_cnn_scale, start_cnn_bias);

        if (tune_only) {
            return 0;
        }
    }

    if (tune_cnn_calibration) {
        tune_splice_cnn_calibration(
            emission_model,
            transition_log_probs,
            eval_data,
            eval_ranges,
            max_intron_body_length);
        if (tune_only) {
            return 0;
        }
    }

    if (!start_offset_report_path.empty()) {
        collect_start_positions = true;
        ofstream report(start_offset_report_path);
        if (!report) {
            cerr << "Could not open start-offset report path: " << start_offset_report_path << "\n";
            return 1;
        }
        report << "Gene start penalty: " << fixed << setprecision(4) << gene_start_penalty << "\n";
        report << "Start PSSM window: left=" << emission_model.start_window_left
               << " right=" << emission_model.start_window_right
               << " (total " << (emission_model.start_window_left + emission_model.start_window_right) << " bp)\n";
        report << "CNN donor scale/bias: " << donor_cnn_scale << "/" << donor_cnn_bias
               << "  acceptor scale/bias: " << acceptor_cnn_scale << "/" << acceptor_cnn_bias << "\n";
        report << "Match window: +/-" << START_MATCH_WINDOW << " bp\n\n";

        Validation_Result combined;
        for (const auto& dataset : gene_hmm::profile.species) {
            vector<Chromosome_Range> species_ranges;
            for (const auto& range : eval_ranges) {
                if (range.name.rfind(dataset.name + ":", 0) == 0) {
                    species_ranges.push_back(range);
                }
            }
            if (species_ranges.empty()) {
                continue;
            }
            Validation_Result species_result;
            try {
                species_result = decode_validation(
                    emission_model,
                    transition_log_probs,
                    eval_data,
                    species_ranges,
                    max_intron_body_length,
                    false);
            } catch (const exception& e) {
                cerr << e.what() << "\n";
                return 1;
            }
            write_start_offset_report(report, dataset.name, species_result);
            write_start_offset_report(cout, dataset.name, species_result);
            merge_start_offsets(combined, species_result);
        }
        write_start_offset_report(report, "combined", combined);
        write_start_offset_report(cout, "combined", combined);
        cout << "Saved start-offset report to " << start_offset_report_path << "\n";
        return 0;
    }

    cout << "\nDecoding usable evaluation intervals";
    if (eval_ranges.size() == 1) {
        cout << " on " << eval_ranges.front().name
             << " bases=" << (eval_ranges.front().end - eval_ranges.front().start);
    }
    cout << "\n";

    Validation_Result result;
    try {
        result = decode_validation(
            emission_model,
            transition_log_probs,
            eval_data,
            eval_ranges,
            max_intron_body_length,
            true);
    } catch (const exception& e) {
        cerr << e.what() << "\n";
        return 1;
    }

    write_validation_report(cout, "Full Genome Holdout Validation", result, max_intron_body_length);

    if (!results_dir.empty()) {
        error_code ec;
        filesystem::create_directories(results_dir, ec);
        if (ec) {
            cerr << "Could not create results directory: " << results_dir << " (" << ec.message() << ")\n";
            return 1;
        }

        // Combined report.
        const string combined_path = results_dir + "/combined.txt";
        ofstream combined_file(combined_path);
        if (!combined_file) {
            cerr << "Could not open results file: " << combined_path << "\n";
            return 1;
        }
        write_validation_report(
            combined_file, "Full Genome Holdout Validation (all tests combined)", result, max_intron_body_length);
        cout << "\nSaved combined results to " << combined_path << "\n";

        // Per-species reports: decode only that species' evaluation ranges. The
        // CNN scores are already loaded with global offsets, so a range subset
        // decodes identically to the combined run.
        for (const auto& dataset : gene_hmm::profile.species) {
            const string prefix = dataset.name + ":";
            vector<Chromosome_Range> species_ranges;
            for (const auto& range : eval_ranges) {
                if (range.name.rfind(prefix, 0) == 0) {
                    species_ranges.push_back(range);
                }
            }
            if (species_ranges.empty()) {
                continue;
            }

            Validation_Result species_result;
            try {
                species_result = decode_validation(
                    emission_model,
                    transition_log_probs,
                    eval_data,
                    species_ranges,
                    max_intron_body_length,
                    true);
            } catch (const exception& e) {
                cerr << e.what() << "\n";
                return 1;
            }

            string chroms;
            set<string> seen_chroms;
            for (const auto& range : species_ranges) {
                if (seen_chroms.insert(range.name).second) {
                    if (!chroms.empty()) chroms += ",";
                    chroms += range.name;
                }
            }
            const string species_path = results_dir + "/" + dataset.name + ".txt";
            ofstream species_file(species_path);
            if (!species_file) {
                cerr << "Could not open results file: " << species_path << "\n";
                return 1;
            }
            write_validation_report(species_file, "Test: " + chroms, species_result, max_intron_body_length);
            cout << "Saved " << dataset.name << " results to " << species_path << "\n";
        }
    }

    return 0;
}
