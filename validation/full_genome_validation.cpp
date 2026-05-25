#include "../src/decoding/Viterbi.hpp"
#include "../src/genome_profiles/Genome_Profile.hpp"
#include "../src/model/Emission_Model.hpp"
#include "../src/model/Transition_Model.hpp"
#include "../src/parsers/FNA_Parser.hpp"
#include "../src/parsers/GFF_Parser.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

using namespace gene_hmm;
using namespace std;

Log_Prob gene_start_penalty = 1.0;

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

void print_binary_metrics_row(const string& label, const Binary_Metrics& metrics) {
    double precision = divide(metrics.tp, metrics.tp + metrics.fp);
    double recall = divide(metrics.tp, metrics.tp + metrics.fn);
    double accuracy = divide(metrics.tp + metrics.tn, metrics.tp + metrics.fp + metrics.fn + metrics.tn);

    cout << fixed << setprecision(4);
    cout << left << setw(12) << label
         << right << setw(12) << precision
         << setw(12) << recall
         << setw(12) << f1(precision, recall)
         << setw(12) << accuracy
         << setw(10) << metrics.tp
         << setw(10) << metrics.fp
         << setw(10) << metrics.fn
         << setw(10) << metrics.tn << "\n";
}

void print_binary_metrics_header() {
    cout << left << setw(12) << "Label"
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
    const string& label,
    size_t exact,
    size_t predicted,
    size_t gold)
{
    cout << fixed << setprecision(4);
    cout << left << setw(12) << label
         << right << setw(12) << divide(exact, predicted)
         << setw(12) << divide(exact, gold)
         << setw(12) << exact
         << setw(12) << predicted
         << setw(12) << gold << "\n";
}

void print_boundary_metrics_header() {
    cout << left << setw(12) << "Boundary"
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

Genome_Profile default_profile() {
    return Genome_Profile::load("src/genome_profiles/yeast.json");
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
        return default_profile();
    }

    Genome_Profile profile = manual_profile();

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            exit(0);
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

    return profile;
}

} // namespace

int main(int argc, char** argv) {
    gene_hmm::profile = parse_args(argc, argv);

    cout << "Loading validation profile: " << gene_hmm::profile.name << "\n";
    cout << "Train FASTA: " << gene_hmm::profile.train_fasta_path << "\n";
    cout << "Train GFF:   " << gene_hmm::profile.train_gff_path << "\n";
    cout << "Test FASTA:  " << gene_hmm::profile.test_fasta_path << "\n";
    cout << "Test GFF:    " << gene_hmm::profile.test_gff_path << "\n";

    vector<Nucleotide> train_nucleotides =
        FNA_Parser::parse_sequence(gene_hmm::profile.train_fasta_path);
    vector<Chromosome_Range> train_chromosomes =
        FNA_Parser::get_chromosome_ranges(gene_hmm::profile.train_fasta_path);
    string train_gff_path = gene_hmm::profile.train_gff_path;
    string train_fasta_path = gene_hmm::profile.train_fasta_path;
    vector<int> train_regions = GFF_Parser::parse_regions(train_gff_path, train_fasta_path);
    vector<State> train_gold_states = GFF_Parser::parse_states(train_regions);

    vector<Nucleotide> eval_nucleotides =
        FNA_Parser::parse_sequence(gene_hmm::profile.test_fasta_path);
    vector<Chromosome_Range> eval_chromosomes =
        FNA_Parser::get_chromosome_ranges(gene_hmm::profile.test_fasta_path);
    string eval_gff_path = gene_hmm::profile.test_gff_path;
    string eval_fasta_path = gene_hmm::profile.test_fasta_path;
    vector<int> eval_regions = GFF_Parser::parse_regions(eval_gff_path, eval_fasta_path);
    vector<State> eval_gold_states = GFF_Parser::parse_states(eval_regions);

    if (train_gold_states.size() != train_nucleotides.size()) {
        cerr << "Train gold state length does not match nucleotide length.\n";
        return 1;
    }
    if (eval_gold_states.size() != eval_nucleotides.size()) {
        cerr << "Test gold state length does not match nucleotide length.\n";
        return 1;
    }

    vector<Chromosome_Range> train_ranges = split_usable_ranges(train_chromosomes, train_regions);
    vector<Chromosome_Range> eval_ranges = split_usable_ranges(eval_chromosomes, eval_regions);
    if (train_ranges.empty()) {
        cerr << "No usable training intervals after applying annotation filters.\n";
        return 1;
    }
    if (eval_ranges.empty()) {
        cerr << "No usable evaluation intervals after applying annotation filters.\n";
        return 1;
    }

    cout << "Training chromosomes: " << train_chromosomes.size() << "\n";
    cout << "Evaluation chromosomes: " << eval_chromosomes.size() << "\n";
    cout << "Usable training intervals: " << train_ranges.size() << "\n";
    cout << "Usable evaluation intervals: " << eval_ranges.size() << "\n";

    auto transition_log_probs = Transition_Model::compute_log_probs(train_gold_states, train_ranges);
    Emission_Model emission_model = train_emissions(train_gold_states, train_nucleotides, train_ranges);
    vector<size_t> train_intron_body_lengths = collect_intron_body_lengths(train_gold_states, train_ranges);
    size_t max_intron_body_length = percentile_length(train_intron_body_lengths, 0.95);

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
    vector<Interval> predicted_gene_intervals;
    vector<Interval> gold_gene_intervals;

    cout << "\nDecoding usable evaluation intervals";
    if (eval_ranges.size() == 1) {
        cout << " on " << eval_ranges.front().name
             << " bases=" << (eval_ranges.front().end - eval_ranges.front().start);
    }
    cout << "\n";

    for (const auto& range : eval_ranges) {
        vector<Nucleotide> chromosome_nucleotides =
            slice_nucleotides(eval_nucleotides, range.start, range.end);
        vector<State> chromosome_gold =
            slice_states(eval_gold_states, range.start, range.end);
        vector<State> chromosome_predicted =
            Viterbi::decode(
                chromosome_nucleotides,
                transition_log_probs,
                emission_model,
                0,
                max_intron_body_length,
                gene_start_penalty);

        if (chromosome_predicted.size() != chromosome_gold.size()) {
            cerr << "Decoded path length mismatch on " << range.name << ".\n";
            return 1;
        }

        add_metrics(coding_total, binary_metrics(chromosome_predicted, chromosome_gold, is_coding));
        add_metrics(intron_total, binary_metrics(chromosome_predicted, chromosome_gold, is_intron));

        total_bases += chromosome_gold.size();
        exact_matches += exact_state_matches(chromosome_predicted, chromosome_gold);
        illegal_transitions += illegal_transition_count(chromosome_predicted);

        start_predicted += count_state(chromosome_predicted, State::START_CODON_1);
        start_gold += count_state(chromosome_gold, State::START_CODON_1);
        start_exact += count_exact_boundary_matches(chromosome_predicted, chromosome_gold, State::START_CODON_1);

        stop_predicted += count_state(chromosome_predicted, State::STOP_CODON_1);
        stop_gold += count_state(chromosome_gold, State::STOP_CODON_1);
        stop_exact += count_exact_boundary_matches(chromosome_predicted, chromosome_gold, State::STOP_CODON_1);

        vector<Interval> pred_genes = collect_intervals(chromosome_predicted, range.name, is_gene);
        vector<Interval> gold_genes = collect_intervals(chromosome_gold, range.name, is_gene);

        predicted_gene_intervals.insert(predicted_gene_intervals.end(), pred_genes.begin(), pred_genes.end());
        gold_gene_intervals.insert(gold_gene_intervals.end(), gold_genes.begin(), gold_genes.end());
    }

    cout << "\n=== Full Genome Holdout Validation ===\n";
    cout << fixed << setprecision(4);
    cout << left << setw(28) << "Metric" << right << setw(14) << "Value" << "\n";
    cout << left << setw(28) << "Evaluated bases" << right << setw(14) << total_bases << "\n";
    cout << left << setw(28) << "Exact 21-state accuracy" << right << setw(14) << divide(exact_matches, total_bases) << "\n";
    cout << left << setw(28) << "Illegal transitions" << right << setw(14) << illegal_transitions << "\n";
    cout << left << setw(28) << "Predicted gene intervals" << right << setw(14) << predicted_gene_intervals.size() << "\n";
    cout << left << setw(28) << "Gold gene intervals" << right << setw(14) << gold_gene_intervals.size() << "\n";
    cout << left << setw(28) << "Intron length cap p95" << right << setw(14) << max_intron_body_length << "\n";
    cout << left << setw(28) << "Gene start penalty" << right << setw(14) << gene_start_penalty << "\n";

    cout << "\nClassification Metrics:\n";
    print_binary_metrics_header();
    print_binary_metrics_row("coding", coding_total);
    print_binary_metrics_row("intron", intron_total);

    cout << "\nBoundary Metrics:\n";
    print_boundary_metrics_header();
    print_boundary_metrics_row("start", start_exact, start_predicted, start_gold);
    print_boundary_metrics_row("stop", stop_exact, stop_predicted, stop_gold);

    return 0;
}
