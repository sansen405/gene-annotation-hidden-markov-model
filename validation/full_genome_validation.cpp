#include "../src/decoding/Viterbi.hpp"
#include "../src/genome_profiles/Genome_Profile.hpp"
#include "../src/model/Emission_Model.hpp"
#include "../src/model/Transition_Model.hpp"
#include "../src/parsers/FNA_Parser.hpp"
#include "../src/parsers/GFF_Parser.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

using namespace gene_hmm;
using namespace std;

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

void print_binary_metrics(const string& label, const Binary_Metrics& metrics) {
    double precision = divide(metrics.tp, metrics.tp + metrics.fp);
    double recall = divide(metrics.tp, metrics.tp + metrics.fn);
    double accuracy = divide(metrics.tp + metrics.tn, metrics.tp + metrics.fp + metrics.fn + metrics.tn);

    cout << fixed << setprecision(4);
    cout << label
         << " precision=" << precision
         << " recall=" << recall
         << " f1=" << f1(precision, recall)
         << " accuracy=" << accuracy
         << " tp=" << metrics.tp
         << " fp=" << metrics.fp
         << " fn=" << metrics.fn
         << " tn=" << metrics.tn << "\n";
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

    model.donor_lp = Emission_Model::compute_pssm_log_probs(
        Emission_Model::count_pssm_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::DONOR_1, State::DONOR_2, State::DONOR_3},
            model.donor_window_left,
            model.donor_window_right));

    model.acceptor_lp = Emission_Model::compute_pssm_log_probs(
        Emission_Model::count_pssm_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::ACCEPTOR_1, State::ACCEPTOR_2, State::ACCEPTOR_3},
            model.acceptor_window_left,
            model.acceptor_window_right));

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

set<string> make_set(const vector<string>& values) {
    return set<string>(values.begin(), values.end());
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

Genome_Profile yeast_profile() {
    Genome_Profile profile;
    profile.name = "yeast";
    profile.fasta_path = "genome_data/yeast_data/GCF_000146045.2_R64_genomic.fna";
    profile.gff_path = "genome_data/yeast_data/genomic.gff";
    profile.test_chromosomes = {"NC_001148.4"};
    profile.exclude_chromosomes = {"NC_001224.1"};
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

} // namespace

namespace gene_hmm {
    Genome_Profile profile;
}

int main(int argc, char** argv) {
    if (argc > 1) {
        cerr << "This standalone validator currently uses the built-in yeast profile defaults.\n";
        cerr << "Usage: " << argv[0] << "\n";
        return 1;
    }

    gene_hmm::profile = yeast_profile();

    cout << "Loading built-in yeast validation profile\n";
    cout << "FASTA: " << gene_hmm::profile.fasta_path << "\n";
    cout << "GFF: " << gene_hmm::profile.gff_path << "\n";

    vector<Nucleotide> nucleotides = FNA_Parser::parse_sequence(gene_hmm::profile.fasta_path);
    vector<Chromosome_Range> ranges = FNA_Parser::get_chromosome_ranges(gene_hmm::profile.fasta_path);

    string gff_path = gene_hmm::profile.gff_path;
    string fasta_path = gene_hmm::profile.fasta_path;
    vector<int> regions = GFF_Parser::parse_regions(gff_path, fasta_path);
    vector<State> gold_states = GFF_Parser::parse_states(regions);

    if (gold_states.size() != nucleotides.size()) {
        cerr << "Gold state length does not match nucleotide length.\n";
        return 1;
    }

    set<string> test_chromosomes = make_set(gene_hmm::profile.test_chromosomes);
    set<string> excluded_chromosomes = make_set(gene_hmm::profile.exclude_chromosomes);
    vector<Chromosome_Range> train_ranges;
    vector<Chromosome_Range> eval_ranges;

    for (const auto& range : ranges) {
        if (excluded_chromosomes.count(range.name)) {
            continue;
        }
        if (test_chromosomes.count(range.name)) {
            eval_ranges.push_back(range);
        } else {
            train_ranges.push_back(range);
        }
    }

    if (eval_ranges.empty()) {
        cerr << "No test chromosomes found in the profile/FASTA.\n";
        return 1;
    }
    if (train_ranges.empty()) {
        cerr << "No training chromosomes found after applying profile split.\n";
        return 1;
    }

    cout << "Training chromosomes: " << train_ranges.size() << "\n";
    cout << "Evaluation chromosomes: " << eval_ranges.size() << "\n";

    auto transition_log_probs = Transition_Model::compute_log_probs(gold_states, train_ranges);
    Emission_Model emission_model = train_emissions(gold_states, nucleotides, train_ranges);

    cout << "\nTraining diagnostics:\n";
    print_state_family_counts("gold all chromosomes", gold_states);
    print_log_prob("log P(EXON_FRAME_3 -> EXON_FRAME_1)", transition_log_probs[idx(State::EXON_FRAME_3)][idx(State::EXON_FRAME_1)]);
    print_log_prob("log P(EXON_FRAME_3 -> DONOR_1)", transition_log_probs[idx(State::EXON_FRAME_3)][idx(State::DONOR_1)]);
    print_log_prob("log P(EXON_FRAME_3 -> STOP_CODON_1)", transition_log_probs[idx(State::EXON_FRAME_3)][idx(State::STOP_CODON_1)]);
    print_log_prob("log P(INTRON_1 -> INTRON_1)", transition_log_probs[idx(State::INTRON_1)][idx(State::INTRON_1)]);
    print_log_prob("log P(INTRON_1 -> ACCEPTOR_1)", transition_log_probs[idx(State::INTRON_1)][idx(State::ACCEPTOR_1)]);

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
    vector<Interval> false_positive_intervals;
    vector<Interval> false_negative_intervals;
    vector<Interval> predicted_gene_intervals;
    vector<Interval> gold_gene_intervals;

    for (const auto& range : eval_ranges) {
        cout << "\nDecoding " << range.name
             << " bases=" << (range.end - range.start) << "\n";

        vector<Nucleotide> chromosome_nucleotides =
            slice_nucleotides(nucleotides, range.start, range.end);
        vector<State> chromosome_gold =
            slice_states(gold_states, range.start, range.end);
        vector<State> chromosome_predicted =
            Viterbi::decode(chromosome_nucleotides, transition_log_probs, emission_model);

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

        vector<Interval> fp = collect_error_intervals(chromosome_predicted, chromosome_gold, range.name, true);
        vector<Interval> fn = collect_error_intervals(chromosome_predicted, chromosome_gold, range.name, false);
        vector<Interval> pred_genes = collect_intervals(chromosome_predicted, range.name, is_gene);
        vector<Interval> gold_genes = collect_intervals(chromosome_gold, range.name, is_gene);

        print_state_family_counts("gold " + range.name, chromosome_gold);
        print_state_family_counts("predicted " + range.name, chromosome_predicted);
        print_average_emission(
            "gold exon emission magnitude " + range.name,
            emission_model,
            chromosome_nucleotides,
            chromosome_gold,
            {State::EXON_FRAME_1, State::EXON_FRAME_2, State::EXON_FRAME_3});
        print_average_emission(
            "gold donor PSSM emission magnitude " + range.name,
            emission_model,
            chromosome_nucleotides,
            chromosome_gold,
            {State::DONOR_1, State::DONOR_2, State::DONOR_3});
        print_average_emission(
            "gold acceptor PSSM emission magnitude " + range.name,
            emission_model,
            chromosome_nucleotides,
            chromosome_gold,
            {State::ACCEPTOR_1, State::ACCEPTOR_2, State::ACCEPTOR_3});

        false_positive_intervals.insert(false_positive_intervals.end(), fp.begin(), fp.end());
        false_negative_intervals.insert(false_negative_intervals.end(), fn.begin(), fn.end());
        predicted_gene_intervals.insert(predicted_gene_intervals.end(), pred_genes.begin(), pred_genes.end());
        gold_gene_intervals.insert(gold_gene_intervals.end(), gold_genes.begin(), gold_genes.end());
    }

    cout << "\n=== Full Genome Holdout Validation ===\n";
    cout << "evaluated bases=" << total_bases << "\n";
    cout << "exact 21-state accuracy=" << divide(exact_matches, total_bases) << "\n";
    print_binary_metrics("coding", coding_total);
    print_binary_metrics("intron", intron_total);

    cout << "\nBoundary metrics:\n";
    cout << "start exact_matches=" << start_exact
         << " predicted=" << start_predicted
         << " gold=" << start_gold
         << " precision=" << divide(start_exact, start_predicted)
         << " recall=" << divide(start_exact, start_gold) << "\n";
    cout << "stop exact_matches=" << stop_exact
         << " predicted=" << stop_predicted
         << " gold=" << stop_gold
         << " precision=" << divide(stop_exact, stop_predicted)
         << " recall=" << divide(stop_exact, stop_gold) << "\n";

    cout << "\nStructure:\n";
    cout << "illegal predicted transitions=" << illegal_transitions << "\n";
    print_intervals("predicted gene intervals", predicted_gene_intervals, 5);
    print_intervals("gold gene intervals", gold_gene_intervals, 5);

    cout << "\nError examples:\n";
    print_intervals("false-positive gene-state runs", false_positive_intervals, 10);
    print_intervals("false-negative gene-state runs", false_negative_intervals, 10);

    return 0;
}
