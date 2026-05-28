#pragma once

#include "../model/transition/Transition_Model.hpp"
#include "../genome_profiles/Genome_Profile.hpp"
#include "Test_Utils.hpp"
#include <cmath>
#include <string>
#include <vector>

namespace gene_hmm {

    using namespace std;

    static bool approx_equal(Log_Prob actual, Log_Prob expected, double eps = 1e-12) {
        return fabs(actual - expected) < eps;
    }

    static vector<State> transition_fixture_states() {
        return {
            State::INTERGENIC,
            State::START_CODON_1,
            State::START_CODON_2,
            State::INTERGENIC,
            State::INTERGENIC,
            State::START_CODON_1
        };
    }

    static vector<Chromosome_Range> transition_fixture_ranges() {
        return {
            {"chr1", 0, 3},
            {"chr2", 3, 6}
        };
    }

    static void test_transition_bigram_counts() {
        cout << "\n[TEST 1] Bigram counts include START and END per chromosome\n";

        auto states = transition_fixture_states();
        auto ranges = transition_fixture_ranges();
        auto counts = Transition_Model::count_bigrams(states, ranges);

        CHECK("START -> INTERGENIC counted for each chromosome", counts[idx(State::START)][idx(State::INTERGENIC)] == 2);
        CHECK("INTERGENIC -> START_CODON_1 counted inside chromosomes", counts[idx(State::INTERGENIC)][idx(State::START_CODON_1)] == 2);
        CHECK("START_CODON_1 -> START_CODON_2 counted once", counts[idx(State::START_CODON_1)][idx(State::START_CODON_2)] == 1);
        CHECK("START_CODON_2 -> END counted at chr1 boundary", counts[idx(State::START_CODON_2)][idx(State::END)] == 1);
        CHECK("START_CODON_1 -> END counted at chr2 boundary", counts[idx(State::START_CODON_1)][idx(State::END)] == 1);
    }

    static void test_transition_no_cross_chromosome_bigrams() {
        cout << "\n[TEST 2] Bigram counts do not cross chromosome boundaries\n";

        auto states = transition_fixture_states();
        auto ranges = transition_fixture_ranges();
        auto counts = Transition_Model::count_bigrams(states, ranges);

        CHECK("last chr1 state is not counted to first chr2 state",
              counts[idx(State::START_CODON_2)][idx(State::INTERGENIC)] == 0);
    }

    static void test_transition_outgoing_counts() {
        cout << "\n[TEST 3] Outgoing row counts match observed source states\n";

        auto states = transition_fixture_states();
        auto ranges = transition_fixture_ranges();
        auto outgoing = Transition_Model::count_outgoing(states, ranges);

        CHECK("START has one outgoing transition per chromosome", outgoing[idx(State::START)] == 2);
        CHECK("INTERGENIC outgoing count includes all intergenic bases", outgoing[idx(State::INTERGENIC)] == 3);
        CHECK("START_CODON_1 outgoing count includes both occurrences", outgoing[idx(State::START_CODON_1)] == 2);
        CHECK("END has no outgoing transitions counted", outgoing[idx(State::END)] == 0);
    }

    static void test_transition_log_probs() {
        cout << "\n[TEST 4] Log probabilities smooth only legal topology transitions\n";

        profile.transition_alpha = 1.0;

        auto states = transition_fixture_states();
        auto ranges = transition_fixture_ranges();
        auto log_probs = Transition_Model::compute_log_probs(states, ranges);

        Log_Prob expected_seen = log((2.0 + 1.0) / (3.0 + Transitions.at(State::INTERGENIC).size()));
        Log_Prob expected_legal_unseen = log((0.0 + 1.0) / (3.0 + Transitions.at(State::INTERGENIC).size()));
        CHECK("seen transition probability is smoothed from counts",
              approx_equal(log_probs[idx(State::INTERGENIC)][idx(State::START_CODON_1)], expected_seen));
        CHECK("legal unseen transition receives alpha mass",
              approx_equal(log_probs[idx(State::INTERGENIC)][idx(State::END)], expected_legal_unseen));
        CHECK("illegal transition remains impossible",
              log_probs[idx(State::INTERGENIC)][idx(State::STOP_CODON_1)] == LOG_ZERO);
    }

    static void run_Transition_Model_tests() {
        cout << "\nRunning Transition_Model tests...\n";

        test_transition_bigram_counts();
        test_transition_no_cross_chromosome_bigrams();
        test_transition_outgoing_counts();
        test_transition_log_probs();

        cout << "\nDone.\n";
    }

}
