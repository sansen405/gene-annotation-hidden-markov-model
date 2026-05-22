#pragma once

#include "../model/Emission_Model.hpp"
#include "../genome_profiles/Genome_Profile.hpp"
#include "Test_Utils.hpp"
#include <cmath>
#include <string>
#include <vector>

namespace gene_hmm {

    using namespace std;

    static bool emission_approx_equal(Log_Prob actual, Log_Prob expected, double eps = 1e-12) {
        return fabs(actual - expected) < eps;
    }

    static size_t test_encode_5mer(const vector<Nucleotide>& nucs, size_t pos) {
        size_t context = 0;
        for (size_t i = 0; i < 5; ++i) {
            context = context * NUM_NUCLEOTIDES + idx(nucs[pos - 5 + i]);
        }
        return context;
    }

    static void test_emission_markov1_counts() {
        cout << "\n[TEST 1] Markov1 emission counts use previous nucleotide context\n";

        vector<Nucleotide> nucs = {
            Nucleotide::A, Nucleotide::C, Nucleotide::G,
            Nucleotide::T, Nucleotide::A
        };
        vector<State> states = {
            State::INTERGENIC, State::INTERGENIC, State::INTRON_1,
            State::INTERGENIC, State::INTERGENIC
        };
        vector<Chromosome_Range> ranges = {{"chr1", 0, nucs.size()}};

        auto counts = Emission_Model::count_markov1_emissions(states, nucs, ranges, {State::INTERGENIC});
        CHECK("A -> C intergenic emission counted", counts[idx(Nucleotide::A)][idx(Nucleotide::C)] == 1);
        CHECK("G -> T intergenic emission counted", counts[idx(Nucleotide::G)][idx(Nucleotide::T)] == 1);
        CHECK("T -> A intergenic emission counted", counts[idx(Nucleotide::T)][idx(Nucleotide::A)] == 1);
        CHECK("C -> G skipped because target state is intron", counts[idx(Nucleotide::C)][idx(Nucleotide::G)] == 0);
    }

    static void test_emission_markov5_counts() {
        cout << "\n[TEST 2] Markov5 emission counts use 5-base context\n";

        vector<Nucleotide> nucs = {
            Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T,
            Nucleotide::A, Nucleotide::C, Nucleotide::G
        };
        vector<State> states = {
            State::INTERGENIC, State::INTERGENIC, State::INTERGENIC,
            State::INTERGENIC, State::INTERGENIC, State::EXON_FRAME_1,
            State::INTERGENIC
        };
        vector<Chromosome_Range> ranges = {{"chr1", 0, nucs.size()}};

        auto counts = Emission_Model::count_markov5_emissions(states, nucs, ranges, {State::EXON_FRAME_1});
        size_t context = test_encode_5mer(nucs, 5);
        CHECK("ACGTA -> C exon emission counted", counts[context][idx(Nucleotide::C)] == 1);
    }

    static void test_emission_pssm_counts() {
        cout << "\n[TEST 3] PSSM emission counts collect window columns\n";

        vector<Nucleotide> nucs = {
            Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T,
            Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T
        };
        vector<State> states(nucs.size(), State::INTERGENIC);
        states[3] = State::DONOR_1;
        states[7] = State::DONOR_1; // Boundary case should be skipped.
        vector<Chromosome_Range> ranges = {{"chr1", 0, nucs.size()}};

        auto counts = Emission_Model::count_pssm_emissions(states, nucs, ranges, {State::DONOR_1}, 2, 3);
        CHECK("PSSM has left + right columns", counts.size() == 5);
        CHECK("column 0 counts C", counts[0][idx(Nucleotide::C)] == 1);
        CHECK("column 1 counts G", counts[1][idx(Nucleotide::G)] == 1);
        CHECK("column 2 counts T", counts[2][idx(Nucleotide::T)] == 1);
        CHECK("column 3 counts A", counts[3][idx(Nucleotide::A)] == 1);
        CHECK("column 4 counts C", counts[4][idx(Nucleotide::C)] == 1);
    }

    static void test_emission_log_probs() {
        cout << "\n[TEST 4] Emission log probabilities use additive smoothing\n";

        profile.emission_alpha = 1.0;

        Emission_Model::Markov1_Count counts = {};
        counts[idx(Nucleotide::A)][idx(Nucleotide::C)] = 3;

        auto log_probs = Emission_Model::compute_markov1_log_probs(counts);
        CHECK("seen Markov1 emission probability is smoothed",
              emission_approx_equal(log_probs[idx(Nucleotide::A)][idx(Nucleotide::C)], log(4.0 / 7.0)));
        CHECK("unseen Markov1 emission probability receives alpha mass",
              emission_approx_equal(log_probs[idx(Nucleotide::A)][idx(Nucleotide::A)], log(1.0 / 7.0)));

        Emission_Model::PSSM_Count pssm_counts(1, array<uint64_t, NUM_NUCLEOTIDES>{});
        pssm_counts[0][idx(Nucleotide::G)] = 2;
        auto pssm_log_probs = Emission_Model::compute_pssm_log_probs(pssm_counts);
        CHECK("PSSM column probabilities are smoothed",
              emission_approx_equal(pssm_log_probs[0][idx(Nucleotide::G)], log(3.0 / 6.0)));
    }

    static void test_emission_deterministic_codons() {
        cout << "\n[TEST 5] Deterministic codon emissions enforce expected bases\n";

        CHECK("START_CODON_1 emits A with log probability 0",
              Emission_Model::get_deterministic_log_prob(State::START_CODON_1, Nucleotide::A) == 0.0);
        CHECK("START_CODON_2 rejects A",
              Emission_Model::get_deterministic_log_prob(State::START_CODON_2, Nucleotide::A) == LOG_ZERO);
        CHECK("STOP_CODON_2 allows A with half probability",
              emission_approx_equal(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_2, Nucleotide::A), log(0.5)));
        CHECK("STOP_CODON_3 rejects C",
              Emission_Model::get_deterministic_log_prob(State::STOP_CODON_3, Nucleotide::C) == LOG_ZERO);
    }

    static void test_emission_dispatch_defaults() {
        cout << "\n[TEST 6] emission_log_prob dispatches by state family\n";

        Emission_Model model;
        vector<Nucleotide> nucs = {
            Nucleotide::A, Nucleotide::T, Nucleotide::G, Nucleotide::C,
            Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T,
            Nucleotide::A, Nucleotide::C
        };

        CHECK("INTERGENIC uses uniform initial emission at t=0",
              emission_approx_equal(model.emission_log_prob(State::INTERGENIC, 0, nucs), log(0.25)));
        CHECK("EXON_FRAME uses uniform fallback before 5-base context exists",
              emission_approx_equal(model.emission_log_prob(State::EXON_FRAME_1, 4, nucs), log(0.25)));
        CHECK("DONOR PSSM averages the full default window when available",
              emission_approx_equal(model.emission_log_prob(State::DONOR_1, 3, nucs), log(0.25)));
        CHECK("DONOR PSSM returns LOG_ZERO near boundary",
              model.emission_log_prob(State::DONOR_1, 1, nucs) == LOG_ZERO);
        CHECK("START_CODON_3 dispatches to deterministic G",
              model.emission_log_prob(State::START_CODON_3, 2, nucs) == 0.0);
    }

    static void run_Emission_Model_tests() {
        cout << "\nRunning Emission_Model tests...\n";

        test_emission_markov1_counts();
        test_emission_markov5_counts();
        test_emission_pssm_counts();
        test_emission_log_probs();
        test_emission_deterministic_codons();
        test_emission_dispatch_defaults();

        cout << "\nDone.\n";
    }

}
