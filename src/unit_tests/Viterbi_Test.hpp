#pragma once

#include "../decoding/Viterbi.hpp"
#include "Test_Utils.hpp"
#include <cmath>
#include <string>
#include <vector>

namespace gene_hmm {

    using namespace std;

    static Transition_Model::Log_Prob_Matrix make_viterbi_log_zero_matrix() {
        Transition_Model::Log_Prob_Matrix matrix = {};
        for (auto& row : matrix) {
            for (auto& value : row) {
                value = LOG_ZERO;
            }
        }
        return matrix;
    }

    static void test_viterbi_empty_sequence() {
        cout << "\n[TEST 1] Empty sequence returns empty path\n";

        Emission_Model model;
        auto transitions = make_viterbi_log_zero_matrix();
        vector<State> path = Viterbi::decode({}, transitions, model);

        CHECK("empty input decodes to empty annotation", path.empty());
    }

    static void test_viterbi_deterministic_start_codon_path() {
        cout << "\n[TEST 2] Deterministic path decodes ATG start codon\n";

        Emission_Model model;
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::T, Nucleotide::G};
        auto transitions = make_viterbi_log_zero_matrix();
        transitions[idx(State::START)][idx(State::START_CODON_1)] = 0.0;
        transitions[idx(State::START_CODON_1)][idx(State::START_CODON_2)] = 0.0;
        transitions[idx(State::START_CODON_2)][idx(State::START_CODON_3)] = 0.0;

        vector<State> path = Viterbi::decode(nucs, transitions, model);
        vector<State> expected = {
            State::START_CODON_1,
            State::START_CODON_2,
            State::START_CODON_3
        };

        CHECK("decoded path length matches sequence length", path.size() == nucs.size());
        CHECK("decoded path follows forced start-codon states", path == expected);
    }

    static void test_viterbi_intergenic_self_loop_path() {
        cout << "\n[TEST 3] Intergenic self-loop decodes repeated intergenic states\n";

        Emission_Model model;
        vector<Nucleotide> nucs = {Nucleotide::C, Nucleotide::G, Nucleotide::T, Nucleotide::A};
        auto transitions = make_viterbi_log_zero_matrix();
        transitions[idx(State::START)][idx(State::INTERGENIC)] = 0.0;
        transitions[idx(State::INTERGENIC)][idx(State::INTERGENIC)] = 0.0;

        vector<State> path = Viterbi::decode(nucs, transitions, model);

        CHECK("decoded path length matches sequence length", path.size() == nucs.size());
        CHECK("all decoded states are INTERGENIC",
              path == vector<State>(nucs.size(), State::INTERGENIC));
    }

    static void test_viterbi_chooses_higher_scoring_path() {
        cout << "\n[TEST 4] Decoder chooses higher-scoring transition path\n";

        Emission_Model model;
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::T, Nucleotide::G};
        auto transitions = make_viterbi_log_zero_matrix();

        transitions[idx(State::START)][idx(State::INTERGENIC)] = log(0.1);
        transitions[idx(State::INTERGENIC)][idx(State::INTERGENIC)] = log(0.1);

        transitions[idx(State::START)][idx(State::START_CODON_1)] = log(0.9);
        transitions[idx(State::START_CODON_1)][idx(State::START_CODON_2)] = log(0.9);
        transitions[idx(State::START_CODON_2)][idx(State::START_CODON_3)] = log(0.9);

        vector<State> path = Viterbi::decode(nucs, transitions, model);

        CHECK("higher-scoring ATG path wins over intergenic alternative",
              path == vector<State>({
                  State::START_CODON_1,
                  State::START_CODON_2,
                  State::START_CODON_3
              }));
    }

    static void run_Viterbi_tests() {
        cout << "\nRunning Viterbi tests...\n";

        test_viterbi_empty_sequence();
        test_viterbi_deterministic_start_codon_path();
        test_viterbi_intergenic_self_loop_path();
        test_viterbi_chooses_higher_scoring_path();

        cout << "\nDone.\n";
    }

}
