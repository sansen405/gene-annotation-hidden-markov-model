#pragma once

#include "../decoding/Forward_Backward.hpp"
#include "Test_Utils.hpp"
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace gene_hmm {

    using namespace std;

    static Transition_Model::Log_Prob_Matrix make_forward_backward_log_zero_matrix() {
        Transition_Model::Log_Prob_Matrix matrix = {};
        for(auto& row : matrix){
            for(auto& value : row){
                value = LOG_ZERO;
            }
        }
        return matrix;
    }

    static void test_forward_backward_empty_sequence() {
        cout << "\n[TEST 1] Empty sequence returns empty posterior matrix\n";

        Emission_Model model;
        auto transitions = make_forward_backward_log_zero_matrix();
        auto posterior = Forward_Backward::posterior_log_probs({}, transitions, model);

        CHECK("empty input has no posterior rows", posterior.empty());
    }

    static void test_forward_backward_forced_path_confidence() {
        cout << "\n[TEST 2] Forced intergenic path has confidence 1 at each base\n";

        Emission_Model model;
        //donor emissions require loaded splice scores; an empty file gives neutral 0.0
        const string splice_scores_path = "/tmp/gene_hmm_fb_test_splice_scores.tsv";
        { ofstream empty_scores(splice_scores_path); }
        model.load_splice_cnn_scores(splice_scores_path, 3);

        vector<Nucleotide> nucs = {Nucleotide::C, Nucleotide::G, Nucleotide::T};
        vector<State> path = {State::INTERGENIC, State::INTERGENIC, State::INTERGENIC};
        auto transitions = make_forward_backward_log_zero_matrix();

        transitions[idx(State::START)][idx(State::INTERGENIC)] = 0.0;
        transitions[idx(State::INTERGENIC)][idx(State::INTERGENIC)] = 0.0;
        transitions[idx(State::INTERGENIC)][idx(State::END)] = 0.0;

        vector<double> confidence = Forward_Backward::confidence(nucs, path, transitions, model);

        CHECK("confidence length matches sequence length", confidence.size() == nucs.size());
        CHECK("all forced states have posterior confidence near 1",
              fabs(confidence[0] - 1.0) < 1e-9 &&
              fabs(confidence[1] - 1.0) < 1e-9 &&
              fabs(confidence[2] - 1.0) < 1e-9);
    }

    static void test_forward_backward_ambiguous_one_base_confidence() {
        cout << "\n[TEST 3] One-base posterior confidence follows transition mass\n";

        Emission_Model model;
        vector<Nucleotide> nucs = {Nucleotide::C};
        vector<State> path = {State::INTRON_1};
        auto transitions = make_forward_backward_log_zero_matrix();

        transitions[idx(State::START)][idx(State::INTERGENIC)] = log(0.25);
        transitions[idx(State::START)][idx(State::INTRON_1)] = log(0.75);
        transitions[idx(State::INTERGENIC)][idx(State::END)] = 0.0;
        transitions[idx(State::INTRON_1)][idx(State::END)] = 0.0;

        vector<double> confidence = Forward_Backward::confidence(nucs, path, transitions, model);

        CHECK("predicted intron state has posterior confidence near 0.75",
              confidence.size() == 1 && fabs(confidence[0] - 0.75) < 1e-9);
    }

    static void run_Forward_Backward_tests() {
        cout << "\nRunning Forward-Backward tests...\n";

        test_forward_backward_empty_sequence();
        test_forward_backward_forced_path_confidence();
        test_forward_backward_ambiguous_one_base_confidence();

        cout << "\nDone.\n";
    }
}

