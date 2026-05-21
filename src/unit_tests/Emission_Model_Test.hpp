#pragma once

#include "../model/Emission_Model.hpp"
#include "test_utils.hpp"
#include <array>
#include <cmath>
#include <vector>

namespace gene_hmm {

    using namespace std;

    static bool approx(double a, double b, double tol = 1e-9) {
        return fabs(a - b) < tol;
    }

    static Chromosome_Range make_chr(size_t start, size_t end) {
        return {"chr", start, end};
    }

    static size_t encode_5mer_ref(const vector<Nucleotide>& nucs, size_t pos) {
        size_t ctx = 0;
        for (size_t i = 0; i < 5; ++i)
            ctx = ctx * NUM_NUCLEOTIDES + (static_cast<size_t>(nucs[pos - 5 + i]) - 1);
        return ctx;
    }

    static void test_markov1_count() {
        cout << "\n[TEST 1] count_markov1_emissions\n";

        {
            vector<State> states = {State::INTERGENIC};
            vector<Nucleotide> nucs = {Nucleotide::A};
            auto counts = Emission_Model::count_markov1_emissions(states, nucs, {}, {State::INTERGENIC});
            bool all_zero = true;
            for (auto& row : counts) for (auto v : row) if (v) all_zero = false;
            CHECK("empty chromosome range → all zeros", all_zero);
        }

        {
            vector<State> states = {State::INTERGENIC};
            vector<Nucleotide> nucs = {Nucleotide::A};
            auto counts = Emission_Model::count_markov1_emissions(states, nucs, {make_chr(0, 1)}, {State::INTERGENIC});
            bool all_zero = true;
            for (auto& row : counts) for (auto v : row) if (v) all_zero = false;
            CHECK("single nucleotide → no counts", all_zero);
        }

        {
            vector<State> states = {State::INTERGENIC, State::INTERGENIC};
            vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C};
            auto counts = Emission_Model::count_markov1_emissions(states, nucs, {make_chr(0, 2)}, {State::INTERGENIC});
            CHECK("A→C: counts[A][C] == 1", counts[0][1] == 1);
            uint64_t total = 0;
            for (auto& row : counts) for (auto v : row) total += v;
            CHECK("A→C: total count == 1", total == 1);
        }

        {
            vector<State> states = {State::INTERGENIC, State::EXON_FRAME_1, State::INTERGENIC};
            vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G};
            auto counts = Emission_Model::count_markov1_emissions(states, nucs, {make_chr(0, 3)}, {State::INTERGENIC});
            // pos 2 (INTERGENIC, target): context = C (idx 1), emit = G (idx 2)
            CHECK("non-target skipped: counts[C][G] == 1", counts[1][2] == 1);
            uint64_t total = 0;
            for (auto& row : counts) for (auto v : row) total += v;
            CHECK("non-target skipped: total == 1", total == 1);
        }

        {
            vector<State> states = {State::INTERGENIC, State::INTERGENIC,
                                     State::INTERGENIC, State::INTERGENIC};
            vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T};
            vector<Chromosome_Range> chrs = {{"chr1", 0, 2}, {"chr2", 2, 4}};
            auto counts = Emission_Model::count_markov1_emissions(states, nucs, chrs, {State::INTERGENIC});
            CHECK("boundary: chr1 A→C counted", counts[0][1] == 1);
            CHECK("boundary: chr2 G→T counted", counts[2][3] == 1);
            uint64_t total = 0;
            for (auto& row : counts) for (auto v : row) total += v;
            CHECK("boundary: cross-chromosome pair not counted (total == 2)", total == 2);
        }

        {
            vector<State> states = {State::INTRON_1, State::INTRON_2};
            vector<Nucleotide> nucs = {Nucleotide::G, Nucleotide::T};
            auto counts = Emission_Model::count_markov1_emissions(
                states, nucs, {make_chr(0, 2)}, {State::INTRON_1, State::INTRON_2});
            CHECK("multiple targets pooled: counts[G][T] == 1", counts[2][3] == 1);
        }

        {
            vector<State> states(6, State::INTERGENIC);
            vector<Nucleotide> nucs(6, Nucleotide::A);
            auto counts = Emission_Model::count_markov1_emissions(states, nucs, {make_chr(0, 6)}, {State::INTERGENIC});
            CHECK("6 A's → 5 A→A transitions", counts[0][0] == 5);
        }
    }


    static void test_markov5_count() {
        cout << "\n[TEST 2] count_markov5_emissions\n";

        {
            auto counts = Emission_Model::count_markov5_emissions({State::EXON_FRAME_1}, {Nucleotide::A}, {}, {State::EXON_FRAME_1});
            bool all_zero = true;
            for (auto& row : counts) for (auto v : row) if (v) all_zero = false;
            CHECK("empty range → all zeros", all_zero);
        }


        {
            vector<State> states(5, State::EXON_FRAME_1);
            vector<Nucleotide> nucs(5, Nucleotide::A);
            auto counts = Emission_Model::count_markov5_emissions(states, nucs, {make_chr(0, 5)}, {State::EXON_FRAME_1});
            bool all_zero = true;
            for (auto& row : counts) for (auto v : row) if (v) all_zero = false;
            CHECK("5 positions → no counts", all_zero);
        }

  
        {
            vector<State> states(6, State::EXON_FRAME_1);
            vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::A, Nucleotide::A,
                                        Nucleotide::A, Nucleotide::A, Nucleotide::C};
            auto counts = Emission_Model::count_markov5_emissions(states, nucs, {make_chr(0, 6)}, {State::EXON_FRAME_1});
            uint64_t total = 0;
            for (auto& row : counts) for (auto v : row) total += v;
            CHECK("6 positions → exactly one count", total == 1);
            CHECK("AAAAAC: context 0 (AAAAA), emit C (idx 1)", counts[0][1] == 1);
        }


        {
            vector<State> states(7, State::EXON_FRAME_1);
            vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T,
                                        Nucleotide::A, Nucleotide::T, Nucleotide::A};
            auto counts = Emission_Model::count_markov5_emissions(states, nucs, {make_chr(0, 7)}, {State::EXON_FRAME_1});
            size_t ctx = encode_5mer_ref(nucs, 5);
            CHECK("ACGTA context: emit T at correct index", counts[ctx][3] == 1);
        }

   
        {
            vector<State> states = {State::EXON_FRAME_1, State::EXON_FRAME_1, State::EXON_FRAME_1,
                                     State::EXON_FRAME_1, State::EXON_FRAME_1, State::INTERGENIC};
            vector<Nucleotide> nucs(6, Nucleotide::A);
            nucs[5] = Nucleotide::C;
            auto counts = Emission_Model::count_markov5_emissions(states, nucs, {make_chr(0, 6)}, {State::EXON_FRAME_1});
            bool all_zero = true;
            for (auto& row : counts) for (auto v : row) if (v) all_zero = false;
            CHECK("emit position is non-target → no count", all_zero);
        }


        {
            vector<State> states(12, State::EXON_FRAME_1);
            vector<Nucleotide> nucs(12, Nucleotide::A);
            nucs[11] = Nucleotide::C;
            vector<Chromosome_Range> chrs = {{"chr1", 0, 6}, {"chr2", 6, 12}};
            auto counts = Emission_Model::count_markov5_emissions(states, nucs, chrs, {State::EXON_FRAME_1});
            CHECK("boundary: chr1 AAAAA→A counted", counts[0][0] == 1);
            CHECK("boundary: chr2 AAAAA→C counted", counts[0][1] == 1);
            uint64_t total = 0;
            for (auto& row : counts) for (auto v : row) total += v;
            CHECK("boundary: total == 2, no cross-chromosome context", total == 2);
        }

        {
            vector<State> states = {State::EXON_FRAME_1, State::EXON_FRAME_1, State::EXON_FRAME_1,
                                     State::EXON_FRAME_1, State::EXON_FRAME_1,
                                     State::EXON_FRAME_1, State::EXON_FRAME_2, State::EXON_FRAME_3};
            vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::A, Nucleotide::A,
                                        Nucleotide::A, Nucleotide::A,
                                        Nucleotide::C, Nucleotide::G, Nucleotide::T};
            auto counts = Emission_Model::count_markov5_emissions(
                states, nucs, {make_chr(0, 8)},
                {State::EXON_FRAME_1, State::EXON_FRAME_2, State::EXON_FRAME_3});
            uint64_t total = 0;
            for (auto& row : counts) for (auto v : row) total += v;
            CHECK("frame_tied: all 3 frame states contribute (total == 3)", total == 3);
        }
    }



    static void test_pssm_count() {
        cout << "\n[TEST 3] count_pssm_emissions\n";

        
        {
            vector<State> states(6, State::INTERGENIC);
            states[2] = State::DONOR_1;
            vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G,
                                        Nucleotide::T, Nucleotide::A, Nucleotide::C};
            auto counts = Emission_Model::count_pssm_emissions(
                states, nucs, {make_chr(0, 6)}, {State::DONOR_1}, 2, 2);
            CHECK("window_size == 4", counts.size() == 4);
            CHECK("col 0 = A", counts[0][0] == 1);
            CHECK("col 1 = C", counts[1][1] == 1);
            CHECK("col 2 = G (target position)", counts[2][2] == 1);
            CHECK("col 3 = T", counts[3][3] == 1);
        }


        {
            vector<State> states = {State::INTERGENIC, State::DONOR_1,
                                     State::INTERGENIC, State::INTERGENIC, State::INTERGENIC, State::INTERGENIC};
            vector<Nucleotide> nucs(6, Nucleotide::A);
            auto counts = Emission_Model::count_pssm_emissions(
                states, nucs, {make_chr(0, 6)}, {State::DONOR_1}, 3, 2);
            bool all_zero = true;
            for (auto& col : counts) for (auto v : col) if (v) all_zero = false;
            CHECK("target too close to chr start → skipped", all_zero);
        }


        {
            vector<State> states(6, State::INTERGENIC);
            states[4] = State::DONOR_1;
            vector<Nucleotide> nucs(6, Nucleotide::A);
            auto counts = Emission_Model::count_pssm_emissions(
                states, nucs, {make_chr(0, 6)}, {State::DONOR_1}, 2, 3);
            bool all_zero = true;
            for (auto& col : counts) for (auto v : col) if (v) all_zero = false;
            CHECK("target too close to chr end → skipped", all_zero);
        }


        {
            vector<State> states(6, State::INTERGENIC);
            states[2] = State::DONOR_1;
            vector<Nucleotide> nucs(6, Nucleotide::T);
            auto counts = Emission_Model::count_pssm_emissions(
                states, nucs, {make_chr(0, 6)}, {State::DONOR_1}, 2, 2);
            uint64_t total = 0;
            for (auto& col : counts) for (auto v : col) total += v;
            CHECK("exactly at left boundary → 4 columns counted", total == 4);
        }


        {
            vector<State> states(8, State::INTERGENIC);
            states[5] = State::DONOR_1;
            vector<Nucleotide> nucs(8, Nucleotide::A);
            vector<Chromosome_Range> chrs = {{"chr1", 0, 4}, {"chr2", 4, 8}};
            auto counts = Emission_Model::count_pssm_emissions(
                states, nucs, chrs, {State::DONOR_1}, 3, 2);
            bool all_zero = true;
            for (auto& col : counts) for (auto v : col) if (v) all_zero = false;
            CHECK("window crosses chromosome boundary → skipped", all_zero);
        }

    
        {
            vector<State> states = {State::INTERGENIC, State::INTERGENIC, State::DONOR_1,
                                     State::INTERGENIC, State::DONOR_1, State::INTERGENIC};
            vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G,
                                        Nucleotide::T, Nucleotide::A, Nucleotide::C};
            auto counts = Emission_Model::count_pssm_emissions(
                states, nucs, {make_chr(0, 6)}, {State::DONOR_1}, 1, 1);
            // pos 2: nucs[1]=C(1), nucs[2]=G(2)
            // pos 4: nucs[3]=T(3), nucs[4]=A(0)
            CHECK("multi: col0 C", counts[0][1] == 1);
            CHECK("multi: col1 G", counts[1][2] == 1);
            CHECK("multi: col0 T", counts[0][3] == 1);
            CHECK("multi: col1 A", counts[1][0] == 1);
        }
    }


    static void test_markov1_log_probs() {
        cout << "\n[TEST 4] compute_markov1_log_probs\n";


        {
            Emission_Model::Markov1_Count counts = {};
            auto lp = Emission_Model::compute_markov1_log_probs(counts, 1.0);
            bool all_uniform = true;
            double expected = log(1.0 / NUM_NUCLEOTIDES);
            for (auto& row : lp) for (auto v : row) if (!approx(v, expected)) all_uniform = false;
            CHECK("zero counts + alpha=1 → uniform log(0.25)", all_uniform);
        }


        {
            Emission_Model::Markov1_Count counts = {};
            counts[0][0] = 10; counts[0][1] = 5; counts[0][2] = 3; counts[0][3] = 2;
            auto lp = Emission_Model::compute_markov1_log_probs(counts, 0.5);
            bool rows_ok = true;
            for (auto& row : lp) {
                double sum = 0;
                for (auto v : row) sum += exp(v);
                if (!approx(sum, 1.0, 1e-9)) rows_ok = false;
            }
            CHECK("all rows sum to 1.0 in probability space", rows_ok);
        }


        {
            Emission_Model::Markov1_Count counts = {};
            counts[0][0] = 3;
            auto lp = Emission_Model::compute_markov1_log_probs(counts, 1.0);
            double expected = log((3.0 + 1.0) / (3.0 + 1.0 * NUM_NUCLEOTIDES));
            CHECK("specific count: log P(A|A) correct", approx(lp[0][0], expected));
        }


        {
            Emission_Model::Markov1_Count counts = {};
            for (auto& row : counts) for (auto& v : row) v = 10;
            auto lp = Emission_Model::compute_markov1_log_probs(counts, 0.0);
            bool all_eq = true;
            for (auto& row : lp) for (auto v : row) if (!approx(v, log(0.25))) all_eq = false;
            CHECK("uniform counts, alpha=0 → log(0.25) everywhere", all_eq);
        }
    }


    static void test_markov5_log_probs() {
        cout << "\n[TEST 5] compute_markov5_log_probs\n";

        {
            Emission_Model::Markov5_Count counts = {};
            auto lp = Emission_Model::compute_markov5_log_probs(counts, 1.0);
            double expected = log(1.0 / NUM_NUCLEOTIDES);
            bool all_uniform = true;
            for (auto& row : lp) for (auto v : row) if (!approx(v, expected)) all_uniform = false;
            CHECK("zero counts + alpha=1 → uniform across all 1024 contexts", all_uniform);
        }

        {
            Emission_Model::Markov5_Count counts = {};
            counts[108][0] = 5; counts[108][1] = 3; counts[108][2] = 2; counts[108][3] = 1;
            auto lp = Emission_Model::compute_markov5_log_probs(counts, 0.5);
            bool all_ok = true;
            for (size_t ctx = 0; ctx < Emission_Model::MARKOV5_CONTEXTS; ++ctx) {
                double sum = 0;
                for (auto v : lp[ctx]) sum += exp(v);
                if (!approx(sum, 1.0, 1e-9)) all_ok = false;
            }
            CHECK("all 1024 rows sum to 1.0 in probability space", all_ok);
        }

        {
            Emission_Model::Markov5_Count counts = {};
            counts[0][0] = 100;
            auto lp = Emission_Model::compute_markov5_log_probs(counts, 0.5);
            double expected = log((100.0 + 0.5) / (100.0 + 0.5 * NUM_NUCLEOTIDES));
            CHECK("specific row: log P(A|AAAAA) correct", approx(lp[0][0], expected));
        }

        {
            Emission_Model::Markov5_Count counts = {};
            counts[0][0] = 100;
            auto lp = Emission_Model::compute_markov5_log_probs(counts, 1.0);
            bool unseen_ok = true;
            for (size_t ctx = 1; ctx < Emission_Model::MARKOV5_CONTEXTS; ++ctx)
                for (auto v : lp[ctx])
                    if (!approx(v, log(0.25))) unseen_ok = false;
            CHECK("unseen contexts → uniform log probs via smoothing", unseen_ok);
        }
    }


    static void test_pssm_log_probs() {
        cout << "\n[TEST 6] compute_pssm_log_probs\n";

        {
            Emission_Model::PSSM_Count counts(9, array<uint64_t, NUM_NUCLEOTIDES>{});
            auto lp = Emission_Model::compute_pssm_log_probs(counts, 1.0);
            bool all_uniform = true;
            for (auto& col : lp) for (auto v : col) if (!approx(v, log(0.25))) all_uniform = false;
            CHECK("zero counts + alpha=1 → uniform log(0.25) per column", all_uniform);
        }

        {
            Emission_Model::PSSM_Count counts(9, array<uint64_t, NUM_NUCLEOTIDES>{});
            counts[0][0] = 10; counts[0][1] = 5;
            counts[3][2] = 7;  counts[3][3] = 3;
            auto lp = Emission_Model::compute_pssm_log_probs(counts, 0.5);
            bool all_ok = true;
            for (auto& col : lp) {
                double sum = 0;
                for (auto v : col) sum += exp(v);
                if (!approx(sum, 1.0, 1e-9)) all_ok = false;
            }
            CHECK("all columns sum to 1.0 in probability space", all_ok);
        }

        {
            Emission_Model::PSSM_Count counts(4, array<uint64_t, NUM_NUCLEOTIDES>{});
            counts[2][0] = 8;
            auto lp = Emission_Model::compute_pssm_log_probs(counts, 1.0);
            double expected = log((8.0 + 1.0) / (8.0 + 1.0 * NUM_NUCLEOTIDES));
            CHECK("specific column: log P(A|col2) correct", approx(lp[2][0], expected));
        }

        {
            Emission_Model::PSSM_Count counts(2, array<uint64_t, NUM_NUCLEOTIDES>{});
            counts[0][0] = 100;
            auto lp = Emission_Model::compute_pssm_log_probs(counts, 1.0);
            bool col1_ok = true;
            for (auto v : lp[1]) if (!approx(v, log(0.25))) col1_ok = false;
            CHECK("column independence: col1 uniform despite col0 counts", col1_ok);
        }

        {
            Emission_Model::PSSM_Count counts(18, array<uint64_t, NUM_NUCLEOTIDES>{});
            auto lp = Emission_Model::compute_pssm_log_probs(counts, 0.5);
            CHECK("output window size preserved (18 = acceptor window)", lp.size() == 18);
        }
    }


    static void test_deterministic_log_prob() {
        cout << "\n[TEST 7] get_deterministic_log_prob\n";

        CHECK("START_CODON_1 + A → 0.0",
              approx(Emission_Model::get_deterministic_log_prob(State::START_CODON_1, Nucleotide::A), 0.0));
        CHECK("START_CODON_1 + C → LOG_ZERO",
              Emission_Model::get_deterministic_log_prob(State::START_CODON_1, Nucleotide::C) == LOG_ZERO);
        CHECK("START_CODON_1 + G → LOG_ZERO",
              Emission_Model::get_deterministic_log_prob(State::START_CODON_1, Nucleotide::G) == LOG_ZERO);
        CHECK("START_CODON_1 + T → LOG_ZERO",
              Emission_Model::get_deterministic_log_prob(State::START_CODON_1, Nucleotide::T) == LOG_ZERO);

        CHECK("START_CODON_2 + T → 0.0",
              approx(Emission_Model::get_deterministic_log_prob(State::START_CODON_2, Nucleotide::T), 0.0));
        CHECK("START_CODON_2 + A → LOG_ZERO",
              Emission_Model::get_deterministic_log_prob(State::START_CODON_2, Nucleotide::A) == LOG_ZERO);

        CHECK("START_CODON_3 + G → 0.0",
              approx(Emission_Model::get_deterministic_log_prob(State::START_CODON_3, Nucleotide::G), 0.0));
        CHECK("START_CODON_3 + A → LOG_ZERO",
              Emission_Model::get_deterministic_log_prob(State::START_CODON_3, Nucleotide::A) == LOG_ZERO);

        CHECK("STOP_CODON_1 + T → 0.0",
              approx(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_1, Nucleotide::T), 0.0));
        CHECK("STOP_CODON_1 + A → LOG_ZERO",
              Emission_Model::get_deterministic_log_prob(State::STOP_CODON_1, Nucleotide::A) == LOG_ZERO);

        CHECK("STOP_CODON_2 + A → log(0.5)",
              approx(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_2, Nucleotide::A), log(0.5)));
        CHECK("STOP_CODON_2 + G → log(0.5)",
              approx(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_2, Nucleotide::G), log(0.5)));
        CHECK("STOP_CODON_2 + C → LOG_ZERO",
              Emission_Model::get_deterministic_log_prob(State::STOP_CODON_2, Nucleotide::C) == LOG_ZERO);
        CHECK("STOP_CODON_2 + T → LOG_ZERO",
              Emission_Model::get_deterministic_log_prob(State::STOP_CODON_2, Nucleotide::T) == LOG_ZERO);

        CHECK("STOP_CODON_3 + A → log(0.5)",
              approx(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_3, Nucleotide::A), log(0.5)));
        CHECK("STOP_CODON_3 + G → log(0.5)",
              approx(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_3, Nucleotide::G), log(0.5)));
        CHECK("STOP_CODON_3 + C → LOG_ZERO",
              Emission_Model::get_deterministic_log_prob(State::STOP_CODON_3, Nucleotide::C) == LOG_ZERO);

        CHECK("INTERGENIC → 0.0 (no-op)",
              approx(Emission_Model::get_deterministic_log_prob(State::INTERGENIC, Nucleotide::A), 0.0));
        CHECK("EXON_FRAME_1 → 0.0 (no-op)",
              approx(Emission_Model::get_deterministic_log_prob(State::EXON_FRAME_1, Nucleotide::C), 0.0));
        CHECK("DONOR_1 → 0.0 (PSSM handled as transition bonus)",
              approx(Emission_Model::get_deterministic_log_prob(State::DONOR_1, Nucleotide::G), 0.0));
    }


    static void run_Emission_Model_tests() {
        cout << "=== Emission_Model Tests ===\n";
        test_markov1_count();
        test_markov5_count();
        test_pssm_count();
        test_markov1_log_probs();
        test_markov5_log_probs();
        test_pssm_log_probs();
        test_deterministic_log_prob();
        cout << "\nDone.\n";
    }

}
