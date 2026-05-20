#include "../src/model/Transition_Table.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

using namespace gene_hmm;
using namespace std;

// ─── helpers ────────────────────────────────────────────────────────────────

static int g_pass = 0, g_fail = 0;

static void check(bool condition, const string& name) {
    if (condition) { cout << "[PASS] " << name << "\n"; ++g_pass; }
    else           { cout << "[FAIL] " << name << "\n"; ++g_fail; }
}

static bool approx(double a, double b, double tol = 1e-9) {
    return fabs(a - b) < tol;
}

static size_t S(State s) { return static_cast<size_t>(s); }

// Verify that every row of log_prob sums to 1.0 in probability space.
// Only checks rows where outgoing_count[i] > 0 or alpha > 0 (defined rows).
static bool rows_sum_to_one(const Transition_Table::Log_Prob_Matrix& lp,
                             const Transition_Table::Row_Sum_Vector& out,
                             double alpha, double tol = 1e-9) {
    for (size_t i = 0; i < NUM_STATES; ++i) {
        if (out[i] == 0 && alpha == 0.0) continue; // 0/0 row — skip
        double row_sum = 0.0;
        for (size_t j = 0; j < NUM_STATES; ++j)
            row_sum += exp(lp[i][j]);
        if (!approx(row_sum, 1.0, tol)) return false;
    }
    return true;
}

// Check the key invariant: outgoing[i] == sum_j bigrams[i][j] for every i.
static bool outgoing_matches_bigram_sums(const Transition_Table::Count_Matrix& bm,
                                          const Transition_Table::Row_Sum_Vector& out) {
    for (size_t i = 0; i < NUM_STATES; ++i) {
        uint64_t row_sum = 0;
        for (size_t j = 0; j < NUM_STATES; ++j) row_sum += bm[i][j];
        if (row_sum != out[i]) return false;
    }
    return true;
}

// ─── count_bigrams ───────────────────────────────────────────────────────────

static void test_bigrams_single_position() {
    // One chromosome, one position.
    // Expected: START→INTERGENIC = 1, INTERGENIC→END = 1, everything else 0.
    vector<State> seq = { State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 1} };
    auto bm = Transition_Table::count_bigrams(seq, chrs);

    check(bm[S(State::START)][S(State::INTERGENIC)] == 1,
          "bigrams | single pos: START→INTERGENIC == 1");
    check(bm[S(State::INTERGENIC)][S(State::END)] == 1,
          "bigrams | single pos: INTERGENIC→END == 1");
    check(bm[S(State::START)][S(State::END)] == 0,
          "bigrams | single pos: START→END == 0");
    check(bm[S(State::INTERGENIC)][S(State::INTERGENIC)] == 0,
          "bigrams | single pos: INTERGENIC→INTERGENIC == 0");
}

static void test_bigrams_internal_transitions() {
    // One chromosome, three positions: INTERGENIC → INTERGENIC → START_CODON_1.
    // Expected:
    //   START → INTERGENIC     = 1
    //   INTERGENIC → INTERGENIC = 1
    //   INTERGENIC → START_CODON_1 = 1
    //   START_CODON_1 → END    = 1
    vector<State> seq = { State::INTERGENIC, State::INTERGENIC, State::START_CODON_1 };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 3} };
    auto bm = Transition_Table::count_bigrams(seq, chrs);

    check(bm[S(State::START)][S(State::INTERGENIC)] == 1,
          "bigrams | internal: START→INTERGENIC == 1");
    check(bm[S(State::INTERGENIC)][S(State::INTERGENIC)] == 1,
          "bigrams | internal: INTERGENIC→INTERGENIC == 1");
    check(bm[S(State::INTERGENIC)][S(State::START_CODON_1)] == 1,
          "bigrams | internal: INTERGENIC→START_CODON_1 == 1");
    check(bm[S(State::START_CODON_1)][S(State::END)] == 1,
          "bigrams | internal: START_CODON_1→END == 1");
    // No direct INTERGENIC→END because the chromosome doesn't end on INTERGENIC.
    check(bm[S(State::INTERGENIC)][S(State::END)] == 0,
          "bigrams | internal: INTERGENIC→END == 0");
}

static void test_bigrams_no_cross_chromosome_transition() {
    // Two chromosomes: [INTERGENIC, EXON_FRAME_1] and [INTERGENIC].
    // The transition EXON_FRAME_1→INTERGENIC must NOT appear; it crosses the boundary.
    vector<State> seq = { State::INTERGENIC, State::EXON_FRAME_1, State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 2}, {"chr2", 2, 3} };
    auto bm = Transition_Table::count_bigrams(seq, chrs);

    check(bm[S(State::EXON_FRAME_1)][S(State::INTERGENIC)] == 0,
          "bigrams | cross-boundary: EXON_FRAME_1→INTERGENIC forbidden");
    // chr1 end: EXON_FRAME_1→END
    check(bm[S(State::EXON_FRAME_1)][S(State::END)] == 1,
          "bigrams | cross-boundary: EXON_FRAME_1→END == 1 (chr1 tail)");
    // chr2 start: START→INTERGENIC (plus the one from chr1)
    check(bm[S(State::START)][S(State::INTERGENIC)] == 2,
          "bigrams | cross-boundary: START→INTERGENIC == 2 (one per chromosome)");
    // chr2 only position: INTERGENIC→END
    check(bm[S(State::INTERGENIC)][S(State::END)] == 1,
          "bigrams | cross-boundary: INTERGENIC→END == 1 (chr2 tail)");
    // chr1 internal: INTERGENIC→EXON_FRAME_1
    check(bm[S(State::INTERGENIC)][S(State::EXON_FRAME_1)] == 1,
          "bigrams | cross-boundary: INTERGENIC→EXON_FRAME_1 == 1 (chr1 internal)");
}

static void test_bigrams_two_chromosomes_different_start_states() {
    // chr1 = [INTERGENIC], chr2 = [EXON_FRAME_1].
    // START accumulates a count for each distinct chromosome head state.
    vector<State> seq = { State::INTERGENIC, State::EXON_FRAME_1 };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 1}, {"chr2", 1, 2} };
    auto bm = Transition_Table::count_bigrams(seq, chrs);

    check(bm[S(State::START)][S(State::INTERGENIC)] == 1,
          "bigrams | two chroms: START→INTERGENIC == 1");
    check(bm[S(State::START)][S(State::EXON_FRAME_1)] == 1,
          "bigrams | two chroms: START→EXON_FRAME_1 == 1");
    check(bm[S(State::INTERGENIC)][S(State::END)] == 1,
          "bigrams | two chroms: INTERGENIC→END == 1");
    check(bm[S(State::EXON_FRAME_1)][S(State::END)] == 1,
          "bigrams | two chroms: EXON_FRAME_1→END == 1");
    // No cross-boundary transition
    check(bm[S(State::INTERGENIC)][S(State::EXON_FRAME_1)] == 0,
          "bigrams | two chroms: INTERGENIC→EXON_FRAME_1 == 0 (no cross-boundary)");
}

static void test_bigrams_repeated_transitions_accumulate() {
    // Four INTERGENIC positions in a row — bigram[IG][IG] should be 3.
    vector<State> seq = { State::INTERGENIC, State::INTERGENIC,
                          State::INTERGENIC, State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 4} };
    auto bm = Transition_Table::count_bigrams(seq, chrs);

    check(bm[S(State::INTERGENIC)][S(State::INTERGENIC)] == 3,
          "bigrams | repeated: INTERGENIC→INTERGENIC == 3");
    check(bm[S(State::INTERGENIC)][S(State::END)] == 1,
          "bigrams | repeated: INTERGENIC→END == 1");
    check(bm[S(State::START)][S(State::INTERGENIC)] == 1,
          "bigrams | repeated: START→INTERGENIC == 1");
}

// ─── count_outgoing ──────────────────────────────────────────────────────────

static void test_outgoing_single_position() {
    vector<State> seq = { State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 1} };
    auto out = Transition_Table::count_outgoing(seq, chrs);

    check(out[S(State::START)] == 1,      "outgoing | single pos: START == 1");
    check(out[S(State::INTERGENIC)] == 1, "outgoing | single pos: INTERGENIC == 1");
    check(out[S(State::END)] == 0,        "outgoing | single pos: END == 0");
}

static void test_outgoing_start_count_equals_num_chromosomes() {
    // START gets exactly one outgoing per chromosome regardless of content.
    vector<State> seq = { State::INTERGENIC, State::INTERGENIC, State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 1}, {"chr2", 1, 2}, {"chr3", 2, 3} };
    auto out = Transition_Table::count_outgoing(seq, chrs);

    check(out[S(State::START)] == 3,
          "outgoing | START == number of chromosomes (3)");
}

static void test_outgoing_consistency_with_bigrams() {
    // Core invariant: outgoing[i] == sum_j bigrams[i][j] for all i.
    // Tested on a multi-chromosome, mixed-state sequence.
    vector<State> seq = {
        State::INTERGENIC, State::INTERGENIC, State::START_CODON_1,
        State::EXON_FRAME_1, State::INTERGENIC
    };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 3}, {"chr2", 3, 5} };
    auto bm  = Transition_Table::count_bigrams(seq, chrs);
    auto out = Transition_Table::count_outgoing(seq, chrs);

    check(outgoing_matches_bigram_sums(bm, out),
          "outgoing | consistency: outgoing[i] == sum_j bigrams[i][j] for all i");
}

static void test_outgoing_consistency_single_chromosome() {
    vector<State> seq = { State::INTERGENIC, State::INTERGENIC,
                          State::INTERGENIC, State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 4} };
    auto bm  = Transition_Table::count_bigrams(seq, chrs);
    auto out = Transition_Table::count_outgoing(seq, chrs);

    check(outgoing_matches_bigram_sums(bm, out),
          "outgoing | consistency: single chromosome repeated state");
}

// ─── compute_log_probs ──────────────────────────────────────────────────────

static void test_log_probs_values_are_non_positive() {
    // Every log probability must be ≤ 0 (i.e., probability ≤ 1).
    vector<State> seq = { State::INTERGENIC, State::INTERGENIC, State::START_CODON_1 };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 3} };
    auto lp = Transition_Table::compute_log_probs(seq, chrs, 1.0);

    bool all_nonpositive = true;
    for (size_t i = 0; i < NUM_STATES; ++i)
        for (size_t j = 0; j < NUM_STATES; ++j)
            if (lp[i][j] > 1e-12) { all_nonpositive = false; break; }

    check(all_nonpositive, "log_probs | all values ≤ 0");
}

static void test_log_probs_rows_normalize_with_smoothing() {
    // With alpha > 0, every row must sum to 1.0 in probability space.
    vector<State> seq = {
        State::INTERGENIC, State::INTERGENIC, State::START_CODON_1,
        State::EXON_FRAME_1, State::INTERGENIC
    };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 3}, {"chr2", 3, 5} };
    double alpha = 1.0;
    auto lp  = Transition_Table::compute_log_probs(seq, chrs, alpha);
    auto out = Transition_Table::count_outgoing(seq, chrs);

    check(rows_sum_to_one(lp, out, alpha),
          "log_probs | rows normalize to 1.0 (alpha=1.0)");
}

static void test_log_probs_rows_normalize_small_alpha() {
    vector<State> seq = { State::INTERGENIC, State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 2} };
    double alpha = 0.01;
    auto lp  = Transition_Table::compute_log_probs(seq, chrs, alpha);
    auto out = Transition_Table::count_outgoing(seq, chrs);

    check(rows_sum_to_one(lp, out, alpha, 1e-9),
          "log_probs | rows normalize to 1.0 (alpha=0.01)");
}

static void test_log_probs_mle_values_alpha_zero() {
    // With alpha=0, check exact MLE values for states that have transitions.
    // seq: INTERGENIC(x3), INTERGENIC(x2) = two chromosomes
    //   chr1: INTERGENIC → INTERGENIC → INTERGENIC
    //   chr2: INTERGENIC → INTERGENIC
    // Bigrams for INTERGENIC:
    //   IG→IG = 1 (chr1 internal) + 1 (chr1 internal) + 1 (chr2 internal) = 3
    //   IG→END = 1 (chr1 tail) + 1 (chr2 tail) = 2
    //   outgoing[IG] = 5
    // P(IG|IG) = 3/5, P(END|IG) = 2/5
    vector<State> seq = { State::INTERGENIC, State::INTERGENIC, State::INTERGENIC,
                          State::INTERGENIC, State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 3}, {"chr2", 3, 5} };
    auto lp = Transition_Table::compute_log_probs(seq, chrs, 0.0);

    check(approx(lp[S(State::INTERGENIC)][S(State::INTERGENIC)], log(3.0/5.0)),
          "log_probs | MLE alpha=0: log P(IG|IG) = log(3/5)");
    check(approx(lp[S(State::INTERGENIC)][S(State::END)], log(2.0/5.0)),
          "log_probs | MLE alpha=0: log P(END|IG) = log(2/5)");
}

static void test_log_probs_laplace_smoothed_values() {
    // Verify the exact Laplace-smoothed formula for a known sequence.
    // seq: [INTERGENIC, INTERGENIC], chr = [{0,2}]
    // Bigrams:  START→IG=1, IG→IG=1, IG→END=1
    // Outgoing: START=1, IG=2
    // With alpha=1:
    //   log P(j|START) = log((bigram[START][j]+1) / (1 + 1*NUM_STATES))
    //   log P(IG|START) = log(2 / (1 + 21)) = log(2/22)
    //   log P(IG|IG)    = log((1+1)/(2+21))  = log(2/23)
    //   log P(END|IG)   = log((1+1)/(2+21))  = log(2/23)
    vector<State> seq = { State::INTERGENIC, State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 2} };
    auto lp = Transition_Table::compute_log_probs(seq, chrs, 1.0);

    double denom_start = 1.0 + 1.0 * NUM_STATES;  // = 22
    double denom_ig    = 2.0 + 1.0 * NUM_STATES;  // = 23

    check(approx(lp[S(State::START)][S(State::INTERGENIC)], log(2.0 / denom_start)),
          "log_probs | Laplace: log P(IG|START) = log(2/22)");
    check(approx(lp[S(State::INTERGENIC)][S(State::INTERGENIC)], log(2.0 / denom_ig)),
          "log_probs | Laplace: log P(IG|IG) = log(2/23)");
    check(approx(lp[S(State::INTERGENIC)][S(State::END)], log(2.0 / denom_ig)),
          "log_probs | Laplace: log P(END|IG) = log(2/23)");
}

static void test_log_probs_no_neg_inf_with_smoothing() {
    // With alpha > 0, no entry should be -infinity (log(0)).
    vector<State> seq = { State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 1} };
    auto lp = Transition_Table::compute_log_probs(seq, chrs, 1.0);

    bool found_neg_inf = false;
    for (size_t i = 0; i < NUM_STATES; ++i)
        for (size_t j = 0; j < NUM_STATES; ++j)
            if (lp[i][j] == -numeric_limits<double>::infinity())
                { found_neg_inf = true; break; }

    check(!found_neg_inf, "log_probs | alpha > 0: no -infinity entries");
}

static void test_log_probs_unseen_transitions_get_neg_inf_alpha_zero() {
    // With alpha=0, an unseen transition must be log(0) = -inf.
    // seq: [INTERGENIC] — only START→IG and IG→END are seen.
    // So START→END is unseen → should be -inf.
    vector<State> seq = { State::INTERGENIC };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 1} };
    auto lp = Transition_Table::compute_log_probs(seq, chrs, 0.0);

    check(lp[S(State::START)][S(State::END)] == -numeric_limits<double>::infinity(),
          "log_probs | alpha=0: unseen START→END is -inf");
}

static void test_log_probs_three_chromosomes_normalization() {
    // Stress test: three chromosomes, varied states, alpha=0.5.
    vector<State> seq = {
        State::INTERGENIC, State::START_CODON_1,
        State::EXON_FRAME_1, State::EXON_FRAME_2,
        State::INTERGENIC
    };
    vector<Chromosome_Range> chrs = { {"chr1", 0, 2}, {"chr2", 2, 4}, {"chr3", 4, 5} };
    double alpha = 0.5;
    auto lp  = Transition_Table::compute_log_probs(seq, chrs, alpha);
    auto out = Transition_Table::count_outgoing(seq, chrs);

    check(rows_sum_to_one(lp, out, alpha),
          "log_probs | three chromosomes: rows normalize to 1.0 (alpha=0.5)");
}

// ─── main ────────────────────────────────────────────────────────────────────

int main() {
    cout << "\n=== count_bigrams ===\n";
    test_bigrams_single_position();
    test_bigrams_internal_transitions();
    test_bigrams_no_cross_chromosome_transition();
    test_bigrams_two_chromosomes_different_start_states();
    test_bigrams_repeated_transitions_accumulate();

    cout << "\n=== count_outgoing ===\n";
    test_outgoing_single_position();
    test_outgoing_start_count_equals_num_chromosomes();
    test_outgoing_consistency_with_bigrams();
    test_outgoing_consistency_single_chromosome();

    cout << "\n=== compute_log_probs ===\n";
    test_log_probs_values_are_non_positive();
    test_log_probs_rows_normalize_with_smoothing();
    test_log_probs_rows_normalize_small_alpha();
    test_log_probs_mle_values_alpha_zero();
    test_log_probs_laplace_smoothed_values();
    test_log_probs_no_neg_inf_with_smoothing();
    test_log_probs_unseen_transitions_get_neg_inf_alpha_zero();
    test_log_probs_three_chromosomes_normalization();

    cout << "\n─────────────────────────────────────\n";
    cout << "Results: " << g_pass << " passed, " << g_fail << " failed\n";
    return g_fail > 0 ? 1 : 0;
}
