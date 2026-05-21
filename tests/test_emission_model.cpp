#include "../src/model/Emission_Model.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gene_hmm;
using namespace std;

static int tests_passed = 0;
static int tests_failed = 0;

static void check(bool condition, const char* name){
    if(condition){
        cout << "  [PASS] " << name << "\n";
        tests_passed++;
    } else {
        cout << "  [FAIL] " << name << "\n";
        tests_failed++;
    }
}

static bool near(double a, double b, double tol = 1e-9){
    return fabs(a - b) < tol;
}

// ── helpers ──────────────────────────────────────────────────────────────────

static Chromosome_Range chr(size_t start, size_t end){
    return {"chr1", start, end};
}

// Encode a 5-mer manually (same formula as encode_5mer in Emission_Model.cpp)
static size_t encode_5mer_manual(const vector<Nucleotide>& nucs, size_t pos){
    size_t ctx = 0;
    for(size_t i = 0; i < 5; ++i)
        ctx = ctx * NUM_NUCLEOTIDES + (static_cast<size_t>(nucs[pos - 5 + i]) - 1);
    return ctx;
}

// ── count_markov1_emissions ───────────────────────────────────────────────────

static void test_markov1_count(){
    cout << "\ncount_markov1_emissions:\n";

    // (1) Empty chromosome range → all zeros
    {
        vector<State> states = {State::INTERGENIC};
        vector<Nucleotide> nucs = {Nucleotide::A};
        auto counts = Emission_Model::count_markov1_emissions(states, nucs, {}, {State::INTERGENIC});
        bool all_zero = true;
        for(auto& row : counts) for(auto v : row) if(v) all_zero = false;
        check(all_zero, "empty chromosome range → all zeros");
    }

    // (2) Single nucleotide (no prior) → no counts
    {
        vector<State> states = {State::INTERGENIC};
        vector<Nucleotide> nucs = {Nucleotide::A};
        vector<Chromosome_Range> chrs = {chr(0, 1)};
        auto counts = Emission_Model::count_markov1_emissions(states, nucs, chrs, {State::INTERGENIC});
        bool all_zero = true;
        for(auto& row : counts) for(auto v : row) if(v) all_zero = false;
        check(all_zero, "single nucleotide → no counts (no prior context)");
    }

    // (3) Two nucleotides A→C → one count at [A][C] = [0][1]
    {
        vector<State> states = {State::INTERGENIC, State::INTERGENIC};
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C};
        vector<Chromosome_Range> chrs = {chr(0, 2)};
        auto counts = Emission_Model::count_markov1_emissions(states, nucs, chrs, {State::INTERGENIC});
        check(counts[0][1] == 1, "A→C: counts[A][C] == 1");
        uint64_t total = 0;
        for(auto& row : counts) for(auto v : row) total += v;
        check(total == 1, "A→C: total count == 1");
    }

    // (4) Only target states are counted, non-target skipped
    {
        // pos 1 = EXON_FRAME_1 (not target), pos 2 = INTERGENIC (target)
        vector<State> states = {State::INTERGENIC, State::EXON_FRAME_1, State::INTERGENIC};
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G};
        vector<Chromosome_Range> chrs = {chr(0, 3)};
        auto counts = Emission_Model::count_markov1_emissions(states, nucs, chrs, {State::INTERGENIC});
        // pos 2 (INTERGENIC, target): context=nucs[1]=C (idx 1), emit=nucs[2]=G (idx 2)
        check(counts[1][2] == 1, "non-target skipped: only INTERGENIC pos counted");
        uint64_t total = 0;
        for(auto& row : counts) for(auto v : row) total += v;
        check(total == 1, "non-target skipped: total == 1");
    }

    // (5) Chromosome boundary: context does not cross chromosomes
    {
        // chr1: positions 0–1, chr2: positions 2–3
        // nucs[1]→nucs[2] would cross boundary → must not count
        vector<State> states = {State::INTERGENIC, State::INTERGENIC,
                                 State::INTERGENIC, State::INTERGENIC};
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T};
        vector<Chromosome_Range> chrs = {{"chr1", 0, 2}, {"chr2", 2, 4}};
        auto counts = Emission_Model::count_markov1_emissions(states, nucs, chrs, {State::INTERGENIC});
        // chr1 pos 1: context=A(0), emit=C(1) → counts[0][1]++
        // chr2 pos 3: context=G(2), emit=T(3) → counts[2][3]++
        // cross-boundary nucs[1]→nucs[2] must NOT appear
        check(counts[0][1] == 1, "boundary: chr1 A→C counted");
        check(counts[2][3] == 1, "boundary: chr2 G→T counted");
        uint64_t total = 0;
        for(auto& row : counts) for(auto v : row) total += v;
        check(total == 2, "boundary: cross-chromosome pair not counted");
    }

    // (6) Multiple target states pooled
    {
        vector<State> states = {State::INTRON_1, State::INTRON_2};
        vector<Nucleotide> nucs = {Nucleotide::G, Nucleotide::T};
        vector<Chromosome_Range> chrs = {chr(0, 2)};
        auto counts = Emission_Model::count_markov1_emissions(
            states, nucs, chrs, {State::INTRON_1, State::INTRON_2});
        // pos 1 (INTRON_2, target): context=G(2), emit=T(3)
        check(counts[2][3] == 1, "multiple target states: INTRON_2 position counted");
    }

    // (7) Correct accumulation over many positions
    {
        // AAAAAA → 5 transitions all A→A: counts[0][0] == 5
        vector<State> states(6, State::INTERGENIC);
        vector<Nucleotide> nucs(6, Nucleotide::A);
        vector<Chromosome_Range> chrs = {chr(0, 6)};
        auto counts = Emission_Model::count_markov1_emissions(states, nucs, chrs, {State::INTERGENIC});
        check(counts[0][0] == 5, "6 A's → 5 A→A transitions");
    }
}

// ── count_markov5_emissions ───────────────────────────────────────────────────

static void test_markov5_count(){
    cout << "\ncount_markov5_emissions:\n";

    // (1) Empty range → all zeros
    {
        vector<State> states = {State::EXON_FRAME_1};
        vector<Nucleotide> nucs = {Nucleotide::A};
        auto counts = Emission_Model::count_markov5_emissions(states, nucs, {}, {State::EXON_FRAME_1});
        bool all_zero = true;
        for(auto& row : counts) for(auto v : row) if(v) all_zero = false;
        check(all_zero, "empty range → all zeros");
    }

    // (2) Fewer than 6 positions → no counts (need 5 context + 1 emit)
    {
        vector<State> states(5, State::EXON_FRAME_1);
        vector<Nucleotide> nucs(5, Nucleotide::A);
        vector<Chromosome_Range> chrs = {chr(0, 5)};
        auto counts = Emission_Model::count_markov5_emissions(states, nucs, chrs, {State::EXON_FRAME_1});
        bool all_zero = true;
        for(auto& row : counts) for(auto v : row) if(v) all_zero = false;
        check(all_zero, "5 positions → no counts (need at least 6)");
    }

    // (3) Exactly 6 positions → exactly one count
    {
        // AAAAAC: context=AAAAA→0, emit=C(1)
        vector<State> states(6, State::EXON_FRAME_1);
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::A, Nucleotide::A,
                                    Nucleotide::A, Nucleotide::A, Nucleotide::C};
        vector<Chromosome_Range> chrs = {chr(0, 6)};
        auto counts = Emission_Model::count_markov5_emissions(states, nucs, chrs, {State::EXON_FRAME_1});
        uint64_t total = 0;
        for(auto& row : counts) for(auto v : row) total += v;
        check(total == 1, "6 positions → exactly one count");
        // AAAAA in base-4 = 0, emit C = idx 1
        check(counts[0][1] == 1, "AAAAAC: context idx 0 (AAAAA), emit C idx 1");
    }

    // (4) Context encoding: verify specific 5-mer index
    {
        // context = ACGT A → emit = T
        // A=0,C=1,G=2,T=3 in 0-indexed
        // index = 0*256 + 1*64 + 2*16 + 3*4 + 0 = 0+64+32+12+0 = 108
        vector<State> states(7, State::EXON_FRAME_1);
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T,
                                    Nucleotide::A, Nucleotide::T, Nucleotide::A};
        vector<Chromosome_Range> chrs = {chr(0, 7)};
        auto counts = Emission_Model::count_markov5_emissions(states, nucs, chrs, {State::EXON_FRAME_1});
        // pos 5: context = nucs[0..4] = ACGTA → check manually
        size_t expected_ctx = encode_5mer_manual(nucs, 5);
        check(counts[expected_ctx][3] == 1, "ACGTA context: emit T at correct index");
    }

    // (5) Only target states counted
    {
        // pos 5 = INTERGENIC (not target) → no count
        vector<State> states = {State::EXON_FRAME_1, State::EXON_FRAME_1, State::EXON_FRAME_1,
                                  State::EXON_FRAME_1, State::EXON_FRAME_1, State::INTERGENIC};
        vector<Nucleotide> nucs(6, Nucleotide::A);
        nucs[5] = Nucleotide::C;
        vector<Chromosome_Range> chrs = {chr(0, 6)};
        auto counts = Emission_Model::count_markov5_emissions(states, nucs, chrs, {State::EXON_FRAME_1});
        bool all_zero = true;
        for(auto& row : counts) for(auto v : row) if(v) all_zero = false;
        check(all_zero, "emit position is non-target → no count");
    }

    // (6) Chromosome boundary: chr.start+5 begins iteration, no cross-boundary context
    {
        // chr1: 0–5, chr2: 6–11. Position 6 (chr2.start+5=11 minimum, so pos=11)
        // chr2 has 6 positions (6..11), so first valid pos = 6+5=11 (last position)
        vector<State> states(12, State::EXON_FRAME_1);
        vector<Nucleotide> nucs(12, Nucleotide::A);
        nucs[11] = Nucleotide::C;
        vector<Chromosome_Range> chrs = {{"chr1", 0, 6}, {"chr2", 6, 12}};
        auto counts = Emission_Model::count_markov5_emissions(states, nucs, chrs, {State::EXON_FRAME_1});
        // chr1: pos 5: context=AAAAA(0), emit=A(0) → counts[0][0]++
        // chr2: pos 11: context=AAAAA(0), emit=C(1) → counts[0][1]++
        // No cross-boundary: nucs[1..5] should not be used as context for pos 6
        check(counts[0][0] == 1, "boundary: chr1 emits one count (AAAAA→A)");
        check(counts[0][1] == 1, "boundary: chr2 emits one count (AAAAA→C)");
        uint64_t total = 0;
        for(auto& row : counts) for(auto v : row) total += v;
        check(total == 2, "boundary: total == 2, no cross-chromosome context used");
    }

    // (7) frame_tied: pooling EXON_FRAME_1/2/3 into one table
    {
        // pos5=EXON_FRAME_1, pos6=EXON_FRAME_2, pos7=EXON_FRAME_3 all contribute
        vector<State> states = {State::EXON_FRAME_1, State::EXON_FRAME_1, State::EXON_FRAME_1,
                                  State::EXON_FRAME_1, State::EXON_FRAME_1,
                                  State::EXON_FRAME_1, State::EXON_FRAME_2, State::EXON_FRAME_3};
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::A, Nucleotide::A,
                                    Nucleotide::A, Nucleotide::A,
                                    Nucleotide::C, Nucleotide::G, Nucleotide::T};
        vector<Chromosome_Range> chrs = {chr(0, 8)};
        auto counts = Emission_Model::count_markov5_emissions(
            states, nucs, chrs, {State::EXON_FRAME_1, State::EXON_FRAME_2, State::EXON_FRAME_3});
        uint64_t total = 0;
        for(auto& row : counts) for(auto v : row) total += v;
        check(total == 3, "frame_tied: all 3 frame states contribute counts");
    }
}

// ── count_pssm_emissions ─────────────────────────────────────────────────────

static void test_pssm_count(){
    cout << "\ncount_pssm_emissions:\n";

    // (1) Target state at valid position with window_left=2, window_right=2 → 4 columns counted
    {
        // chr: positions 0–5. State at pos=2: window = [0,4)
        vector<State> states(6, State::INTERGENIC);
        states[2] = State::DONOR_1;
        // nucs[0]=A, [1]=C, [2]=G, [3]=T, [4]=A, [5]=C
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G,
                                    Nucleotide::T, Nucleotide::A, Nucleotide::C};
        vector<Chromosome_Range> chrs = {chr(0, 6)};
        auto counts = Emission_Model::count_pssm_emissions(
            states, nucs, chrs, {State::DONOR_1}, 2, 2);
        // window: nucs[0]=A(0), nucs[1]=C(1), nucs[2]=G(2), nucs[3]=T(3)
        check(counts.size() == 4, "pssm: window_size == 4");
        check(counts[0][0] == 1, "pssm: col 0 = A");
        check(counts[1][1] == 1, "pssm: col 1 = C");
        check(counts[2][2] == 1, "pssm: col 2 = G");
        check(counts[3][3] == 1, "pssm: col 3 = T");
    }

    // (2) Target too close to chromosome start → skipped
    {
        // window_left=3, pos=1 → pos < chr.start + window_left (1 < 3) → skip
        vector<State> states = {State::INTERGENIC, State::DONOR_1, State::INTERGENIC,
                                  State::INTERGENIC, State::INTERGENIC, State::INTERGENIC};
        vector<Nucleotide> nucs(6, Nucleotide::A);
        vector<Chromosome_Range> chrs = {chr(0, 6)};
        auto counts = Emission_Model::count_pssm_emissions(
            states, nucs, chrs, {State::DONOR_1}, 3, 2);
        bool all_zero = true;
        for(auto& col : counts) for(auto v : col) if(v) all_zero = false;
        check(all_zero, "pssm: target too close to chr start → skipped");
    }

    // (3) Target too close to chromosome end → skipped
    {
        // window_right=3, chr.end=6, pos=4 → pos+3=7 > 6 → skip
        vector<State> states(6, State::INTERGENIC);
        states[4] = State::DONOR_1;
        vector<Nucleotide> nucs(6, Nucleotide::A);
        vector<Chromosome_Range> chrs = {chr(0, 6)};
        auto counts = Emission_Model::count_pssm_emissions(
            states, nucs, chrs, {State::DONOR_1}, 2, 3);
        bool all_zero = true;
        for(auto& col : counts) for(auto v : col) if(v) all_zero = false;
        check(all_zero, "pssm: target too close to chr end → skipped");
    }

    // (4) Exactly at boundary (pos = chr.start + window_left) → counted
    {
        // window_left=2, pos=2 (== chr.start+2), window_right=2, chr.end=6 (pos+2=4 <= 6) → ok
        vector<State> states(6, State::INTERGENIC);
        states[2] = State::DONOR_1;
        vector<Nucleotide> nucs(6, Nucleotide::T);
        vector<Chromosome_Range> chrs = {chr(0, 6)};
        auto counts = Emission_Model::count_pssm_emissions(
            states, nucs, chrs, {State::DONOR_1}, 2, 2);
        uint64_t total = 0;
        for(auto& col : counts) for(auto v : col) total += v;
        check(total == 4, "pssm: exactly at boundary → 4 columns counted");
    }

    // (5) Chromosome boundary: pos in chr2 that would need chr1 context → skipped
    {
        // chr1: 0–4, chr2: 4–8. State at pos=5 (chr2), window_left=3 → needs pos 2,3,4
        // pos=5, chr2.start=4, chr2.start+window_left=7 → 5 < 7 → skip
        vector<State> states(8, State::INTERGENIC);
        states[5] = State::DONOR_1;
        vector<Nucleotide> nucs(8, Nucleotide::A);
        vector<Chromosome_Range> chrs = {{"chr1", 0, 4}, {"chr2", 4, 8}};
        auto counts = Emission_Model::count_pssm_emissions(
            states, nucs, chrs, {State::DONOR_1}, 3, 2);
        bool all_zero = true;
        for(auto& col : counts) for(auto v : col) if(v) all_zero = false;
        check(all_zero, "pssm: window would cross chromosome boundary → skipped");
    }

    // (6) Multiple target positions → counts accumulate across columns
    {
        // Two DONOR_1 positions at pos=2 and pos=4, window_left=1, window_right=1
        vector<State> states = {State::INTERGENIC, State::INTERGENIC, State::DONOR_1,
                                  State::INTERGENIC, State::DONOR_1, State::INTERGENIC};
        // nucs: A C G T A C
        vector<Nucleotide> nucs = {Nucleotide::A, Nucleotide::C, Nucleotide::G,
                                    Nucleotide::T, Nucleotide::A, Nucleotide::C};
        vector<Chromosome_Range> chrs = {chr(0, 6)};
        auto counts = Emission_Model::count_pssm_emissions(
            states, nucs, chrs, {State::DONOR_1}, 1, 1);
        // pos=2: window nucs[1]=C(1), nucs[2]=G(2) → col0[C]++, col1[G]++
        // pos=4: window nucs[3]=T(3), nucs[4]=A(0) → col0[T]++, col1[A]++
        check(counts[0][1] == 1, "pssm multi: col0 C count");
        check(counts[1][2] == 1, "pssm multi: col1 G count");
        check(counts[0][3] == 1, "pssm multi: col0 T count");
        check(counts[1][0] == 1, "pssm multi: col1 A count");
    }
}

// ── compute_markov1_log_probs ─────────────────────────────────────────────────

static void test_markov1_log_probs(){
    cout << "\ncompute_markov1_log_probs:\n";

    // (1) All zeros + alpha=1 → uniform log(0.25) for all entries
    {
        Emission_Model::Markov1_Count counts = {};
        auto lp = Emission_Model::compute_markov1_log_probs(counts, 1.0);
        bool all_uniform = true;
        double expected = log(1.0 / NUM_NUCLEOTIDES);
        for(auto& row : lp)
            for(auto v : row)
                if(!near(v, expected)) all_uniform = false;
        check(all_uniform, "zero counts + alpha=1 → uniform log(0.25)");
    }

    // (2) Each row sums to 1 in probability space (exp then sum ≈ 1)
    {
        Emission_Model::Markov1_Count counts = {};
        counts[0][0] = 10; counts[0][1] = 5; counts[0][2] = 3; counts[0][3] = 2;
        auto lp = Emission_Model::compute_markov1_log_probs(counts, 0.5);
        bool rows_sum_to_one = true;
        for(auto& row : lp){
            double sum = 0;
            for(auto v : row) sum += exp(v);
            if(!near(sum, 1.0, 1e-9)) rows_sum_to_one = false;
        }
        check(rows_sum_to_one, "all rows sum to 1.0 in probability space");
    }

    // (3) Specific count → verify exact log prob
    {
        Emission_Model::Markov1_Count counts = {};
        // row 0 (context A): 3 A, 0 C, 0 G, 0 T. alpha=1
        counts[0][0] = 3;
        auto lp = Emission_Model::compute_markov1_log_probs(counts, 1.0);
        double expected = log((3.0 + 1.0) / (3.0 + 1.0 * NUM_NUCLEOTIDES));
        check(near(lp[0][0], expected), "specific count: log P(A|A) correct");
    }

    // (4) Alpha=0 and all identical counts → still normalized
    {
        Emission_Model::Markov1_Count counts = {};
        for(auto& row : counts) for(auto& v : row) v = 10;
        auto lp = Emission_Model::compute_markov1_log_probs(counts, 0.0);
        double expected = log(0.25);
        bool all_eq = true;
        for(auto& row : lp) for(auto v : row) if(!near(v, expected)) all_eq = false;
        check(all_eq, "alpha=0, uniform counts → log(0.25) everywhere");
    }
}

// ── compute_markov5_log_probs ─────────────────────────────────────────────────

static void test_markov5_log_probs(){
    cout << "\ncompute_markov5_log_probs:\n";

    // (1) All zeros + alpha=1 → uniform for all 1024 rows
    {
        Emission_Model::Markov5_Count counts = {};
        auto lp = Emission_Model::compute_markov5_log_probs(counts, 1.0);
        double expected = log(1.0 / NUM_NUCLEOTIDES);
        bool all_uniform = true;
        for(auto& row : lp)
            for(auto v : row)
                if(!near(v, expected)) all_uniform = false;
        check(all_uniform, "zero counts + alpha=1 → uniform across all 1024 contexts");
    }

    // (2) All 1024 rows normalized (sum to 1 in probability space)
    {
        Emission_Model::Markov5_Count counts = {};
        // Give some varied counts to row 108 (ACGTA)
        counts[108][0] = 5; counts[108][1] = 3; counts[108][2] = 2; counts[108][3] = 1;
        auto lp = Emission_Model::compute_markov5_log_probs(counts, 0.5);
        bool all_normalized = true;
        for(size_t ctx = 0; ctx < Emission_Model::MARKOV5_CONTEXTS; ++ctx){
            double sum = 0;
            for(auto v : lp[ctx]) sum += exp(v);
            if(!near(sum, 1.0, 1e-9)) all_normalized = false;
        }
        check(all_normalized, "all 1024 rows sum to 1.0 in probability space");
    }

    // (3) Specific row with known counts → verify log prob
    {
        Emission_Model::Markov5_Count counts = {};
        // row 0 (AAAAA): 100 A, 0 others. alpha=0.5
        counts[0][0] = 100;
        auto lp = Emission_Model::compute_markov5_log_probs(counts, 0.5);
        double expected = log((100.0 + 0.5) / (100.0 + 0.5 * NUM_NUCLEOTIDES));
        check(near(lp[0][0], expected), "specific row: log P(A|AAAAA) correct");
    }

    // (4) Unseen contexts (zero counts) still get valid log probs via smoothing
    {
        Emission_Model::Markov5_Count counts = {};
        counts[0][0] = 100; // only row 0 touched
        auto lp = Emission_Model::compute_markov5_log_probs(counts, 1.0);
        // All other rows should be uniform log(0.25)
        bool unseen_uniform = true;
        double expected = log(0.25);
        for(size_t ctx = 1; ctx < Emission_Model::MARKOV5_CONTEXTS; ++ctx)
            for(auto v : lp[ctx])
                if(!near(v, expected)) unseen_uniform = false;
        check(unseen_uniform, "unseen contexts get uniform log probs via alpha smoothing");
    }
}

// ── compute_pssm_log_probs ────────────────────────────────────────────────────

static void test_pssm_log_probs(){
    cout << "\ncompute_pssm_log_probs:\n";

    // (1) All zeros + alpha=1 → uniform log(0.25) per column
    {
        Emission_Model::PSSM_Count counts(9, array<uint64_t, NUM_NUCLEOTIDES>{});
        auto lp = Emission_Model::compute_pssm_log_probs(counts, 1.0);
        double expected = log(0.25);
        bool all_uniform = true;
        for(auto& col : lp) for(auto v : col) if(!near(v, expected)) all_uniform = false;
        check(all_uniform, "zero counts + alpha=1 → uniform log(0.25) per column");
    }

    // (2) Each column sums to 1 in probability space
    {
        Emission_Model::PSSM_Count counts(9, array<uint64_t, NUM_NUCLEOTIDES>{});
        counts[0][0] = 10; counts[0][1] = 5;
        counts[3][2] = 7; counts[3][3] = 3;
        auto lp = Emission_Model::compute_pssm_log_probs(counts, 0.5);
        bool all_normalized = true;
        for(auto& col : lp){
            double sum = 0;
            for(auto v : col) sum += exp(v);
            if(!near(sum, 1.0, 1e-9)) all_normalized = false;
        }
        check(all_normalized, "all columns sum to 1.0 in probability space");
    }

    // (3) Specific column with known counts → verify exact value
    {
        Emission_Model::PSSM_Count counts(4, array<uint64_t, NUM_NUCLEOTIDES>{});
        counts[2][0] = 8; // col 2: 8 A, rest 0. alpha=1
        auto lp = Emission_Model::compute_pssm_log_probs(counts, 1.0);
        double expected = log((8.0 + 1.0) / (8.0 + 1.0 * NUM_NUCLEOTIDES));
        check(near(lp[2][0], expected), "specific column: log P(A|col2) correct");
    }

    // (4) Columns are independent: count in col 0 doesn't affect col 1
    {
        Emission_Model::PSSM_Count counts(2, array<uint64_t, NUM_NUCLEOTIDES>{});
        counts[0][0] = 100; // col 0: all A
        // col 1: all zero
        auto lp = Emission_Model::compute_pssm_log_probs(counts, 1.0);
        double col1_expected = log(0.25);
        bool col1_uniform = true;
        for(auto v : lp[1]) if(!near(v, col1_expected)) col1_uniform = false;
        check(col1_uniform, "column independence: col1 still uniform despite col0 counts");
    }

    // (5) Window size preserved in output
    {
        Emission_Model::PSSM_Count counts(18, array<uint64_t, NUM_NUCLEOTIDES>{});
        auto lp = Emission_Model::compute_pssm_log_probs(counts, 0.5);
        check(lp.size() == 18, "PSSM output preserves window size (18 for acceptor)");
    }
}

// ── get_deterministic_log_prob ────────────────────────────────────────────────

static void test_deterministic_log_prob(){
    cout << "\nget_deterministic_log_prob:\n";

    // (1) START_CODON_1 = A (first base of ATG)
    check(near(Emission_Model::get_deterministic_log_prob(State::START_CODON_1, Nucleotide::A), 0.0),
          "START_CODON_1 + A → 0.0");
    check(Emission_Model::get_deterministic_log_prob(State::START_CODON_1, Nucleotide::C) == LOG_ZERO,
          "START_CODON_1 + C → LOG_ZERO");
    check(Emission_Model::get_deterministic_log_prob(State::START_CODON_1, Nucleotide::G) == LOG_ZERO,
          "START_CODON_1 + G → LOG_ZERO");
    check(Emission_Model::get_deterministic_log_prob(State::START_CODON_1, Nucleotide::T) == LOG_ZERO,
          "START_CODON_1 + T → LOG_ZERO");

    // (2) START_CODON_2 = T
    check(near(Emission_Model::get_deterministic_log_prob(State::START_CODON_2, Nucleotide::T), 0.0),
          "START_CODON_2 + T → 0.0");
    check(Emission_Model::get_deterministic_log_prob(State::START_CODON_2, Nucleotide::A) == LOG_ZERO,
          "START_CODON_2 + A → LOG_ZERO");

    // (3) START_CODON_3 = G
    check(near(Emission_Model::get_deterministic_log_prob(State::START_CODON_3, Nucleotide::G), 0.0),
          "START_CODON_3 + G → 0.0");
    check(Emission_Model::get_deterministic_log_prob(State::START_CODON_3, Nucleotide::A) == LOG_ZERO,
          "START_CODON_3 + A → LOG_ZERO");

    // (4) STOP_CODON_1 = T (all stop codons TAA/TAG/TGA start with T)
    check(near(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_1, Nucleotide::T), 0.0),
          "STOP_CODON_1 + T → 0.0");
    check(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_1, Nucleotide::A) == LOG_ZERO,
          "STOP_CODON_1 + A → LOG_ZERO");

    // (5) STOP_CODON_2 = {A, G} with log(0.5) each (TAA→A, TAG→A, TGA→G)
    check(near(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_2, Nucleotide::A), log(0.5)),
          "STOP_CODON_2 + A → log(0.5)");
    check(near(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_2, Nucleotide::G), log(0.5)),
          "STOP_CODON_2 + G → log(0.5)");
    check(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_2, Nucleotide::C) == LOG_ZERO,
          "STOP_CODON_2 + C → LOG_ZERO");
    check(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_2, Nucleotide::T) == LOG_ZERO,
          "STOP_CODON_2 + T → LOG_ZERO");

    // (6) STOP_CODON_3 = {A, G} with log(0.5) each (TAA→A, TAG→G, TGA→A)
    check(near(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_3, Nucleotide::A), log(0.5)),
          "STOP_CODON_3 + A → log(0.5)");
    check(near(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_3, Nucleotide::G), log(0.5)),
          "STOP_CODON_3 + G → log(0.5)");
    check(Emission_Model::get_deterministic_log_prob(State::STOP_CODON_3, Nucleotide::C) == LOG_ZERO,
          "STOP_CODON_3 + C → LOG_ZERO");

    // (7) Non-deterministic states return 0.0 (log-prob no-op)
    check(near(Emission_Model::get_deterministic_log_prob(State::INTERGENIC, Nucleotide::A), 0.0),
          "INTERGENIC → 0.0 (no-op)");
    check(near(Emission_Model::get_deterministic_log_prob(State::EXON_FRAME_1, Nucleotide::C), 0.0),
          "EXON_FRAME_1 → 0.0 (no-op)");
    check(near(Emission_Model::get_deterministic_log_prob(State::INTRON_2, Nucleotide::T), 0.0),
          "INTRON_2 → 0.0 (no-op)");
    check(near(Emission_Model::get_deterministic_log_prob(State::DONOR_1, Nucleotide::G), 0.0),
          "DONOR_1 → 0.0 (no-op, emission handled via PSSM transition bonus)");
}

// ─────────────────────────────────────────────────────────────────────────────

int main(){
    cout << "=== Emission_Model Tests ===\n";
    test_markov1_count();
    test_markov5_count();
    test_pssm_count();
    test_markov1_log_probs();
    test_markov5_log_probs();
    test_pssm_log_probs();
    test_deterministic_log_prob();

    cout << "\n--- Results ---\n";
    cout << "Passed: " << tests_passed << "\n";
    cout << "Failed: " << tests_failed << "\n";
    return tests_failed == 0 ? 0 : 1;
}
