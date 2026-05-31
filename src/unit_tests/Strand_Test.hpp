#pragma once

#include "../parsers/FNA_Parser.hpp"
#include "../parsers/GFF_Parser.hpp"
#include "../topology/Topology.hpp"
#include "Test_Utils.hpp"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace gene_hmm {

    using namespace std;

    static const string STRAND_FNA_PATH = "src/unit_tests/tmp_strand_test.fna";
    static const string STRAND_GFF_PATH = "src/unit_tests/tmp_strand_test.gff";

    static void write_strand_fixture() {
        ofstream fna(STRAND_FNA_PATH);
        fna << ">chr1 strand fixture\n";
        fna << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

        ofstream gff(STRAND_GFF_PATH);
        gff << "chr1\ttest\tCDS\t11\t19\t.\t-\t.\tParent=gene1\n";
        gff << "chr1\ttest\tCDS\t42\t50\t.\t-\t.\tParent=gene1\n";
    }

    static void cleanup_strand_fixture() {
        remove(STRAND_FNA_PATH.c_str());
        remove(STRAND_GFF_PATH.c_str());
    }

    static void test_reverse_complement() {
        cout << "\n[TEST 1] reverse_complement correctness and involution\n";

        vector<Nucleotide> seq = {
            Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T, Nucleotide::A};
        vector<Nucleotide> expected = {
            Nucleotide::T, Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T};

        vector<Nucleotide> rc = FNA_Parser::reverse_complement(seq);
        CHECK("reverse_complement(ACGTA) == TACGT", rc == expected);

        vector<Nucleotide> twice = FNA_Parser::reverse_complement(rc);
        CHECK("reverse_complement is its own inverse", twice == seq);
    }

    static size_t count_state(const vector<State>& states, State target) {
        return static_cast<size_t>(count(states.begin(), states.end(), target));
    }

    static void test_minus_strand_labeling() {
        cout << "\n[TEST 2] Minus-strand gene labels in revcomp coordinates\n";

        string gff = STRAND_GFF_PATH;
        string fna = STRAND_FNA_PATH;
        vector<int> regions = GFF_Parser::parse_regions(gff, fna, '-');
        vector<State> states = GFF_Parser::parse_states(regions);

        CHECK("state vector covers all 60 bases", states.size() == 60);
        CHECK("exactly one start codon (START_CODON_1)", count_state(states, State::START_CODON_1) == 1);
        CHECK("exactly one stop codon (STOP_CODON_3)", count_state(states, State::STOP_CODON_3) == 1);
        CHECK("one donor and one acceptor",
              count_state(states, State::DONOR_1) +
              count_state(states, State::DONOR_2) +
              count_state(states, State::DONOR_3) == 1 &&
              count_state(states, State::ACCEPTOR_1) +
              count_state(states, State::ACCEPTOR_2) +
              count_state(states, State::ACCEPTOR_3) == 1);

        CHECK("start codon at revcomp index 10", states[10] == State::START_CODON_1);
        CHECK("stop codon closes at revcomp index 49", states[49] == State::STOP_CODON_3);

        bool legal = true;
        for (size_t i = 0; i + 1 < states.size(); ++i) {
            const auto& allowed = Transitions.at(states[i]);
            if (find(allowed.begin(), allowed.end(), states[i + 1]) == allowed.end()) {
                legal = false;
                break;
            }
        }
        CHECK("all minus-strand transitions are legal", legal);
    }

    static void test_plus_strand_ignores_minus_gene() {
        cout << "\n[TEST 3] Plus-strand track ignores the minus-strand gene\n";

        string gff = STRAND_GFF_PATH;
        string fna = STRAND_FNA_PATH;
        vector<int> regions = GFF_Parser::parse_regions(gff, fna, '+');
        vector<State> states = GFF_Parser::parse_states(regions);

        CHECK("plus track has no start codons", count_state(states, State::START_CODON_1) == 0);
        CHECK("plus track has no CDS regions", count(regions.begin(), regions.end(), GFF_Parser::CDS_REGION) == 0);
    }

    static void run_Strand_tests() {
        cout << "\nRunning Strand (minus-strand) tests...\n";
        write_strand_fixture();

        test_reverse_complement();
        test_minus_strand_labeling();
        test_plus_strand_ignores_minus_gene();

        cleanup_strand_fixture();
        cout << "\nDone.\n";
    }

}
