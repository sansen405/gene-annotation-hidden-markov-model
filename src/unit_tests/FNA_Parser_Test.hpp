#pragma once

#include "../parsers/FNA_Parser.hpp"
#include "test_utils.hpp"
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace gene_hmm {

    using namespace std;

    static const string FNA_TEST_PATH = "src/unit_tests/tmp_fna_parser_test.fna";

    static void write_fna_fixture() {
        ofstream file(FNA_TEST_PATH);
        file << ">chr1 first chromosome\n";
        file << "ACgt\n";
        file << "TA\n";
        file << ">chr2 second chromosome\n";
        file << "ccGG\n";
    }

    static void cleanup_fna_fixture() {
        remove(FNA_TEST_PATH.c_str());
    }

    static void test_fna_parse_sequence() {
        cout << "\n[TEST 1] Parse sequence converts FASTA bases\n";

        vector<Nucleotide> seq = FNA_Parser::parse_sequence(FNA_TEST_PATH);
        vector<Nucleotide> expected = {
            Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T,
            Nucleotide::T, Nucleotide::A, Nucleotide::C, Nucleotide::C,
            Nucleotide::G, Nucleotide::G
        };

        CHECK("parsed sequence length matches fixture bases", seq.size() == expected.size());
        CHECK("parsed nucleotides preserve FASTA order and uppercase lowercase bases", seq == expected);
    }

    static void test_fna_sequence_length() {
        cout << "\n[TEST 2] Sequence length ignores headers\n";

        size_t length = FNA_Parser::get_sequence_length(FNA_TEST_PATH);
        CHECK("sequence length counts only base characters", length == 10);
    }

    static void test_fna_chromosome_offsets() {
        cout << "\n[TEST 3] Chromosome offsets use first header token\n";

        auto offsets = FNA_Parser::get_chromosome_offsets(FNA_TEST_PATH);
        CHECK("two chromosome offsets parsed", offsets.size() == 2);
        CHECK("chr1 starts at global offset 0", offsets["chr1"] == 0);
        CHECK("chr2 starts after chr1 bases", offsets["chr2"] == 6);
    }

    static void test_fna_chromosome_ranges() {
        cout << "\n[TEST 4] Chromosome ranges cover concatenated sequence\n";

        vector<Chromosome_Range> ranges = FNA_Parser::get_chromosome_ranges(FNA_TEST_PATH);
        CHECK("two chromosome ranges parsed", ranges.size() == 2);

        if (ranges.size() == 2) {
            CHECK("chr1 range is [0, 6)", ranges[0].name == "chr1" && ranges[0].start == 0 && ranges[0].end == 6);
            CHECK("chr2 range is [6, 10)", ranges[1].name == "chr2" && ranges[1].start == 6 && ranges[1].end == 10);
        }
    }

    static void test_fna_missing_file() {
        cout << "\n[TEST 5] Missing FASTA file throws\n";

        bool threw = false;
        try {
            FNA_Parser::get_sequence_length("src/unit_tests/does_not_exist.fna");
        } catch (const runtime_error&) {
            threw = true;
        }

        CHECK("missing FASTA path raises runtime_error", threw);
    }

    static void run_FNA_Parser_tests() {
        cout << "\nRunning FNA_Parser tests...\n";
        write_fna_fixture();

        test_fna_parse_sequence();
        test_fna_sequence_length();
        test_fna_chromosome_offsets();
        test_fna_chromosome_ranges();
        test_fna_missing_file();

        cleanup_fna_fixture();
        cout << "\nDone.\n";
    }

}
