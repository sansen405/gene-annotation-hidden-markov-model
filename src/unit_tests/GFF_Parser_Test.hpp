#pragma once

#include "../parsers/FNA_Parser.hpp"
#include "../parsers/GFF_Parser.hpp"
#include "../topology/Topology.hpp"
#include "test_utils.hpp"
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

namespace gene_hmm {

    using namespace std;

    static string state_name(State s) {
        switch (s) {
            case State::START:          return "START";
            case State::INTERGENIC:     return "INTERGENIC";
            case State::START_CODON_1:  return "START_CODON_1";
            case State::START_CODON_2:  return "START_CODON_2";
            case State::START_CODON_3:  return "START_CODON_3";
            case State::EXON_FRAME_1:   return "EXON_FRAME_1";
            case State::EXON_FRAME_2:   return "EXON_FRAME_2";
            case State::EXON_FRAME_3:   return "EXON_FRAME_3";
            case State::DONOR_1:        return "DONOR_1";
            case State::DONOR_2:        return "DONOR_2";
            case State::DONOR_3:        return "DONOR_3";
            case State::INTRON_1:       return "INTRON_1";
            case State::INTRON_2:       return "INTRON_2";
            case State::INTRON_3:       return "INTRON_3";
            case State::ACCEPTOR_1:     return "ACCEPTOR_1";
            case State::ACCEPTOR_2:     return "ACCEPTOR_2";
            case State::ACCEPTOR_3:     return "ACCEPTOR_3";
            case State::STOP_CODON_1:   return "STOP_CODON_1";
            case State::STOP_CODON_2:   return "STOP_CODON_2";
            case State::STOP_CODON_3:   return "STOP_CODON_3";
            case State::END:            return "END";
            default:                    return "UNKNOWN";
        }
    }

    static void test_sequence_length(const string& fna, const vector<Nucleotide>& nuc_seq) {
        cout << "\n[TEST 1] Sequence length sanity\n";
        size_t reported = FNA_Parser::get_sequence_length(fna);
        CHECK("parse_sequence length matches get_sequence_length", nuc_seq.size() == reported);
        CHECK("sequence is non-empty", !nuc_seq.empty());
        cout << "  total bases: " << nuc_seq.size() << "\n";
    }

    static void test_regions(const vector<int>& regions, const vector<Nucleotide>& nuc_seq) {
        cout << "\n[TEST 2] Region vector coverage\n";
        CHECK("regions length == nuc_seq length", regions.size() == nuc_seq.size());

        size_t n_intergenic = count(regions.begin(), regions.end(), 0);
        size_t n_cds        = count(regions.begin(), regions.end(), 1);
        size_t n_intron     = count(regions.begin(), regions.end(), 2);
        size_t unexpected   = regions.size() - n_intergenic - n_cds - n_intron;

        CHECK("no unexpected region codes", unexpected == 0);
        CHECK("at least some CDS bases present", n_cds > 0);
        CHECK("at least some intergenic bases present", n_intergenic > 0);

        cout << "  intergenic: " << n_intergenic
             << "  CDS: " << n_cds
             << "  intron: " << n_intron << "\n";
    }

    static void test_states_basic(const vector<State>& states, const vector<Nucleotide>& nuc_seq) {
        cout << "\n[TEST 3] State vector basic sanity\n";
        CHECK("states length == nuc_seq length", states.size() == nuc_seq.size());

        size_t n_intergenic = 0, n_exon = 0, n_intron = 0,
               n_start = 0, n_stop = 0, n_donor = 0, n_acceptor = 0;
        for (State s : states) {
            switch (s) {
                case State::INTERGENIC:                               n_intergenic++; break;
                case State::EXON_FRAME_1:
                case State::EXON_FRAME_2:
                case State::EXON_FRAME_3:                            n_exon++;       break;
                case State::INTRON_1:
                case State::INTRON_2:
                case State::INTRON_3:                                n_intron++;     break;
                case State::START_CODON_1:
                case State::START_CODON_2:
                case State::START_CODON_3:                           n_start++;      break;
                case State::STOP_CODON_1:
                case State::STOP_CODON_2:
                case State::STOP_CODON_3:                            n_stop++;       break;
                case State::DONOR_1:
                case State::DONOR_2:
                case State::DONOR_3:                                 n_donor++;      break;
                case State::ACCEPTOR_1:
                case State::ACCEPTOR_2:
                case State::ACCEPTOR_3:                              n_acceptor++;   break;
                default: break;
            }
        }

        CHECK("start codon bases are a multiple of 3", n_start % 3 == 0);
        CHECK("stop codon bases are a multiple of 3",  n_stop  % 3 == 0);
        CHECK("equal number of start and stop codon bases", n_start == n_stop);
        CHECK("at least some gene states present", n_exon > 0);

        cout << "  intergenic: " << n_intergenic
             << "  exon: " << n_exon
             << "  intron: " << n_intron
             << "  donor: " << n_donor
             << "  acceptor: " << n_acceptor
             << "  start codons: " << n_start
             << "  stop codons: " << n_stop << "\n";
    }

    static void test_all_transitions_legal(const vector<State>& states) {
        cout << "\n[TEST 4] All consecutive transitions are legal\n";

        size_t illegal_count = 0;
        size_t first_illegal_pos = 0;
        State  first_from = State::INTERGENIC, first_to = State::INTERGENIC;

        for (size_t i = 0; i + 1 < states.size(); ++i) {
            const auto& allowed = Transitions.at(states[i]);
            if (find(allowed.begin(), allowed.end(), states[i + 1]) == allowed.end()) {
                if (illegal_count == 0) {
                    first_illegal_pos = i;
                    first_from = states[i];
                    first_to   = states[i + 1];
                }
                illegal_count++;
            }
        }

        if (illegal_count == 0) {
            ut_pass("zero illegal transitions in full state sequence");
        } else {
            ut_fail("zero illegal transitions in full state sequence");
            cout << "  " << illegal_count << " illegal transition(s) found\n";
            cout << "  first at pos " << first_illegal_pos
                 << ": " << state_name(first_from)
                 << " -> " << state_name(first_to) << "\n";
        }
    }

    static void test_gene_structure(const vector<State>& states) {
        cout << "\n[TEST 5] Per-gene structural checks\n";

        size_t total_genes    = 0;
        size_t bad_start      = 0;
        size_t bad_stop       = 0;
        size_t bad_frame_seq  = 0;
        size_t donor_acceptor_mismatch = 0;

        size_t i = 0;
        while (i < states.size()) {
            if (states[i] == State::INTERGENIC) { ++i; continue; }

            size_t gene_start = i;
            while (i < states.size() && states[i] != State::INTERGENIC) ++i;
            size_t gene_end = i;

            total_genes++;
            const auto gene = vector<State>(states.begin() + gene_start,
                                            states.begin() + gene_end);

            if (gene.front() != State::START_CODON_1) bad_start++;
            if (gene.back()  != State::STOP_CODON_3)  bad_stop++;

            int expected_frame = 1;
            bool frame_ok = true;
            for (State s : gene) {
                if (s == State::START_CODON_3) { expected_frame = 1; continue; }
                if (s == State::ACCEPTOR_1 || s == State::ACCEPTOR_2 || s == State::ACCEPTOR_3)
                    continue;
                if (s == State::EXON_FRAME_1 || s == State::EXON_FRAME_2 || s == State::EXON_FRAME_3) {
                    int actual = (s == State::EXON_FRAME_1) ? 1 :
                                 (s == State::EXON_FRAME_2) ? 2 : 3;
                    if (actual != expected_frame) { frame_ok = false; break; }
                    expected_frame = (expected_frame % 3) + 1;
                }
                if (s == State::STOP_CODON_1 && expected_frame != 1) { frame_ok = false; break; }
            }
            if (!frame_ok) bad_frame_seq++;

            int last_donor_frame = -1;
            bool da_ok = true;
            for (State s : gene) {
                if (s == State::DONOR_1)    last_donor_frame = 1;
                if (s == State::DONOR_2)    last_donor_frame = 2;
                if (s == State::DONOR_3)    last_donor_frame = 3;
                if (s == State::ACCEPTOR_1 && last_donor_frame != 1) { da_ok = false; break; }
                if (s == State::ACCEPTOR_2 && last_donor_frame != 2) { da_ok = false; break; }
                if (s == State::ACCEPTOR_3 && last_donor_frame != 3) { da_ok = false; break; }
            }
            if (!da_ok) donor_acceptor_mismatch++;
        }

        cout << "  total gene loci scanned: " << total_genes << "\n";
        CHECK("all genes open with START_CODON_1",       bad_start == 0);
        CHECK("all genes close with STOP_CODON_3",       bad_stop  == 0);
        CHECK("all genes have valid exon frame cycling", bad_frame_seq == 0);
        CHECK("all donor/acceptor frames match",         donor_acceptor_mismatch == 0);

        if (bad_start)               cout << "  " << bad_start << " gene(s) missing proper start\n";
        if (bad_stop)                cout << "  " << bad_stop  << " gene(s) missing proper stop\n";
        if (bad_frame_seq)           cout << "  " << bad_frame_seq << " gene(s) with broken frame cycling\n";
        if (donor_acceptor_mismatch) cout << "  " << donor_acceptor_mismatch << " gene(s) with mismatched donor/acceptor frames\n";
    }

    static void test_triplet_balance(const vector<State>& states) {
        cout << "\n[TEST 6] Triplet balance across all gene loci\n";

        unordered_map<int, size_t> start_cnt, stop_cnt, donor_cnt, acceptor_cnt, intron_cnt;
        for (State s : states) {
            switch (s) {
                case State::START_CODON_1:  start_cnt[1]++;     break;
                case State::START_CODON_2:  start_cnt[2]++;     break;
                case State::START_CODON_3:  start_cnt[3]++;     break;
                case State::STOP_CODON_1:   stop_cnt[1]++;      break;
                case State::STOP_CODON_2:   stop_cnt[2]++;      break;
                case State::STOP_CODON_3:   stop_cnt[3]++;      break;
                case State::DONOR_1:        donor_cnt[1]++;     break;
                case State::DONOR_2:        donor_cnt[2]++;     break;
                case State::DONOR_3:        donor_cnt[3]++;     break;
                case State::ACCEPTOR_1:     acceptor_cnt[1]++;  break;
                case State::ACCEPTOR_2:     acceptor_cnt[2]++;  break;
                case State::ACCEPTOR_3:     acceptor_cnt[3]++;  break;
                case State::INTRON_1:       intron_cnt[1]++;    break;
                case State::INTRON_2:       intron_cnt[2]++;    break;
                case State::INTRON_3:       intron_cnt[3]++;    break;
                default: break;
            }
        }

        CHECK("START_CODON _1/_2/_3 counts equal",
              start_cnt[1] == start_cnt[2] && start_cnt[2] == start_cnt[3]);
        CHECK("STOP_CODON _1/_2/_3 counts equal",
              stop_cnt[1] == stop_cnt[2] && stop_cnt[2] == stop_cnt[3]);
        CHECK("DONOR and ACCEPTOR counts equal per frame",
              donor_cnt[1] == acceptor_cnt[1] &&
              donor_cnt[2] == acceptor_cnt[2] &&
              donor_cnt[3] == acceptor_cnt[3]);

        cout << "  starts: " << start_cnt[1] << "/" << start_cnt[2] << "/" << start_cnt[3] << "\n";
        cout << "  stops:  " << stop_cnt[1]  << "/" << stop_cnt[2]  << "/" << stop_cnt[3]  << "\n";
        cout << "  donors:    " << donor_cnt[1] << "/" << donor_cnt[2] << "/" << donor_cnt[3] << "\n";
        cout << "  acceptors: " << acceptor_cnt[1] << "/" << acceptor_cnt[2] << "/" << acceptor_cnt[3] << "\n";
    }

    static void run_GFF_Parser_tests(const string& fna, const string& gff) {
        cout << "Loading genome: " << fna << "\n";
        vector<Nucleotide> nuc_seq = FNA_Parser::parse_sequence(fna);

        cout << "Parsing regions...\n";
        string gff_str = gff, fna_str = fna;
        vector<int> regions = GFF_Parser::parse_regions(gff_str, fna_str);

        cout << "Parsing states...\n";
        vector<State> states = GFF_Parser::parse_states(regions);

        test_sequence_length(fna, nuc_seq);
        test_regions(regions, nuc_seq);
        test_states_basic(states, nuc_seq);
        test_all_transitions_legal(states);
        test_gene_structure(states);
        test_triplet_balance(states);

        cout << "\nDone.\n";
    }

}
