#pragma once

#include <cstdint>
#include <vector>
#include <limits>
#include <unordered_map>
#include <string>

namespace gene_hmm {

    using namespace std;

    using LogProb = double;
    const LogProb LOG_INF = numeric_limits<LogProb>::infinity();
    const LogProb LOG_ZERO = -numeric_limits<LogProb>::infinity();

    enum class Nucleotide: uint8_t{
        A = 1,
        C = 2,
        G = 3,
        T = 4,
    };

    enum class State: uint8_t{
        START = 0, 
        INTERGENIC = 1,
        START_CODON_1 = 2, START_CODON_2 = 3, START_CODON_3 = 4,
        EXON_FRAME_1 = 5, EXON_FRAME_2 = 6, EXON_FRAME_3 = 7,
        DONOR_1 = 8, DONOR_2 = 9, DONOR_3 = 10,
        INTRON_1 = 11, INTRON_2 = 12, INTRON_3 = 13,
        ACCEPTOR_1 = 14, ACCEPTOR_2 = 15, ACCEPTOR_3 = 16,
        STOP_CODON_1 = 17, STOP_CODON_2 = 18, STOP_CODON_3 = 19,
        END = 20,
    };

    const size_t NUM_STATES = 21;
    const size_t NUM_NUCLEOTIDES = 4;
    const size_t MEMORY_WINDOW = 5;

    const unordered_map<State, vector<State>> Transitions = {
        {State::START, {State::INTERGENIC, State::START_CODON_1}},
        {State::INTERGENIC, {State::INTERGENIC, State::START_CODON_1, State::END}},
        {State::START_CODON_1, {State::START_CODON_2}},
        {State::START_CODON_2, {State::START_CODON_3}},
        {State::START_CODON_3, {State::EXON_FRAME_1, State::DONOR_1}},
        {State::EXON_FRAME_1, {State::EXON_FRAME_2, State::DONOR_2}},
        {State::EXON_FRAME_2, {State::EXON_FRAME_3, State::DONOR_3}},
        {State::EXON_FRAME_3, {State::EXON_FRAME_1, State::DONOR_1, State::STOP_CODON_1}},
        {State::DONOR_1, {State::INTRON_1}},
        {State::DONOR_2, {State::INTRON_2}},
        {State::DONOR_3, {State::INTRON_3}},
        {State::INTRON_1, {State::INTRON_1, State::ACCEPTOR_1}},
        {State::INTRON_2, {State::INTRON_2, State::ACCEPTOR_2}},
        {State::INTRON_3, {State::INTRON_3, State::ACCEPTOR_3}},
        {State::ACCEPTOR_1, {State::EXON_FRAME_1}},
        {State::ACCEPTOR_2, {State::EXON_FRAME_2}},
        {State::ACCEPTOR_3, {State::EXON_FRAME_3}},
        {State::STOP_CODON_1, {State::STOP_CODON_2}},
        {State::STOP_CODON_2, {State::STOP_CODON_3}},
        {State::STOP_CODON_3, {State::INTERGENIC, State::END}},
        {State::END, {}}
    };
}