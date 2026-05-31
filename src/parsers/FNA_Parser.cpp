#include "FNA_Parser.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cctype>

namespace gene_hmm {
    using namespace std;

    size_t FNA_Parser::get_sequence_length(const string& fasta_path) {
        ifstream file(fasta_path);
        if (!file.is_open()) {
            throw runtime_error("Unable to open FASTA file: " + fasta_path);
        }
        size_t total_length = 0;
        string line;
        while (getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') continue;
            for (char c : line) {
                if (c != '\n' && c != '\r') {
                    total_length++;
                }
            }
        }
        return total_length;
    }

    vector<Chromosome_Range> FNA_Parser::get_chromosome_ranges(const string& fasta_path) {
        ifstream file(fasta_path);
        if (!file.is_open()) {
            throw runtime_error("Unable to open FASTA file: " + fasta_path);
        }
        vector<Chromosome_Range> ranges;
        size_t cumulative = 0;
        string line;
        while (getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                //Close the previous range now that we know where it ended
                if (!ranges.empty()) ranges.back().end = cumulative;
                size_t end = 1;
                while (end < line.size() &&
                       !isspace(static_cast<unsigned char>(line[end]))) {
                    end++;
                }
                ranges.push_back({line.substr(1, end - 1), cumulative, 0});
                continue;
            }
            for (char c : line) {
                if (c != '\n' && c != '\r') cumulative++;
            }
        }
        if (!ranges.empty()) ranges.back().end = cumulative;
        return ranges;
    }

    unordered_map<string, size_t> FNA_Parser::get_chromosome_offsets(const string& fasta_path) {
        ifstream file(fasta_path);
        if (!file.is_open()) {
            throw runtime_error("Unable to open FASTA file: " + fasta_path);
        }
        unordered_map<string, size_t> offsets;
        size_t cumulative = 0;
        string line;
        while (getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                size_t end = 1;
                while (end < line.size() &&
                       !isspace(static_cast<unsigned char>(line[end]))) {
                    end++;
                }
                string chrom = line.substr(1, end - 1);
                offsets[chrom] = cumulative;
                continue;
            }
            for (char c : line) {
                if (c != '\n' && c != '\r') cumulative++;
            }
        }
        return offsets;
    }

    vector<Nucleotide> FNA_Parser::reverse_complement(const vector<Nucleotide>& nucleotides) {
        vector<Nucleotide> result(nucleotides.size());
        const size_t n = nucleotides.size();
        for (size_t i = 0; i < n; i++) {
            // 5 - value maps A<->T (1<->4) and C<->G (2<->3); reverse the order.
            uint8_t complement = static_cast<uint8_t>(5 - static_cast<uint8_t>(nucleotides[n - 1 - i]));
            result[i] = static_cast<Nucleotide>(complement);
        }
        return result;
    }

    vector<Nucleotide> FNA_Parser::parse_sequence(const string& fasta_path) {
        ifstream file(fasta_path);
        if (!file.is_open()) {
            throw runtime_error("Unable to Open Sequence File");
        }
        vector<Nucleotide> nuc_sequence;
        string line;
        while (getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') continue;
            for (char c : line) {
                c = toupper(static_cast<unsigned char>(c));
                switch (c) {
                    case 'A': nuc_sequence.push_back(Nucleotide::A); break;
                    case 'C': nuc_sequence.push_back(Nucleotide::C); break;
                    case 'G': nuc_sequence.push_back(Nucleotide::G); break;
                    case 'T': nuc_sequence.push_back(Nucleotide::T); break;
                    default:
                        // Preserve FASTA/GFF coordinate alignment for ambiguous bases.
                        nuc_sequence.push_back(Nucleotide::A);
                        break;
                }
            }
        }
        return nuc_sequence;
    }
}
