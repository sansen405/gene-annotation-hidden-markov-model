#include "src/decoding/Viterbi.hpp"
#include "src/model/Model_IO.hpp"
#include "src/parsers/FNA_Parser.hpp"
#include <nlohmann/json.hpp>
#include <iostream>
#include <vector>

using namespace gene_hmm;
using namespace std;
using json = nlohmann::json;

namespace {

bool is_coding(State s) {
    return s == State::START_CODON_1 || s == State::START_CODON_2 || s == State::START_CODON_3 ||
           s == State::EXON_FRAME_1  || s == State::EXON_FRAME_2  || s == State::EXON_FRAME_3  ||
           s == State::STOP_CODON_1  || s == State::STOP_CODON_2  || s == State::STOP_CODON_3;
}

bool is_intron(State s) {
    return s == State::DONOR_1    || s == State::DONOR_2    || s == State::DONOR_3    ||
           s == State::INTRON_1   || s == State::INTRON_2   || s == State::INTRON_3   ||
           s == State::ACCEPTOR_1 || s == State::ACCEPTOR_2 || s == State::ACCEPTOR_3;
}

bool is_gene(State s) { return is_coding(s) || is_intron(s); }

struct Gene {
    size_t start = 0, end = 0;
    vector<pair<size_t,size_t>> exons;
    vector<pair<size_t,size_t>> introns;
};

vector<Gene> extract_genes(const vector<State>& states) {
    vector<Gene> genes;
    size_t i = 0;
    while (i < states.size()) {
        if (!is_gene(states[i])) { i++; continue; }
        Gene g;
        g.start = i;
        while (i < states.size() && is_gene(states[i])) i++;
        g.end = i;
        // Partition into consecutive exon/intron segments
        size_t j = g.start;
        while (j < g.end) {
            if (is_coding(states[j])) {
                size_t s = j;
                while (j < g.end && is_coding(states[j])) j++;
                g.exons.push_back({s, j});
            } else {
                size_t s = j;
                while (j < g.end && is_intron(states[j])) j++;
                g.introns.push_back({s, j});
            }
        }
        genes.push_back(g);
    }
    return genes;
}

} // namespace

int main(int argc, char** argv) {
    string model_path, fasta_path;
    for (int i = 1; i + 1 < argc; i += 2) {
        string flag = argv[i];
        if (flag == "--model") model_path = argv[i + 1];
        else if (flag == "--fasta") fasta_path = argv[i + 1];
    }

    if (model_path.empty() || fasta_path.empty()) {
        cerr << "Usage: " << argv[0] << " --model <model.json> --fasta <genome.fna>\n";
        return 1;
    }

    Trained_Model model;
    try {
        model = Model_IO::load(model_path);
    } catch (const exception& e) {
        cerr << "Failed to load model: " << e.what() << "\n";
        return 1;
    }

    vector<Nucleotide> nucleotides;
    vector<Chromosome_Range> ranges;
    try {
        nucleotides = FNA_Parser::parse_sequence(fasta_path);
        ranges      = FNA_Parser::get_chromosome_ranges(fasta_path);
    } catch (const exception& e) {
        cerr << "Failed to parse FASTA: " << e.what() << "\n";
        return 1;
    }

    json output;
    output["genome"]      = model.genome_name;
    output["chromosomes"] = json::array();
    size_t total_genes    = 0;

    for (const auto& range : ranges) {
        vector<Nucleotide> chr_nucs(
            nucleotides.begin() + range.start,
            nucleotides.begin() + range.end);

        vector<State> predicted = Viterbi::decode(
            chr_nucs,
            model.transition_log_probs,
            model.emission_model,
            0,
            model.max_intron_body_length,
            model.gene_start_penalty);

        auto genes = extract_genes(predicted);
        total_genes += genes.size();

        json chr_json;
        chr_json["name"]   = range.name;
        chr_json["length"] = range.end - range.start;
        chr_json["genes"]  = json::array();

        for (size_t g = 0; g < genes.size(); ++g) {
            const auto& gene = genes[g];
            json gj;
            gj["id"]    = "gene_" + to_string(g + 1);
            gj["start"] = gene.start;
            gj["end"]   = gene.end;

            gj["exons"] = json::array();
            for (const auto& [s, e] : gene.exons)
                gj["exons"].push_back({{"start", s}, {"end", e}});

            gj["introns"] = json::array();
            for (const auto& [s, e] : gene.introns)
                gj["introns"].push_back({{"start", s}, {"end", e}});

            chr_json["genes"].push_back(gj);
        }
        output["chromosomes"].push_back(chr_json);
    }

    output["summary"] = {
        {"total_genes",       total_genes},
        {"total_chromosomes", ranges.size()}
    };

    cout << output.dump(2) << "\n";
    return 0;
}
