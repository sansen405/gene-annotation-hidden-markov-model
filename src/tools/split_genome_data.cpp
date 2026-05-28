#include "../genome_profiles/Genome_Profile.hpp"
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

using namespace std;
using namespace gene_hmm;

string chrom_name_from_fasta_header(const string& line) {
    size_t end = 1;
    while (end < line.size() && !isspace(static_cast<unsigned char>(line[end]))) {
        end++;
    }
    return line.substr(1, end - 1);
}

string chrom_name_from_gff_line(const string& line) {
    size_t tab = line.find('\t');
    if (tab == string::npos) {
        return "";
    }
    return line.substr(0, tab);
}

bool is_gff_feature_line(const string& line) {
    if (line.empty() || line[0] == '#') {
        return false;
    }
    size_t tabs = 0;
    for (char c : line) {
        if (c == '\t') {
            tabs++;
        }
    }
    return tabs >= 7;
}

string sequence_region_chrom(const string& line) {
    const string prefix = "##sequence-region ";
    if (line.rfind(prefix, 0) != 0) {
        return "";
    }
    size_t start = prefix.size();
    size_t end = start;
    while (end < line.size() && line[end] != ' ') {
        end++;
    }
    return line.substr(start, end - start);
}

void write_fasta_split(
    const string& source_path,
    const string& output_path,
    const set<string>& include_chromosomes)
{
    ifstream input(source_path);
    if (!input.is_open()) {
        throw runtime_error("Cannot open FASTA: " + source_path);
    }

    ofstream output(output_path);
    if (!output.is_open()) {
        throw runtime_error("Cannot write FASTA: " + output_path);
    }

    string line;
    bool writing = false;
    while (getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            string chrom = chrom_name_from_fasta_header(line);
            writing = include_chromosomes.count(chrom) > 0;
        }
        if (writing) {
            output << line << '\n';
        }
    }
}

void write_gff_split(
    const string& source_path,
    const string& output_path,
    const set<string>& include_chromosomes)
{
    ifstream input(source_path);
    if (!input.is_open()) {
        throw runtime_error("Cannot open GFF: " + source_path);
    }

    ofstream output(output_path);
    if (!output.is_open()) {
        throw runtime_error("Cannot write GFF: " + output_path);
    }

    string line;
    while (getline(input, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            string region_chrom = sequence_region_chrom(line);
            if (!region_chrom.empty()) {
                if (include_chromosomes.count(region_chrom) > 0) {
                    output << line << '\n';
                }
                continue;
            }
            output << line << '\n';
            continue;
        }

        if (is_gff_feature_line(line)) {
            string chrom = chrom_name_from_gff_line(line);
            if (include_chromosomes.count(chrom) > 0) {
                output << line << '\n';
            }
            continue;
        }

        output << line << '\n';
    }
}

set<string> train_chromosomes_from_dataset(const Species_Dataset& dataset) {
    ifstream fasta(dataset.source_fasta_path);
    if (!fasta.is_open()) {
        throw runtime_error("Cannot open source FASTA: " + dataset.source_fasta_path);
    }

    set<string> test_set(dataset.test_chromosomes.begin(), dataset.test_chromosomes.end());
    set<string> excluded_set(dataset.excluded_chromosomes.begin(), dataset.excluded_chromosomes.end());
    set<string> train_set;

    string line;
    while (getline(fasta, line)) {
        if (line.empty() || line[0] != '>') {
            continue;
        }
        string chrom = chrom_name_from_fasta_header(line);
        if (test_set.count(chrom) || excluded_set.count(chrom)) {
            continue;
        }
        train_set.insert(chrom);
    }

    if (train_set.empty()) {
        throw runtime_error("No training chromosomes found in source FASTA.");
    }
    return train_set;
}

set<string> test_chromosomes_from_dataset(const Species_Dataset& dataset) {
    return set<string>(dataset.test_chromosomes.begin(), dataset.test_chromosomes.end());
}

} // namespace

int main(int argc, char** argv) {
    try {
        string profile_path = argc > 1 ? argv[1] : "src/genome_profiles/fission_yeasts.json";
        Genome_Profile profile = Genome_Profile::load(profile_path);

        if (profile.species.empty()) {
            throw runtime_error("Profile must define dataset.species.");
        }

        cout << "Splitting " << profile.name << " genome data\n";
        for (const auto& dataset : profile.species) {
            if (dataset.source_fasta_path.empty() || dataset.source_gff_path.empty()) {
                throw runtime_error("Species dataset must define source_fasta and source_gff: " + dataset.name);
            }
            if (dataset.test_chromosomes.empty()) {
                throw runtime_error("Species dataset must define test_chromosomes: " + dataset.name);
            }

            set<string> train_chromosomes = train_chromosomes_from_dataset(dataset);
            set<string> test_chromosomes = test_chromosomes_from_dataset(dataset);

            cout << "  " << dataset.name << "\n";
            cout << "    source FASTA: " << dataset.source_fasta_path << "\n";
            cout << "    source GFF:   " << dataset.source_gff_path << "\n";
            cout << "    train chromosomes: " << train_chromosomes.size() << "\n";
            cout << "    test chromosomes:  " << test_chromosomes.size() << "\n";
            cout << "    excluded:          " << dataset.excluded_chromosomes.size() << "\n";

            write_fasta_split(dataset.source_fasta_path, dataset.train_fasta_path, train_chromosomes);
            write_gff_split(dataset.source_gff_path, dataset.train_gff_path, train_chromosomes);
            write_fasta_split(dataset.source_fasta_path, dataset.test_fasta_path, test_chromosomes);
            write_gff_split(dataset.source_gff_path, dataset.test_gff_path, test_chromosomes);

            cout << "    wrote train FASTA: " << dataset.train_fasta_path << "\n";
            cout << "    wrote train GFF:   " << dataset.train_gff_path << "\n";
            cout << "    wrote test FASTA:  " << dataset.test_fasta_path << "\n";
            cout << "    wrote test GFF:    " << dataset.test_gff_path << "\n";
        }
        return 0;
    } catch (const exception& error) {
        cerr << error.what() << "\n";
        return 1;
    }
}
