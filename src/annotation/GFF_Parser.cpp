#include "GFF_Parser.hpp"
#include "Sequence_Parser.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>

namespace gene_hmm {
    using namespace std;

    struct ExonFragment {
        int start;
        int end;
        bool operator<(const ExonFragment& other) const { return start < other.start; }
    };

    vector<int> GFF_Parser::parse_regions(string& gff_path, string& fna_path) {
        // regions: {0 -> Intergenic, 1 -> CDS, 2 -> Intron}

        size_t sequence_length = Sequence_Parser::get_sequence_length(fna_path);
        auto chrom_offsets = Sequence_Parser::get_chromosome_offsets(fna_path);

        ifstream file(gff_path);
        if (!file.is_open()) {
            throw runtime_error("Unable to Open GFF File");
        }

        vector<int> regions_sequence(sequence_length, 0);
        string curr_line = "";

        vector<string> tokens;
        string t;

        map<string, vector<ExonFragment>> gene_builder;

        //group fragments by parent_id (namespaced by chromosome)
        while (getline(file, curr_line)) {
            if (curr_line.empty() || curr_line[0] == '#') {
                continue;
            }

            tokens.clear();
            stringstream line_stream(curr_line);

            while (getline(line_stream, t, '\t')) {
                tokens.push_back(t);
            }

            if (tokens.size() != 9 || tokens[2] != "CDS" || tokens[6] != "+") {
                continue;
            }

            // Translate per-chromosome GFF coords into the global concatenated vector.
            auto off_it = chrom_offsets.find(tokens[0]);
            if (off_it == chrom_offsets.end()) continue;
            int chrom_off = static_cast<int>(off_it->second);

            int CDS_curr_start = stoi(tokens[3]) - 1 + chrom_off;
            int CDS_curr_end = stoi(tokens[4]) - 1 + chrom_off;

            string attributes = tokens[8];
            string CDS_curr_parent_id = "";
            size_t parent_pos = attributes.find("Parent=");
            if (parent_pos != string::npos) {
                size_t start_val = parent_pos + 7;
                size_t end_val = attributes.find(";", start_val);
                CDS_curr_parent_id = attributes.substr(start_val, end_val - start_val);
            }

            if (!CDS_curr_parent_id.empty()) {
                gene_builder[tokens[0] + "|" + CDS_curr_parent_id]
                    .push_back({CDS_curr_start, CDS_curr_end});
            }
        }

        //map back to regions
        for (auto& [parent_id, fragments] : gene_builder) {
            sort(fragments.begin(), fragments.end());

            for (size_t i = 0; i < fragments.size(); ++i) {
                for (int index = fragments[i].start; index <= fragments[i].end; ++index) {
                    if (index >= 0 && static_cast<size_t>(index) < sequence_length) {
                        regions_sequence[index] = 1;
                    }
                }

                if (i + 1 < fragments.size()) {
                    int intron_start = fragments[i].end + 1;
                    int intron_end = fragments[i+1].start - 1;

                    for (int index = intron_start; index <= intron_end; ++index) {
                        if (index >= 0 && static_cast<size_t>(index) < sequence_length) {
                            regions_sequence[index] = 2;
                        }
                    }
                }
            }
        }
        return regions_sequence;
    }

    vector<State> GFF_Parser::parse_states(vector<int> region_sequence) {
        vector<State> state_sequence(region_sequence.size(), State::INTERGENIC);
        int frame_counter = 0; //0 -> exon_1, etc.
        int active_intron = 0;
        size_t index = 0;
        while(index < region_sequence.size()){
            //CASE 0: INTERGENIC STATES
            if(region_sequence[index] == 0){
                state_sequence[index] = State::INTERGENIC;
                frame_counter = 0;
                index++;
                continue;
            }

            //CASE 1: CDS STATES
            if(region_sequence[index] == 1){
                //CASE 1A: START CODON
                if(index == 0 || region_sequence[index-1] == 0){ //INTERGENIC -> START_CODONS
                    if (index + 2 < region_sequence.size()
                        && region_sequence[index+1] == 1
                        && region_sequence[index+2] == 1) {
                        state_sequence[index] = State::START_CODON_1;
                        state_sequence[index+1] = State::START_CODON_2;
                        state_sequence[index+2] = State::START_CODON_3;
                        index += 3;
                        frame_counter = 0;
                        continue;
                    }
                    // First CDS fragment is split by an intron (rare); fall through
                    // and treat this base as a normal exon-frame_1 base.
                }

                //CASE 1B: END_CODON
                if(index + 2 < region_sequence.size()
                   && region_sequence[index+1] == 1
                   && region_sequence[index+2] == 1
                   && (index + 3 >= region_sequence.size() || region_sequence[index+3] == 0)){ //CDS -> intergenic (or end)
                    state_sequence[index] = State::STOP_CODON_1;
                    state_sequence[index+1] = State::STOP_CODON_2;
                    state_sequence[index+2] = State::STOP_CODON_3;

                    index += 3;
                    continue;
                }

                //CASE 1: STANDARD CDS EXON READING FRAMES
                if(frame_counter == 0){
                    state_sequence[index] = State::EXON_FRAME_1;
                }
                if(frame_counter == 1){
                    state_sequence[index] = State::EXON_FRAME_2;
                }
                if(frame_counter == 2){
                    state_sequence[index] = State::EXON_FRAME_3;
                }
                
                frame_counter = (frame_counter+1)%3;
                index++;
                continue;
            }

            //CASE 2: INTRON STATES
            if (region_sequence[index] == 2){
                bool is_donor = (index == 0 || region_sequence[index-1] == 1);
                bool is_acceptor = (index + 1 < region_sequence.size() && region_sequence[index+1] == 1);

                //Case 2A: DONOR_1,2,3 states
                if(is_donor){
                    if(frame_counter == 0){//Exon 3 Spun off
                        state_sequence[index] = State::DONOR_1;
                        active_intron = frame_counter +1;
                    }
                    if(frame_counter == 1){ // Exon1 Spun oFF
                        state_sequence[index] = State::DONOR_2;
                        active_intron = frame_counter +1;
                    }
                    if(frame_counter == 2){// Expon2 Spun off
                        state_sequence[index] = State::DONOR_3;
                        active_intron = frame_counter + 1;
                    }
                    index++;
                    continue;
                }
                //Case 2B: Acceptor 1,2,3
                if(is_acceptor){
                    if(active_intron == 1){
                        state_sequence[index] = State::ACCEPTOR_1;
                    }
                    else if(active_intron == 2){
                        state_sequence[index] = State::ACCEPTOR_2;
                    }
                    else if(active_intron == 3){
                        state_sequence[index] = State::ACCEPTOR_3;
                    }
                    index++;
                    continue;
                }
                //CADE2C: Intron 
                if(active_intron == 1){
                    state_sequence[index] = State::INTRON_1;
                }
                else if(active_intron == 2){
                    state_sequence[index] = State::INTRON_2;
        
                }
                else if(active_intron == 3){
                    state_sequence[index] = State::INTRON_3;
                }

                index++;
            }

        }

        return state_sequence;
    }
}