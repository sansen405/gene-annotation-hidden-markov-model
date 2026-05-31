#include "GFF_Parser.hpp"
#include "FNA_Parser.hpp"
#include "../genome_profiles/Genome_Profile.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>

namespace gene_hmm {
    using namespace std;

    struct Exon_Fragment {
        int start;
        int end;
        bool operator<(const Exon_Fragment& other) const { return start < other.start; }
    };

    struct CDS_Group {
        vector<Exon_Fragment> fragments;
        bool supported = true;
    };

    static bool group_passes_filters(const vector<Exon_Fragment>& fragments) {
        if (fragments.empty()) {
            return false;
        }

        size_t total_cds_bp = 0;
        for (const auto& fragment : fragments) {
            total_cds_bp += static_cast<size_t>(fragment.end - fragment.start + 1);
        }

        size_t first_cds_bp = static_cast<size_t>(fragments.front().end - fragments.front().start + 1);
        size_t last_cds_bp = static_cast<size_t>(fragments.back().end - fragments.back().start + 1);
        if (first_cds_bp < profile.min_first_cds_bp || last_cds_bp < profile.min_last_cds_bp) {
            return false;
        }

        if (profile.require_3n_cds && total_cds_bp % 3 != 0) {
            return false;
        }

        for (size_t i = 0; i + 1 < fragments.size(); ++i) {
            int intron_bp = fragments[i + 1].start - fragments[i].end - 1;
            if (intron_bp > 0 && static_cast<size_t>(intron_bp) < profile.min_intron_bp) {
                return false;
            }
        }

        return true;
    }

    static void paint_region(
        vector<int>& regions,
        int start,
        int end,
        int region)
    {
        for (int index = start; index <= end; index++) {
            if (index >= 0 && static_cast<size_t>(index) < regions.size()) {
                regions[index] = region;
            }
        }
    }

    vector<int> GFF_Parser::parse_regions(string& gff_path, string& fna_path) {
        return parse_regions(gff_path, fna_path, '+');
    }

    vector<int> GFF_Parser::parse_regions(string& gff_path, string& fna_path, char strand) {
        size_t sequence_length = FNA_Parser::get_sequence_length(fna_path);
        auto chrom_offsets = FNA_Parser::get_chromosome_offsets(fna_path);

        //pivot K maps a forward global position p to revcomp position K - p within the chromosome
        unordered_map<string, int> revcomp_pivot;
        for (const auto& range : FNA_Parser::get_chromosome_ranges(fna_path)) {
            revcomp_pivot[range.name] =
                static_cast<int>(range.start) + static_cast<int>(range.end) - 1;
        }

        ifstream file(gff_path);
        if (!file.is_open()) {
            throw runtime_error("Unable to Open GFF File");
        }

        vector<int> regions_sequence(sequence_length, GFF_Parser::INTERGENIC_REGION);
        string curr_line = "";

        vector<string> tokens;
        string t;

        map<string, CDS_Group> gene_builder;

        while (getline(file, curr_line)) {
            if (curr_line.empty() || curr_line[0] == '#') {
                continue;
            }

            tokens.clear();
            stringstream line_stream(curr_line);

            while (getline(line_stream, t, '\t')) {
                tokens.push_back(t);
            }

            if (tokens.size() != 9 || tokens[2] != "CDS") {
                continue;
            }

            auto off_it = chrom_offsets.find(tokens[0]);
            if (off_it == chrom_offsets.end()) continue;
            int chrom_off = static_cast<int>(off_it->second);

            int CDS_curr_start = stoi(tokens[3]) - 1 + chrom_off;
            int CDS_curr_end = stoi(tokens[4]) - 1 + chrom_off;

            //in minus mode every fragment is mirrored within its chromosome and start/end swap
            if (strand == '-') {
                int pivot = revcomp_pivot.at(tokens[0]);
                int mirrored_start = pivot - CDS_curr_end;
                int mirrored_end = pivot - CDS_curr_start;
                CDS_curr_start = mirrored_start;
                CDS_curr_end = mirrored_end;
            }

            string attributes = tokens[8];
            string CDS_curr_parent_id = "";
            size_t parent_pos = attributes.find("Parent=");
            if (parent_pos != string::npos) {
                size_t start_val = parent_pos + 7;
                size_t end_val = attributes.find(";", start_val);
                CDS_curr_parent_id = attributes.substr(start_val, end_val - start_val);
            }

            if (!CDS_curr_parent_id.empty()) {
                string key = tokens[0] + "|" + CDS_curr_parent_id;
                gene_builder[key].fragments.push_back({CDS_curr_start, CDS_curr_end});
                //opposite-strand transcripts are marked unsupported so they are painted as ignored
                string required_strand(1, strand);
                if (tokens[6] != required_strand) {
                    gene_builder[key].supported = false;
                }
            }
        }

        for (auto& [parent_id, group] : gene_builder) {
            sort(group.fragments.begin(), group.fragments.end());
            group.supported = group.supported && group_passes_filters(group.fragments);
        }

        for (auto& [parent_id, group] : gene_builder) {
            if (group.supported) {
                continue;
            }

            for (size_t i = 0; i < group.fragments.size(); ++i) {
                paint_region(
                    regions_sequence,
                    group.fragments[i].start,
                    group.fragments[i].end,
                    GFF_Parser::IGNORED_REGION);

                if (i + 1 < group.fragments.size()) {
                    paint_region(
                        regions_sequence,
                        group.fragments[i].end + 1,
                        group.fragments[i + 1].start - 1,
                        GFF_Parser::IGNORED_REGION);
                }
            }
        }

        for (auto& [parent_id, group] : gene_builder) {
            if (!group.supported) {
                continue;
            }

            for (size_t i = 0; i < group.fragments.size(); ++i) {
                paint_region(
                    regions_sequence,
                    group.fragments[i].start,
                    group.fragments[i].end,
                    GFF_Parser::CDS_REGION);

                if (i + 1 < group.fragments.size()) {
                    paint_region(
                        regions_sequence,
                        group.fragments[i].end + 1,
                        group.fragments[i + 1].start - 1,
                        GFF_Parser::INTRON_REGION);
                }
            }
        }
        return regions_sequence;
    }

    vector<State> GFF_Parser::parse_states(vector<int> region_sequence) {
        vector<State> state_sequence(region_sequence.size(), State::INTERGENIC);
        int frame_counter = 0;
        int active_intron = 0;
        size_t index = 0;
        while(index < region_sequence.size()){
            if(region_sequence[index] == GFF_Parser::INTERGENIC_REGION ||
               region_sequence[index] == GFF_Parser::IGNORED_REGION){
                state_sequence[index] = State::INTERGENIC;
                frame_counter = 0;
                index++;
                continue;
            }

            if(region_sequence[index] == GFF_Parser::CDS_REGION){
                if(index == 0 || region_sequence[index-1] == GFF_Parser::INTERGENIC_REGION ||
                   region_sequence[index-1] == GFF_Parser::IGNORED_REGION){
                    if (index + 2 < region_sequence.size()
                        && region_sequence[index+1] == GFF_Parser::CDS_REGION
                        && region_sequence[index+2] == GFF_Parser::CDS_REGION) {
                        state_sequence[index] = State::START_CODON_1;
                        state_sequence[index+1] = State::START_CODON_2;
                        state_sequence[index+2] = State::START_CODON_3;
                        index += 3;
                        frame_counter = 0;
                        continue;
                    }
                }

                if(index + 2 < region_sequence.size()
                   && region_sequence[index+1] == GFF_Parser::CDS_REGION
                   && region_sequence[index+2] == GFF_Parser::CDS_REGION
                   && (index + 3 >= region_sequence.size() ||
                       region_sequence[index+3] == GFF_Parser::INTERGENIC_REGION ||
                       region_sequence[index+3] == GFF_Parser::IGNORED_REGION)){
                    state_sequence[index] = State::STOP_CODON_1;
                    state_sequence[index+1] = State::STOP_CODON_2;
                    state_sequence[index+2] = State::STOP_CODON_3;

                    index += 3;
                    continue;
                }

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

            if (region_sequence[index] == GFF_Parser::INTRON_REGION){
                bool is_donor = (index == 0 || region_sequence[index-1] == GFF_Parser::CDS_REGION);
                bool is_acceptor = (index + 1 < region_sequence.size() && region_sequence[index+1] == GFF_Parser::CDS_REGION);

                if(is_donor){
                    if(frame_counter == 0){
                        state_sequence[index] = State::DONOR_1;
                        active_intron = frame_counter + 1;
                    }
                    if(frame_counter == 1){
                        state_sequence[index] = State::DONOR_2;
                        active_intron = frame_counter + 1;
                    }
                    if(frame_counter == 2){
                        state_sequence[index] = State::DONOR_3;
                        active_intron = frame_counter + 1;
                    }
                    index++;
                    continue;
                }

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

        //reset any gene containing an illegal transition back to intergenic
        int most_recent_intergenic = -1;
        int num_illegal_transitions = 0;
        bool is_illegal_transition = false;


        for(size_t state_index = 0; state_index < state_sequence.size()-1; state_index++){
            if(state_sequence[state_index] == State::INTERGENIC){
                if((int)state_index-most_recent_intergenic > 1 && is_illegal_transition){
                    for(size_t i = most_recent_intergenic+1; i < state_index; i++){
                        state_sequence[i] = State::INTERGENIC;
                    }
                    is_illegal_transition = false;
                    num_illegal_transitions++;
                }
                most_recent_intergenic = state_index;
            }
            if(find(Transitions.at(state_sequence[state_index]).begin(), 
                Transitions.at(state_sequence[state_index]).end(), 
                state_sequence[state_index+1]) == Transitions.at(state_sequence[state_index]).end()){
                is_illegal_transition = true;
            }
        }
        if(is_illegal_transition)
            for(size_t i = most_recent_intergenic+1; i < state_sequence.size(); i++)
                state_sequence[i] = State::INTERGENIC;

        return state_sequence;
    }
}