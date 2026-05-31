#include "../decoding/Forward_Backward.hpp"
#include "../decoding/Viterbi.hpp"
#include "../genome_profiles/Genome_Profile.hpp"
#include "../model/emission/Emission_Model.hpp"
#include "../model/transition/Transition_Model.hpp"
#include "../parsers/FNA_Parser.hpp"
#include "../parsers/GFF_Parser.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

using namespace gene_hmm;
using namespace std;

Log_Prob gene_start_penalty = 1.0;

struct Training_Data {
    vector<Nucleotide> nucleotides;
    vector<State> states;
    vector<int> regions;
    vector<Chromosome_Range> chromosomes;
};

bool is_coding(State state) {
    return state == State::START_CODON_1 ||
           state == State::START_CODON_2 ||
           state == State::START_CODON_3 ||
           state == State::EXON_FRAME_1 ||
           state == State::EXON_FRAME_2 ||
           state == State::EXON_FRAME_3 ||
           state == State::STOP_CODON_1 ||
           state == State::STOP_CODON_2 ||
           state == State::STOP_CODON_3;
}

bool is_intron(State state) {
    return state == State::DONOR_1 ||
           state == State::DONOR_2 ||
           state == State::DONOR_3 ||
           state == State::INTRON_1 ||
           state == State::INTRON_2 ||
           state == State::INTRON_3 ||
           state == State::ACCEPTOR_1 ||
           state == State::ACCEPTOR_2 ||
           state == State::ACCEPTOR_3;
}

bool is_intron_body(State state) {
    return state == State::INTRON_1 ||
           state == State::INTRON_2 ||
           state == State::INTRON_3;
}

bool is_gene(State state) {
    return is_coding(state) || is_intron(state);
}

bool is_usable_region(int region) {
    return region != GFF_Parser::IGNORED_REGION;
}

string state_family(State state) {
    if(state == State::START_CODON_1 ||
       state == State::START_CODON_2 ||
       state == State::START_CODON_3) {
        return "Start";
    }
    if(state == State::STOP_CODON_1 ||
       state == State::STOP_CODON_2 ||
       state == State::STOP_CODON_3) {
        return "Stop";
    }
    if(state == State::EXON_FRAME_1 ||
       state == State::EXON_FRAME_2 ||
       state == State::EXON_FRAME_3) {
        return "UTR";
    }
    if(is_intron(state)) {
        return "Intron";
    }
    return "Intergenic";
}

char nucleotide_char(Nucleotide nucleotide) {
    switch(nucleotide){
        case Nucleotide::A: return 'A';
        case Nucleotide::C: return 'C';
        case Nucleotide::G: return 'G';
        case Nucleotide::T: return 'T';
    }
    return 'N';
}

string sequence_preview(const vector<Nucleotide>& nucleotides, size_t start, size_t end) {
    string sequence;
    for(size_t i = start; i < end; i++){
        sequence.push_back(nucleotide_char(nucleotides[i]));
    }
    return sequence;
}

vector<Nucleotide> slice_nucleotides(
    const vector<Nucleotide>& nucleotides,
    size_t start,
    size_t end)
{
    return vector<Nucleotide>(nucleotides.begin() + start, nucleotides.begin() + end);
}

vector<Chromosome_Range> split_usable_ranges(
    const vector<Chromosome_Range>& ranges,
    const vector<int>& regions)
{
    vector<Chromosome_Range> usable_ranges;

    for(const auto& range : ranges){
        size_t pos = range.start;
        while(pos < range.end){
            while(pos < range.end && !is_usable_region(regions[pos])){
                pos++;
            }
            if(pos >= range.end){
                break;
            }

            size_t start = pos;
            while(pos < range.end && is_usable_region(regions[pos])){
                pos++;
            }
            usable_ranges.push_back({range.name, start, pos});
        }
    }

    return usable_ranges;
}

void append_training_dataset(Training_Data& combined, const Species_Dataset& dataset) {
    size_t offset = combined.nucleotides.size();

    vector<Nucleotide> nucleotides = FNA_Parser::parse_sequence(dataset.train_fasta_path);
    vector<Chromosome_Range> chromosomes = FNA_Parser::get_chromosome_ranges(dataset.train_fasta_path);
    string gff_path = dataset.train_gff_path;
    string fasta_path = dataset.train_fasta_path;
    vector<int> regions = GFF_Parser::parse_regions(gff_path, fasta_path);
    vector<State> states = GFF_Parser::parse_states(regions);

    if (states.size() != nucleotides.size()) {
        throw runtime_error("Training state length does not match nucleotide length for " + dataset.name + ".");
    }

    combined.nucleotides.insert(combined.nucleotides.end(), nucleotides.begin(), nucleotides.end());
    combined.states.insert(combined.states.end(), states.begin(), states.end());
    combined.regions.insert(combined.regions.end(), regions.begin(), regions.end());
    for (auto range : chromosomes) {
        range.name = dataset.name + ":" + range.name;
        range.start += offset;
        range.end += offset;
        combined.chromosomes.push_back(range);
    }
}

Training_Data load_training_data(const Genome_Profile& profile) {
    Training_Data combined;
    for (const auto& dataset : profile.species) {
        append_training_dataset(combined, dataset);
    }
    return combined;
}

Emission_Model train_emissions(
    const vector<State>& states,
    const vector<Nucleotide>& nucleotides,
    const vector<Chromosome_Range>& train_ranges)
{
    Emission_Model model;

    model.intergenic_lp = Emission_Model::compute_markov1_log_probs(
        Emission_Model::count_markov1_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::INTERGENIC}));

    model.intron_lp = Emission_Model::compute_markov1_log_probs(
        Emission_Model::count_markov1_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::INTRON_1, State::INTRON_2, State::INTRON_3}));

    model.exon_lp = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::EXON_FRAME_1, State::EXON_FRAME_2, State::EXON_FRAME_3}));
    model.exon_frame_lp[0] = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::EXON_FRAME_1}));
    model.exon_frame_lp[1] = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::EXON_FRAME_2}));
    model.exon_frame_lp[2] = Emission_Model::compute_markov5_log_probs(
        Emission_Model::count_markov5_emissions(
            states,
            nucleotides,
            train_ranges,
            {State::EXON_FRAME_3}));

    auto start_targets = vector<State>{State::START_CODON_1};
    model.start_codon_lp = Emission_Model::compute_pssm_log_odds(
        Emission_Model::count_pssm_emissions(
            states,
            nucleotides,
            train_ranges,
            start_targets,
            model.start_window_left,
            model.start_window_right),
        Emission_Model::count_pssm_background_emissions(
            states,
            nucleotides,
            train_ranges,
            start_targets,
            model.start_window_left,
            model.start_window_right,
            Splice_Signal::START_CODON));

    auto donor_targets = vector<State>{State::DONOR_1, State::DONOR_2, State::DONOR_3};
    model.donor_lp = Emission_Model::compute_pssm_log_odds(
        Emission_Model::count_pssm_emissions(
            states,
            nucleotides,
            train_ranges,
            donor_targets,
            model.donor_window_left,
            model.donor_window_right),
        Emission_Model::count_pssm_background_emissions(
            states,
            nucleotides,
            train_ranges,
            donor_targets,
            model.donor_window_left,
            model.donor_window_right,
            Splice_Signal::DONOR));

    auto acceptor_targets = vector<State>{State::ACCEPTOR_1, State::ACCEPTOR_2, State::ACCEPTOR_3};
    model.acceptor_lp = Emission_Model::compute_pssm_log_odds(
        Emission_Model::count_pssm_emissions(
            states,
            nucleotides,
            train_ranges,
            acceptor_targets,
            model.acceptor_window_left,
            model.acceptor_window_right),
        Emission_Model::count_pssm_background_emissions(
            states,
            nucleotides,
            train_ranges,
            acceptor_targets,
            model.acceptor_window_left,
            model.acceptor_window_right,
            Splice_Signal::ACCEPTOR));

    return model;
}

json intervals_for_family(
    const vector<State>& states,
    size_t offset,
    size_t start,
    size_t end,
    bool (*predicate)(State))
{
    json intervals = json::array();
    size_t pos = start;
    while(pos < end){
        if(!predicate(states[pos])){
            pos++;
            continue;
        }

        size_t interval_start = pos;
        while(pos < end && predicate(states[pos])){
            pos++;
        }
        intervals.push_back({
            {"start", interval_start - offset + 1},
            {"end", pos - offset}
        });
    }
    return intervals;
}

json confidence_for_range(
    const vector<State>& states,
    const vector<double>& confidence,
    size_t offset,
    size_t start,
    size_t end)
{
    json values = json::array();
    for(size_t pos = start; pos < end; pos++){
        values.push_back({
            {"position", pos - offset + 1},
            {"state", state_family(states[pos])},
            {"confidence", confidence[pos - offset]}
        });
    }
    return values;
}

// A decoded gene held in forward-genome coordinates so plus- and minus-strand
// predictions can be merged. score is the Viterbi path log-probability of the gene
// sub-path, used to resolve overlaps between strands.
struct Gene_Prediction {
    size_t fwd_start = 0; // 0-based inclusive, forward coordinates
    size_t fwd_end = 0;   // 0-based exclusive, forward coordinates
    Log_Prob score = 0.0;
    json payload;         // prediction JSON without the id
};

// Map exon/intron intervals from reverse-complement coordinates back to forward
// coordinates. Input intervals are 1-based inclusive in revcomp-local space; output
// is 1-based inclusive in forward-local space, sorted ascending by start.
json map_intervals_to_forward(const json& revcomp_intervals, size_t length) {
    json out = json::array();
    for(const auto& interval : revcomp_intervals){
        size_t rev_start = interval["start"].get<size_t>(); // 1-based inclusive
        size_t rev_end = interval["end"].get<size_t>();      // 1-based inclusive
        // revcomp 0-based [rev_start-1, rev_end) -> forward 0-based [L-rev_end, L-rev_start+1)
        size_t fwd_start = length - rev_end;       // 0-based inclusive
        size_t fwd_end = length - (rev_start - 1); // 0-based exclusive
        out.push_back({{"start", fwd_start + 1}, {"end", fwd_end}});
    }
    sort(out.begin(), out.end(), [](const json& a, const json& b){
        return a["start"].get<size_t>() < b["start"].get<size_t>();
    });
    return out;
}

// Per-base confidence for a minus-strand gene, reported at forward coordinates.
json confidence_for_range_minus(
    const vector<State>& states,
    const vector<double>& confidence,
    size_t start,
    size_t end,
    size_t length)
{
    json values = json::array();
    for(size_t pos = start; pos < end; pos++){
        size_t fwd_position = length - 1 - pos; // 0-based forward
        values.push_back({
            {"position", fwd_position + 1},
            {"state", state_family(states[pos])},
            {"confidence", confidence[pos]}
        });
    }
    sort(values.begin(), values.end(), [](const json& a, const json& b){
        return a["position"].get<size_t>() < b["position"].get<size_t>();
    });
    return values;
}

void collect_plus_genes(
    vector<Gene_Prediction>& genes,
    const vector<Nucleotide>& chromosome_nucleotides,
    const vector<State>& chromosome_states,
    const vector<double>& chromosome_confidence,
    const string& scaffold,
    const Transition_Model::Log_Prob_Matrix& transition_log_probs,
    const Emission_Model& emission_model)
{
    size_t pos = 0;
    while(pos < chromosome_states.size()){
        if(!is_gene(chromosome_states[pos])){
            pos++;
            continue;
        }

        size_t gene_start = pos;
        while(pos < chromosome_states.size() && is_gene(chromosome_states[pos])){
            pos++;
        }
        size_t gene_end = pos;

        json prediction;
        prediction["scaffold"] = scaffold;
        prediction["strand"] = "+";
        prediction["start"] = gene_start + 1;
        prediction["end"] = gene_end;
        prediction["exons"] = intervals_for_family(chromosome_states, 0, gene_start, gene_end, is_coding);
        prediction["introns"] = intervals_for_family(chromosome_states, 0, gene_start, gene_end, is_intron);
        prediction["intron_count"] = intervals_for_family(
            chromosome_states, 0, gene_start, gene_end, is_intron_body).size();
        prediction["sequence"] = sequence_preview(chromosome_nucleotides, gene_start, gene_end);
        prediction["confidence"] = confidence_for_range(
            chromosome_states, chromosome_confidence, 0, gene_start, gene_end);

        Log_Prob score = Viterbi::path_log_prob(
            chromosome_states, chromosome_nucleotides, transition_log_probs,
            emission_model, gene_start_penalty, gene_start, gene_end);

        genes.push_back({gene_start, gene_end, score, prediction});
    }
}

// chromosome_states/nucleotides/confidence are in reverse-complement coordinates of a
// chromosome of length `length`. Genes are mapped to forward coordinates and tagged
// strand "-"; the reported sequence is the coding (mRNA-sense) sequence.
void collect_minus_genes(
    vector<Gene_Prediction>& genes,
    const vector<Nucleotide>& chromosome_nucleotides,
    const vector<State>& chromosome_states,
    const vector<double>& chromosome_confidence,
    const string& scaffold,
    size_t length,
    const Transition_Model::Log_Prob_Matrix& transition_log_probs,
    const Emission_Model& emission_model)
{
    size_t pos = 0;
    while(pos < chromosome_states.size()){
        if(!is_gene(chromosome_states[pos])){
            pos++;
            continue;
        }

        size_t gene_start = pos;
        while(pos < chromosome_states.size() && is_gene(chromosome_states[pos])){
            pos++;
        }
        size_t gene_end = pos;

        size_t fwd_start = length - gene_end; // 0-based inclusive (forward)
        size_t fwd_end = length - gene_start; // 0-based exclusive (forward)

        json prediction;
        prediction["scaffold"] = scaffold;
        prediction["strand"] = "-";
        prediction["start"] = fwd_start + 1;
        prediction["end"] = fwd_end;
        prediction["exons"] = map_intervals_to_forward(
            intervals_for_family(chromosome_states, 0, gene_start, gene_end, is_coding), length);
        prediction["introns"] = map_intervals_to_forward(
            intervals_for_family(chromosome_states, 0, gene_start, gene_end, is_intron), length);
        prediction["intron_count"] = intervals_for_family(
            chromosome_states, 0, gene_start, gene_end, is_intron_body).size();
        prediction["sequence"] = sequence_preview(chromosome_nucleotides, gene_start, gene_end);
        prediction["confidence"] = confidence_for_range_minus(
            chromosome_states, chromosome_confidence, gene_start, gene_end, length);

        Log_Prob score = Viterbi::path_log_prob(
            chromosome_states, chromosome_nucleotides, transition_log_probs,
            emission_model, gene_start_penalty, gene_start, gene_end);

        genes.push_back({fwd_start, fwd_end, score, prediction});
    }
}

// Resolve overlapping genes (only possible across strands) by keeping the higher
// Viterbi path log-probability, mirroring how a single dual-strand HMM would let
// likelihood pick one strand per locus.
vector<Gene_Prediction> merge_strand_genes(vector<Gene_Prediction> genes) {
    sort(genes.begin(), genes.end(), [](const Gene_Prediction& a, const Gene_Prediction& b){
        if(a.fwd_start != b.fwd_start) return a.fwd_start < b.fwd_start;
        return a.fwd_end < b.fwd_end;
    });

    vector<Gene_Prediction> kept;
    for(auto& gene : genes){
        bool drop = false;
        while(!kept.empty() && kept.back().fwd_end > gene.fwd_start){
            if(gene.score > kept.back().score){
                kept.pop_back();
            } else {
                drop = true;
                break;
            }
        }
        if(!drop){
            kept.push_back(gene);
        }
    }
    return kept;
}

string value_after_arg(int argc, char** argv, const string& name, const string& fallback) {
    for(int i = 1; i + 1 < argc; i++){
        if(argv[i] == name){
            return argv[i + 1];
        }
    }
    return fallback;
}

} // namespace

namespace gene_hmm {
    extern Genome_Profile profile;
}

int main(int argc, char** argv) {
    try {
        string input_fna = value_after_arg(argc, argv, "--fna", "");
        string profile_path = value_after_arg(argc, argv, "--profile", "");
        string splice_cnn_scores_path = value_after_arg(argc, argv, "--splice-cnn-scores", "");
        string start_cnn_scores_path = value_after_arg(argc, argv, "--start-cnn-scores", "");
        string splice_cnn_scores_minus_path = value_after_arg(argc, argv, "--splice-cnn-scores-minus", "");
        string start_cnn_scores_minus_path = value_after_arg(argc, argv, "--start-cnn-scores-minus", "");
        if(input_fna.empty()){
            throw runtime_error("--fna PATH is required.");
        }
        if(profile_path.empty()){
            throw runtime_error("--profile PATH is required.");
        }

        gene_hmm::profile = Genome_Profile::load(profile_path);

        Training_Data training = load_training_data(gene_hmm::profile);
        vector<Chromosome_Range> train_ranges = split_usable_ranges(training.chromosomes, training.regions);

        auto transition_log_probs = Transition_Model::compute_log_probs(training.states, train_ranges);
        Emission_Model emission_model = train_emissions(training.states, training.nucleotides, train_ranges);

        vector<Nucleotide> input_nucleotides = FNA_Parser::parse_sequence(input_fna);
        vector<Chromosome_Range> input_ranges = FNA_Parser::get_chromosome_ranges(input_fna);
        if(splice_cnn_scores_path.empty()){
            throw runtime_error("--splice-cnn-scores PATH is required for CNN donor/acceptor emissions.");
        }
        emission_model.load_splice_cnn_scores(splice_cnn_scores_path, input_nucleotides.size());
        emission_model.set_splice_cnn_calibration(
            gene_hmm::profile.splice_cnn.donor_scale,
            gene_hmm::profile.splice_cnn.donor_bias,
            gene_hmm::profile.splice_cnn.acceptor_scale,
            gene_hmm::profile.splice_cnn.acceptor_bias);

        // Translation-start CNN scores are required for START_CODON_1 emissions.
        if(start_cnn_scores_path.empty()){
            throw runtime_error("--start-cnn-scores PATH is required for CNN start emissions.");
        }
        emission_model.load_start_cnn_scores(start_cnn_scores_path, input_nucleotides.size());
        emission_model.set_start_cnn_calibration(
            gene_hmm::profile.start_cnn.start_scale,
            gene_hmm::profile.start_cnn.start_bias);

        // Optional minus-strand decode: requires the profile flag plus revcomp CNN
        // score files for this input. The emission model is copied so the minus pass
        // reads its own scores without disturbing the plus pass.
        bool dual_strand = gene_hmm::profile.include_minus_strand &&
            !splice_cnn_scores_minus_path.empty() &&
            !start_cnn_scores_minus_path.empty();
        Emission_Model emission_model_minus = emission_model;
        if(dual_strand){
            emission_model_minus.load_splice_cnn_scores(splice_cnn_scores_minus_path, input_nucleotides.size());
            emission_model_minus.set_splice_cnn_calibration(
                gene_hmm::profile.splice_cnn.donor_scale,
                gene_hmm::profile.splice_cnn.donor_bias,
                gene_hmm::profile.splice_cnn.acceptor_scale,
                gene_hmm::profile.splice_cnn.acceptor_bias);
            emission_model_minus.load_start_cnn_scores(start_cnn_scores_minus_path, input_nucleotides.size());
            emission_model_minus.set_start_cnn_calibration(
                gene_hmm::profile.start_cnn.start_scale,
                gene_hmm::profile.start_cnn.start_bias);
        }

        json result;
        result["summary"] = {
            {"inputFile", input_fna.substr(input_fna.find_last_of("/\\") + 1)},
            {"totalBases", input_nucleotides.size()},
            {"scaffolds", input_ranges.size()},
            {"genes", 0},
            {"exons", 0},
            {"introns", 0}
        };
        result["scaffolds"] = json::array();
        result["predictions"] = json::array();
        result["confidenceByScaffold"] = json::object();

        vector<Gene_Prediction> all_genes;
        for(const auto& range : input_ranges){
            emission_model.set_splice_cnn_position_offset(range.start);
            emission_model.set_start_cnn_position_offset(range.start);
            vector<Nucleotide> chromosome_nucleotides = slice_nucleotides(input_nucleotides, range.start, range.end);
            vector<State> chromosome_states = Viterbi::decode(
                chromosome_nucleotides,
                transition_log_probs,
                emission_model,
                0,
                numeric_limits<size_t>::max(),
                gene_start_penalty);
            vector<double> chromosome_confidence = Forward_Backward::confidence(
                chromosome_nucleotides,
                chromosome_states,
                transition_log_probs,
                emission_model,
                gene_start_penalty);

            result["scaffolds"].push_back({
                {"name", range.name},
                {"length", range.end - range.start}
            });
            result["confidenceByScaffold"][range.name] = confidence_for_range(
                chromosome_states,
                chromosome_confidence,
                0,
                0,
                chromosome_states.size());

            vector<Gene_Prediction> scaffold_genes;
            collect_plus_genes(
                scaffold_genes,
                chromosome_nucleotides,
                chromosome_states,
                chromosome_confidence,
                range.name,
                transition_log_probs,
                emission_model);

            if(dual_strand){
                emission_model_minus.set_splice_cnn_position_offset(range.start);
                emission_model_minus.set_start_cnn_position_offset(range.start);
                vector<Nucleotide> revcomp_nucleotides = FNA_Parser::reverse_complement(chromosome_nucleotides);
                vector<State> minus_states = Viterbi::decode(
                    revcomp_nucleotides,
                    transition_log_probs,
                    emission_model_minus,
                    0,
                    numeric_limits<size_t>::max(),
                    gene_start_penalty);
                vector<double> minus_confidence = Forward_Backward::confidence(
                    revcomp_nucleotides,
                    minus_states,
                    transition_log_probs,
                    emission_model_minus,
                    gene_start_penalty);
                collect_minus_genes(
                    scaffold_genes,
                    revcomp_nucleotides,
                    minus_states,
                    minus_confidence,
                    range.name,
                    revcomp_nucleotides.size(),
                    transition_log_probs,
                    emission_model_minus);
            }

            vector<Gene_Prediction> merged = merge_strand_genes(scaffold_genes);
            all_genes.insert(all_genes.end(), merged.begin(), merged.end());
        }

        for(auto& gene : all_genes){
            string id = "pred_" + string(4 - min<size_t>(4, to_string(result["predictions"].size() + 1).size()), '0') +
                        to_string(result["predictions"].size() + 1);
            gene.payload["id"] = id;
            result["predictions"].push_back(gene.payload);
        }

        size_t exon_count = 0;
        size_t intron_count = 0;
        for(const auto& prediction : result["predictions"]){
            exon_count += prediction["exons"].size();
            intron_count += prediction["intron_count"].get<size_t>();
        }
        result["summary"]["genes"] = result["predictions"].size();
        result["summary"]["exons"] = exon_count;
        result["summary"]["introns"] = intron_count;

        cout << result.dump();
        return 0;
    } catch(const exception& error) {
        cerr << error.what() << "\n";
        return 1;
    }
}

