#include "../genome_profiles/Genome_Profile.hpp"
#include <nlohmann/json.hpp>
#include <string>
#include <fstream>
#include <stdexcept>

using json = nlohmann::json;

namespace gene_hmm {

    Genome_Profile profile;

    namespace {
        vector<string> parse_string_list(const json& value) {
            vector<string> paths;
            if (value.is_string()) {
                paths.push_back(value.get<string>());
                return paths;
            }
            for (const auto& item : value) {
                paths.push_back(item.get<string>());
            }
            return paths;
        }

        Species_Dataset parse_species_dataset(const json& dataset) {
            Species_Dataset species;
            species.name = dataset.value("name", "");
            species.source_fasta_path = dataset["source_fasta"];
            species.source_gff_path   = dataset["source_gff"];
            species.train_fasta_path  = dataset["train_fasta"];
            species.train_gff_path    = dataset["train_gff"];
            species.test_fasta_path   = dataset["test_fasta"];
            species.test_gff_path     = dataset["test_gff"];
            for (const auto& c : dataset.value("test_chromosomes", json::array()))
                species.test_chromosomes.push_back(c);
            for (const auto& c : dataset.value("excluded_chromosomes", json::array()))
                species.excluded_chromosomes.push_back(c);
            return species;
        }
    }

    Genome_Profile Genome_Profile::load(const string& json_path) {
        ifstream file(json_path);
        if (!file.is_open())
            throw runtime_error("Cannot open genome profile: " + json_path);

        json j = json::parse(file);

        Genome_Profile p;

        p.name = j["name"];

        const auto& dataset = j["dataset"];
        if (dataset.contains("species")) {
            for (const auto& species_json : dataset["species"])
                p.species.push_back(parse_species_dataset(species_json));
        } else {
            p.species.push_back(parse_species_dataset(dataset));
        }

        if (!p.species.empty()) {
            const auto& first = p.species.front();
            p.source_fasta_path = first.source_fasta_path;
            p.source_gff_path   = first.source_gff_path;
            p.train_fasta_path  = first.train_fasta_path;
            p.train_gff_path    = first.train_gff_path;
            p.test_fasta_path   = first.test_fasta_path;
            p.test_gff_path     = first.test_gff_path;
            p.test_chromosomes  = first.test_chromosomes;
            p.excluded_chromosomes = first.excluded_chromosomes;
        }

        if (j.contains("splice_cnn")) {
            const auto& splice_cnn = j["splice_cnn"];
            p.splice_cnn.model_path = splice_cnn.value("model", "");
            if (splice_cnn.contains("train_scores")) {
                p.splice_cnn.train_score_paths = parse_string_list(splice_cnn["train_scores"]);
            }
            if (splice_cnn.contains("test_scores")) {
                p.splice_cnn.test_score_paths = parse_string_list(splice_cnn["test_scores"]);
            }
        }

        p.min_first_cds_bp     = j["filters"]["min_first_cds_bp"];
        p.min_last_cds_bp      = j["filters"]["min_last_cds_bp"];
        p.min_intron_bp        = j["filters"]["min_intron_bp"];
        p.require_3n_cds       = j["filters"]["require_3n_cds"];
        p.include_minus_strand = j["filters"]["include_minus_strand"];

        for (const auto& [name, val] : j["emissions"].items()) {
            Emission_Config cfg;
            string type = val["type"];
            if (type == "markov") {
                cfg.type       = Emission_Family::MARKOV;
                cfg.order      = val.value("order",      1);
                cfg.frame_tied = val.value("frame_tied", false);
            } else if (type == "pssm") {
                cfg.type         = Emission_Family::PSSM;
                cfg.window_left  = val["window_left"];
                cfg.window_right = val["window_right"];
            } else {
                cfg.type = Emission_Family::DETERMINISTIC;
            }
            p.emissions[name] = cfg;
        }

        p.transition_alpha = j["smoothing"]["transition_alpha"];
        p.emission_alpha   = j["smoothing"]["emission_alpha"];

        return p;
    }

}
