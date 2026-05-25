#include "../genome_profiles/Genome_Profile.hpp"
#include <nlohmann/json.hpp>
#include <string>
#include <fstream>
#include <stdexcept>

using json = nlohmann::json;

namespace gene_hmm {

    Genome_Profile profile;

    Genome_Profile Genome_Profile::load(const string& json_path) {
        ifstream file(json_path);
        if (!file.is_open())
            throw runtime_error("Cannot open genome profile: " + json_path);

        json j = json::parse(file);

        Genome_Profile p;

        p.name = j["name"];

        const auto& dataset = j["dataset"];
        p.source_fasta_path = dataset["source_fasta"];
        p.source_gff_path   = dataset["source_gff"];
        p.train_fasta_path  = dataset["train_fasta"];
        p.train_gff_path    = dataset["train_gff"];
        p.test_fasta_path   = dataset["test_fasta"];
        p.test_gff_path     = dataset["test_gff"];
        for (const auto& c : dataset["test_chromosomes"])
            p.test_chromosomes.push_back(c);
        for (const auto& c : dataset["excluded_chromosomes"])
            p.excluded_chromosomes.push_back(c);

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
