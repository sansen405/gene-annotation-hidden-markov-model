#pragma once

#include "../topology/Topology.hpp"
#include "Emission_Model.hpp"
#include "Transition_Model.hpp"
#include <cstddef>
#include <string>

namespace gene_hmm {

    struct Trained_Model {
        std::string genome_name;
        Transition_Model::Log_Prob_Matrix transition_log_probs{};
        Emission_Model emission_model;
        size_t max_intron_body_length = 0;
        double gene_start_penalty = 1.0;
    };

    class Model_IO {
    public:
        static void save(const std::string& path, const Trained_Model& model);
        static Trained_Model load(const std::string& path);
    };

} // namespace gene_hmm
