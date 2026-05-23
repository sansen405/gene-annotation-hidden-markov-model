#include "Model_IO.hpp"
#include <nlohmann/json.hpp>
#include <cmath>
#include <fstream>
#include <stdexcept>

using json = nlohmann::json;
using namespace gene_hmm;
using namespace std;

namespace {

// -inf log probs are stored as this sentinel so JSON stays valid.
constexpr double JSON_NEG_INF = -1e300;

double to_jval(double v) {
    return (isinf(v) && v < 0) ? JSON_NEG_INF : v;
}

double from_jval(double v) {
    return (v <= JSON_NEG_INF * 0.5) ? LOG_ZERO : v;
}

json markov1_to_json(const Emission_Model::Markov1_Log_Prob& t) {
    json j = json::array();
    for (const auto& row : t) {
        json r = json::array();
        for (double v : row) r.push_back(to_jval(v));
        j.push_back(r);
    }
    return j;
}

Emission_Model::Markov1_Log_Prob markov1_from_json(const json& j) {
    Emission_Model::Markov1_Log_Prob t{};
    for (size_t i = 0; i < NUM_NUCLEOTIDES; ++i)
        for (size_t k = 0; k < NUM_NUCLEOTIDES; ++k)
            t[i][k] = from_jval(j[i][k].get<double>());
    return t;
}

json markov5_to_json(const Emission_Model::Markov5_Log_Prob& t) {
    json j = json::array();
    for (const auto& row : t) {
        json r = json::array();
        for (double v : row) r.push_back(to_jval(v));
        j.push_back(r);
    }
    return j;
}

Emission_Model::Markov5_Log_Prob markov5_from_json(const json& j) {
    Emission_Model::Markov5_Log_Prob t{};
    for (size_t i = 0; i < Emission_Model::MARKOV5_CONTEXTS; ++i)
        for (size_t k = 0; k < NUM_NUCLEOTIDES; ++k)
            t[i][k] = from_jval(j[i][k].get<double>());
    return t;
}

json pssm_to_json(const Emission_Model::PSSM_Log_Prob& t) {
    json j = json::array();
    for (const auto& row : t) {
        json r = json::array();
        for (double v : row) r.push_back(to_jval(v));
        j.push_back(r);
    }
    return j;
}

Emission_Model::PSSM_Log_Prob pssm_from_json(const json& j) {
    Emission_Model::PSSM_Log_Prob t;
    for (const auto& row : j) {
        array<double, NUM_NUCLEOTIDES> arr{};
        for (size_t k = 0; k < NUM_NUCLEOTIDES; ++k)
            arr[k] = from_jval(row[k].get<double>());
        t.push_back(arr);
    }
    return t;
}

json transition_to_json(const Transition_Model::Log_Prob_Matrix& m) {
    json j = json::array();
    for (const auto& row : m) {
        json r = json::array();
        for (double v : row) r.push_back(to_jval(v));
        j.push_back(r);
    }
    return j;
}

Transition_Model::Log_Prob_Matrix transition_from_json(const json& j) {
    Transition_Model::Log_Prob_Matrix m{};
    for (size_t i = 0; i < NUM_STATES; ++i)
        for (size_t k = 0; k < NUM_STATES; ++k)
            m[i][k] = from_jval(j[i][k].get<double>());
    return m;
}

} // namespace

namespace gene_hmm {

void Model_IO::save(const string& path, const Trained_Model& model) {
    json j;
    j["genome"]                 = model.genome_name;
    j["max_intron_body_length"] = model.max_intron_body_length;
    j["gene_start_penalty"]     = model.gene_start_penalty;
    j["transition_log_probs"]   = transition_to_json(model.transition_log_probs);

    json em;
    em["intergenic_lp"]  = markov1_to_json(model.emission_model.intergenic_lp);
    em["intron_lp"]      = markov1_to_json(model.emission_model.intron_lp);
    em["exon_lp"]        = markov5_to_json(model.emission_model.exon_lp);

    json frame_lp = json::array();
    for (size_t f = 0; f < 3; ++f)
        frame_lp.push_back(markov5_to_json(model.emission_model.exon_frame_lp[f]));
    em["exon_frame_lp"]  = frame_lp;

    em["start_codon_lp"] = pssm_to_json(model.emission_model.start_codon_lp);
    em["donor_lp"]       = pssm_to_json(model.emission_model.donor_lp);
    em["acceptor_lp"]    = pssm_to_json(model.emission_model.acceptor_lp);

    em["window_sizes"] = {
        {"start_left",     model.emission_model.start_window_left},
        {"start_right",    model.emission_model.start_window_right},
        {"donor_left",     model.emission_model.donor_window_left},
        {"donor_right",    model.emission_model.donor_window_right},
        {"acceptor_left",  model.emission_model.acceptor_window_left},
        {"acceptor_right", model.emission_model.acceptor_window_right}
    };
    j["emission"] = em;

    ofstream file(path);
    if (!file) throw runtime_error("Cannot open model file for writing: " + path);
    file << j.dump(2);
}

Trained_Model Model_IO::load(const string& path) {
    ifstream file(path);
    if (!file) throw runtime_error("Cannot open model file: " + path);

    json j;
    file >> j;

    Trained_Model model;
    model.genome_name             = j["genome"].get<string>();
    model.max_intron_body_length  = j["max_intron_body_length"].get<size_t>();
    model.gene_start_penalty      = j["gene_start_penalty"].get<double>();
    model.transition_log_probs    = transition_from_json(j["transition_log_probs"]);

    const json& em = j["emission"];
    model.emission_model.intergenic_lp  = markov1_from_json(em["intergenic_lp"]);
    model.emission_model.intron_lp      = markov1_from_json(em["intron_lp"]);
    model.emission_model.exon_lp        = markov5_from_json(em["exon_lp"]);

    for (size_t f = 0; f < 3; ++f)
        model.emission_model.exon_frame_lp[f] = markov5_from_json(em["exon_frame_lp"][f]);

    model.emission_model.start_codon_lp = pssm_from_json(em["start_codon_lp"]);
    model.emission_model.donor_lp       = pssm_from_json(em["donor_lp"]);
    model.emission_model.acceptor_lp    = pssm_from_json(em["acceptor_lp"]);

    const json& ws = em["window_sizes"];
    model.emission_model.start_window_left    = ws["start_left"].get<size_t>();
    model.emission_model.start_window_right   = ws["start_right"].get<size_t>();
    model.emission_model.donor_window_left    = ws["donor_left"].get<size_t>();
    model.emission_model.donor_window_right   = ws["donor_right"].get<size_t>();
    model.emission_model.acceptor_window_left = ws["acceptor_left"].get<size_t>();
    model.emission_model.acceptor_window_right= ws["acceptor_right"].get<size_t>();

    return model;
}

} // namespace gene_hmm
