#include "globals.h"
#include "neuron.h"

const std::map<std::string, NeuronParams> Neuron::NeuronTypes{
    {"RS",  {.02, .2,  -65, 8,   false}},
    {"IB",  {.02, .2,  -55, 4,   false}},
    {"CH",  {.02, .2,  -50, 2,   false}},
    {"FS",  {.1,  .2,  -65, 2,   true }},
    {"LTS", {.02, .25, -65, 2,   true }},
    {"TC",  {.02, .25, -65, .05, false}},
    {"RZ",  {.1,  .26, -65, 2,   false}}
};

double Neuron::firing_thresh = _FIRING_TH_;

Neuron::Neuron() : _poten(_REST_VAL_) {
    set_type("RS");
}

void Neuron::set_params(const NeuronParams &np, double noise) {
    params = np;
    if (std::abs(noise)>1e-8) {
        if (params.inhib) {
            params.a *= 1-_AVAR_*noise;
            params.b *= 1+_BVAR_*noise;
        } else {
            noise *= noise;
            params.c *= 1-_CVAR_*noise;
            params.d *= 1-_DVAR_*noise;
        }
    }
    _recov = params.b*_poten;
}

void Neuron::set_type(std::string typ) {    
    if (!type_exists(typ)) typ = "RS";
    _type = NeuronTypes.find(typ);
}

void Neuron::set_default_params(const std::string &typ, const double noise) {
    set_type(typ);
    set_params(_type->second, noise);
}

bool Neuron::is_type(const std::string &_t) {
    return (_type != NeuronTypes.end()) && (_type->first == _t);
}

void Neuron::step() {
    _poten += 0.5*(0.04*_poten*_poten+5*_poten+140-_recov+_input);
    _poten += 0.5*(0.04*_poten*_poten+5*_poten+140-_recov+_input);
    _recov += params.a*(params.b*_poten-_recov);

}

std::string Neuron::formatted_params() const {
    std::stringstream ss;
    ss << _type->first << '\t'
       << params.a << '\t'
       << params.b << '\t'
       << params.c << '\t'
       << params.d << '\t'
       << (int)params.inhib;
    return ss.str();
}

std::string Neuron::formatted_values() const {
    std::stringstream ss;
    ss << _poten << '\t' << _recov << '\t' << _input;
    return ss.str();
}
