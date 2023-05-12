#include "models.hh"

#include "Chili/Model/Model.hh"
#include "Chili/Channel/Integrand.hh"
#include "Chili/Channel/ChannelNode.hh"
#include "Chili/Channel/Channel.hh"
#include "Chili/Channel/MultiChannel.hh"

#include <cmath>

using namespace bcnutau;

std::vector<int> combine(int i, int j) {
    // Handle gluon-quark vertex
    if(i == 21 && std::abs(j) < 7) return { j };
    else if(j == 21 && std::abs(i) < 7) return { i };
    else if(std::abs(i) < 7 && i == -j) return { 21 };
    // Handle W vertices
    else if(i == 24 && j == -4) return { -5 };
    else if(i == 24 && j == 5) return { 4 };
    else if(j == 24 && i == -4) return { -5 };
    else if(j == 24 && i == 5) return { 4 };
    else if(i == -24 && j == 4) return { 5 };
    else if(i == -24 && j == -5) return { -4 };
    else if(j == -24 && i == 4) return { 5 };
    else if(j == -24 && i == -5) return { -4 };

    else if(i == 24 && j == 15) return { 16 };
    else if(i == 24 && j == -16) return { -15 };
    else if(j == 24 && i == 15) return { 16 };
    else if(j == 24 && i == -16) return { -15 };
    else if(i == -24 && j == -15) return { -16 };
    else if(i == -24 && j == 16) return { 15 };
    else if(j == -24 && i == -15) return { -16 };
    else if(j == -24 && i == 16) return { 15 };

    else if(i == 4 && j == -5) return { 24 };
    else if(j == 4 && i == -5) return { 24 };
    else if(i == -4 && j == 5) return { -24 };
    else if(j == -4 && i == 5) return { -24 };
    else if(i == 15 && j == -16) return { -24 };
    else if(j == 15 && i == -16) return { -24 };
    else if(i == -15 && j == 16) return { 24 };
    else if(j == -15 && i == 16) return { 24 };

    // Handle BSM
    else if(i == 9999 && j == 4) return { -15 };
    else if(j == 9999 && i == 4) return { -15 };
    else if(i == 9999 && j == 15) return { -4 };
    else if(j == 9999 && i == 15) return { -4 };
    else if(i == 9999 && j == 5) return { -16 };
    else if(j == 9999 && i == 5) return { -16 };
    else if(i == 9999 && j == 16) return { -5 };
    else if(j == 9999 && i == 16) return { -5 };
    else if(i == 9999 && j == -4) return { 15 };
    else if(j == 9999 && i == -4) return { 15 };
    else if(i == 9999 && j == -15) return { 4 };
    else if(j == 9999 && i == -15) return { 4 };
    else if(i == 9999 && j == -5) return { 16 };
    else if(j == 9999 && i == -5) return { 16 };
    else if(i == 9999 && j == -16) return { 5 };
    else if(j == 9999 && i == -16) return { 5 };

    else if(i == 4 && j == -15) return { 9999 };
    else if(j == 4 && i == -15) return { 9999 };
    else if(i == 5 && j == -16) return { 9999 };
    else if(j == 5 && i == -16) return { 9999 };
    else if(i == -4 && j == 15) return { 9999 };
    else if(j == -4 && i == 15) return { 9999 };
    else if(i == -5 && j == 16) return { 9999 };
    else if(j == -5 && i == 16) return { 9999 };
    
    // Default (error)
    return { };
}

bool PreProcess(const std::vector<chili::FourVector> &mom) {
    if(std::isnan(mom[0][0])) {
        return false;
    }
    if((mom[0]+mom[1]).Mass2() > 13000*13000) {
        return false;
    }
    return true;
}

bool PostProcess(const std::vector<chili::FourVector>&, double) { return true; }

int main() {
    chili::Model model(combine);
    model.Mass(4) = 0;
    model.Mass(-4) = 0;
    model.Mass(5) = 172;
    model.Mass(-5) = 172;
    model.Mass(15) = 0;
    model.Mass(-15) = 0;
    model.Mass(16) = 0;
    model.Mass(-16) = 0;
    model.Mass(24) = 80.385;
    model.Mass(21) = 0;
    model.Mass(9999) = 100;
    model.Width(4) = 0;
    model.Width(-4) = 0;
    model.Width(5) = 0;
    model.Width(-5) = 0;
    model.Width(15) = 0;
    model.Width(-15) = 0;
    model.Width(16) = 0;
    model.Width(-16) = 0;
    model.Width(24) = 2.045;
    model.Width(21) = 0;
    model.Width(9999) = 1;

    spdlog::set_level(spdlog::level::info);
    std::vector<int> process{21, 4, 5, 16, -15};

    // Setup integrator
    auto mappings = chili::ConstructChannels(13000, process, model, 2);
    std::cout << mappings.size() << std::endl;
    chili::Integrand<chili::FourVector> integrand;
    for(auto &mapping : mappings) {
        chili::Channel<chili::FourVector> channel;
        channel.mapping = std::move(mapping);
        // Initializer takes the number of integration dimensions
        // and the number of bins for vegas to start with
        chili::AdaptiveMap map(channel.mapping -> NDims(), 2);
        // Initializer takes adaptive map and settings (found in struct VegasParams)
        channel.integrator = chili::Vegas(map, chili::VegasParams{});
        integrand.AddChannel(std::move(channel));
    }

    // Initialize the multichannel integrator
    // Takes the number of dimensions, the number of channels, and options
    // The options can be found in the struct MultiChannelParams
    chili::MultiChannel integrator{integrand.NDims(), integrand.NChannels(), {}};

    // To integrate a function you need to pass it in and tell it to optimize
    // Summary will print out a summary of the results including the values of alpha
    bcnutau::StandardModel model_sm("CT18NNLO"); ;
    auto func = [&](const std::vector<chili::FourVector> &mom) {
        return model_sm.Evaluate(mom);
    };
    integrand.Function() = func;
    integrand.PreProcess() = PreProcess;
    integrand.PostProcess() = PostProcess;
    spdlog::info("Starting optimization");
    integrator.Optimize(integrand);
    integrator.Summary(std::cout);

    spdlog::info("Saving trained integrator");
    integrator.SaveAs(integrand);

    return 0;
}
