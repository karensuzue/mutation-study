#include <iostream>
#include <sstream>
#include <string>

#include "emp/base/vector.hpp"

#include "Organism.hpp"
#include "Population.hpp"

// TODO: SettingsManager!
int main(/*int argc, char * argv[]*/) { 

    // EXPERIMENT 1 - CONSTANT MUTATION
    // U = mutations per genome per generation (genome-wide rate)
    const emp::vector<double> U_rates = {
        1.0e-5,
        2.6827e-5,
        7.1969e-5,
        1.9307e-4,
        5.1795e-4,
        1.3895e-3,
        3.7276e-3,
        1.0e-2,
        2.6827e-2,
        7.1969e-2,
        1.9307e-1,
        5.1795e-1,
        1.3895,
        3.7276,
        10.0
    };

    for (double  U : U_rates) {
        Population pop;

        // Ensure contant mutation-rate inheritance
        if (!pop.IsConstantMutation()) pop.ToggleConstantMutation();

        // Convert genome-wide U to per-site mu
        const double mu = U / static_cast<double>(pop.GetGenomeSize());
        pop.SetInitMutation(mu);

        // Formatting
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(4) << U;
        const std::string tag = oss.str();

        std::cout << "Running U=" << tag << " (mu=" 
            << std::scientific << std::setprecision(4) << mu << ")\n";
        pop.MultiRun(tag);
    }

}