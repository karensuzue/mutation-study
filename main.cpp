#include <iostream>

#include "emp/tools/String.hpp"

#include "Organism.hpp"
#include "Population.hpp"

int main(/*int argc, char * argv[]*/) {
    emp::Random random{0};

    Population pop;
    pop.InitializeUniform(random);
    pop.SetPopulationMutation(0.5);
    pop.EvaluateFitness();
    std::cout << pop << std::endl;
    // if (!pop.IsConstantMutation()) pop.ToggleConstantMutation();
    pop.RunOneGeneration(random);
    std::cout << pop << std::endl;
    
    // pop.MultiRun();
}