#include <iostream>

#include "emp/tools/String.hpp"

#include "Organism.hpp"
#include "Population.hpp"

int main(/*int argc, char * argv[]*/) {
    Organism test_org(5);
    std::cout << test_org << std::endl;

    emp::BitVector genome(5);
    emp::Random random;
    for (size_t i = 0; i < genome.GetSize(); ++i) {
        genome.Set(i, random.P(0.5)); // 50-50
    }
    Organism test_org_2(genome);
    std::cout << test_org_2 << std::endl;

    // emp::String config_name = "Mutation.cfg"; // default
    // if (argc > 1) config_name = argv[1]; // through command line input
    // emp::SettingsManager settings;

    // bool success = settings.Load(config_name);
    // if (!success) {
    //     emp::PrintLn(settings.GetError());
    //     exit(1);
    // }

    Population pop;
    // pop.SetupConfig(settings);
    // pop.InitializeUniform();
    pop.InitializeRandom();
    pop.PrintPopulation(std::cout);
    

}