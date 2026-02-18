#include <iostream>

#include "emp/tools/String.hpp"

#include "Organism.hpp"
#include "Population.hpp"

int main(/*int argc, char * argv[]*/) {
    emp::Random random{11};

    Population pop;
    pop.InitializeUniform(random);
    pop.MultiRun();
}