#include <iostream>

#include "DiagOrganism.hpp"
#include "DiagPopulation.hpp"

int main() {
    emp::vector<double> genome = {0, 0, 10, 9, 8};
    DiagOrganism org(genome);
    org.SetPhenotype(genome);
    org.SetStartIndex(2);
    std::cout << org;
    org.UpdateFitnessFromPhenotype();
    std::cout << org;
}