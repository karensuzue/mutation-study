/*
 * Fitness evaluators taken from 
 *      Hernandez, J. G., Lalejini, A., & Ofria, C. (2022),
 *      A suite of diagnostic metrics for characterizing selection schemes
*/
#pragma once

#include <cassert>
#include <cstdint>

#include "DiagOrganism.hpp"

struct DiagFitness {
    virtual ~Fitness() = default;
    virtual double Evaluate(const DiagOrganism &) const = 0; // higher is better!
};

// Try to set every gene to 1
struct OneMaxFitness : Fitness {
    double Evaluate(const DiagOrganism & org) const override {
        return static_cast<double>(org.GetGenome().CountOnes());
    }
};
