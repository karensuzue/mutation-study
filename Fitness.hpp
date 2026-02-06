#pragma once

#include <cassert>

#include "Organism.hpp"

struct Fitness {
    virtual ~Fitness() = default;
    virtual double Evaluate(const Organism &) const = 0; // higher is better!
};

// Try to set every gene to 1
struct OneMaxFitness : Fitness {
    double Evaluate(const Organism & org) const override {
        return static_cast<double>(org.GetGenome().CountOnes());
    }
};

// Given a target unsigned integer, make sure genes sum up to it
struct IntMatchFitness : Fitness {
    size_t target;
    size_t genome_len;

    IntMatchFitness(size_t target, size_t genome_len) 
        : target(target), genome_len(genome_len) {
            assert(genome_len > 0 && "Genome length can not be 0!");
        }

    double Evaluate(const Organism & org) const override {
        const size_t ones = org.GetGenome().CountOnes();
        const size_t dist = (ones > target) ? (ones - target) : (target - ones);
        // normalize (best = 1, worst = 0)
        return 1.0 - static_cast<double>(dist) / static_cast<double>(genome_len);
    }
};