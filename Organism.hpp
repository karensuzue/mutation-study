#pragma once

#include <algorithm>
#include <cmath>

#include "emp/bits/Bits.hpp"
#include "emp/math/Random.hpp"

// Minimum and maximum mutation rates
static constexpr double MUT_MIN = 1e-12; // if 0 can't be shifted
static constexpr double MUT_MAX = 1.0;

class Organism {
private:
    emp::BitVector genome{};
    double fitness = 0.0; // normalized and maximized
    double mut_rate = MUT_MIN; // per gene

public:
    Organism() = default;
    Organism(const Organism &) = default;

    Organism(const emp::BitVector & init_genome)
      : genome(init_genome) {}

    Organism(size_t genome_size)
      : genome(emp::BitVector (genome_size)) {}

    friend std::ostream & operator<<(std::ostream & os, const Organism & org) {
        os << "Genome=" << org.GetGenome() 
           << ", Fitness=" << org.fitness
           << ", Mutation=" << org.mut_rate;
        return os;
    }

    emp::BitVector & GetGenome() { return genome; }
    const emp::BitVector & GetGenome() const { return genome; }

    double GetFitness() const { return fitness; }
    void SetFitness(double f) { fitness = f; }

    double GetMutationRate() const { return mut_rate; }
    void SetMutationRate(double m) { mut_rate = m; }

    // Return a mutated COPY
    Organism Mutate(emp::Random & random, bool const_mutation,
                    double sigma=0.138, double b=0.0, double pi=0.8) const {
        // Mutate genome using current mutation rate
        emp::BitVector new_genome = genome;
        for (size_t i = 0; i < new_genome.size(); ++i) {
            if (random.P(mut_rate)) new_genome.Toggle(i);
        }

        // Inherit parent's mutation rate by default
        double mut_child = mut_rate;

        // With probability 'pi', optionally shift mutation rate (mutator mutation)
        if (!const_mutation && random.P(pi)) { 
            double mean = b * sigma * sigma; // for upward bias, if b > 0

            // x = log2(mut_child / mut_parent) ~ N(mean, sigma^2)
            // For sigma=0.138, a 1-sigma draw is about a 1.10x increase (or 0.91x decrease)
            double x = random.GetNormal(mean, sigma);
            mut_child *= std::pow(2.0, x); // multiplicative shift
            mut_child = std::clamp(mut_child, MUT_MIN, MUT_MAX);
        }

        Organism offspring(new_genome);
        offspring.SetFitness(0.0); // remove parent's fitness
        offspring.SetMutationRate(mut_child);
        return offspring;
    }

};