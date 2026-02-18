#pragma once

#include "emp/bits/Bits.hpp"
#include "emp/math/Random.hpp"

class Organism {
private:
    emp::BitVector genome{};
    double fitness = 0.0; // normalized and maximized
    double mut_rate = 0.0; // per gene

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
    Organism Mutate(emp::Random & random) const {
        emp::BitVector new_genome = genome;
        for (size_t i = 0; i < new_genome.size(); ++i) {
            if (random.P(mut_rate)) {
                new_genome.Toggle(i);
            }
        }

        Organism offspring(new_genome);
        offspring.SetFitness(0.0); // remove parent's fitness
        offspring.SetMutationRate(random.GetDouble());
        return offspring;
    }

};