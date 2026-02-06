#pragma once


#include "emp/base/vector.hpp"

#include "Organism.hpp"

struct Selector {
    virtual ~Selector() = default;
    // Define one round of selection, return an index
    virtual size_t Select(const emp::vector<Organism> &, emp::Random &) = 0;
};


struct TournamentSelector : Selector {
    size_t tour_size;
    TournamentSelector(size_t tour_size) : tour_size(tour_size) {}

    // Select a single parent 
    size_t Select(const emp::vector<Organism> & pop, emp::Random & random) override {
        size_t best_idx = random.GetUInt(pop.GetSize());

        // We already picked one, so start at i = 1
        for (int i = 1; i < tour_size; ++i) {
            size_t cand_idx = random.GetUInt(pop.GetSize());
            if (pop[cand_idx].GetFitness() > pop[best_idx].GetFitness()) best_idx = cand_idx;
        }

        return best_idx;
    }
}