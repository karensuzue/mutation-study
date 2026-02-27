#pragma once

#include "emp/base/vector.hpp"

#include "DiagOrganism.hpp"

struct DiagSelector {
    virtual ~Selector() = default;
    // Define one round of selection, return an index
    virtual size_t Select(const emp::vector<DiagOrganism> &, emp::Random &) = 0;
};


struct TournamentSelector : DiagSelector {
    size_t tour_size;
    TournamentSelector(size_t tour_size) : tour_size(tour_size) {}

    // Select a single parent 
    size_t Select(const emp::vector<DiagOrganism> & pop, emp::Random & random) override {
        size_t best_idx = random.GetUInt(pop.size());

        // We already picked one, so start at i = 1
        for (size_t i = 1; i < tour_size; ++i) {
            size_t cand_idx = random.GetUInt(pop.size());
            if (pop[cand_idx].GetFitness() > pop[best_idx].GetFitness()) best_idx = cand_idx;
        }

        return best_idx;
    }
};