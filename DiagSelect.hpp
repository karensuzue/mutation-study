#pragma once

#include "emp/base/vector.hpp"

#include "DiagOrganism.hpp"

struct DiagSelect {
    // Select a single parent index
    size_t Tournament(const emp::vector<DiagOrganism> & pop, 
                      size_t tour_size, 
                      emp::Random & random) {
        size_t best_idx = random.GetUInt(pop.size());
        // We already picked one, so start at i = 1
        for (size_t i = 1; i < tour_size; ++i) {
            size_t cand_idx = random.GetUInt(pop.size());
            if (pop[cand_idx].GetFitness() > pop[best_idx].GetFitness()) best_idx = cand_idx;
        }
        return best_idx;
    }
};
