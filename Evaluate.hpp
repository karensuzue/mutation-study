#pragma once

#include "Organism.hpp"

struct Evaluate {
    using phenotype_t = emp::vector<double>;
    using genome_t = emp::vector<double>;

    /* ------- GENOTYPE TO PHENOTYPE TRANSLATORS -------*/
    phenotype_t ExploitationRate(const genome_t & g) {
        phenotype_t p = g;
        return p;
    }

    /*
     * Taken from: 
     * https://github.com/jgh9094/ECJ-2023-Suite-Of-Diagnostic-Metrics-For-Characterizing-Selection-Schemes
     * 
     * Multi-valley Crossing function
     *
     * Solutions are pressured to cross valleys of different widths at each gene.
     * We use a peaks vector that is supplied to calculate the penalty at each valley.
     * We use the difference between the gene value and penalty value to determine
     * the reduction when calculating the trait value.
     *
     * phenotype[i] = peaks[ floor(genome[i]) ] -  (genome[i] - peaks[ floor(genome[i]) ])
     *
     * @param g genome from organism being evaluated.
     * @param peaks vector of peaks for each floored gene value
     * @param dips_start gene value where peaks begin
     * @param dips_end gene value where dips end
     *
     * @return phenotype vector that is calculated from 'g'.
     */
    phenotype_t MultiValleyCrossing(const phenotype_t & p, 
                                    const phenotype_t & peaks, 
                                    const double & dips_start, 
                                    const double & dips_end) {

        // Quick checks
        emp_assert(peaks.size() > 0);

        // Intialize vector with size p
        phenotype_t phenotype(p.size());
        for(size_t i = 0; i < p.size(); ++i)
        {
            if (p[i] <= dips_start || p[i] >= dips_end)
            {
                phenotype[i] = p[i];
            }
            else
            {
                phenotype[i] = 2.0 * peaks[static_cast<size_t>(p[i])] - p[i];
            }
        }
        return phenotype;
    }

    /* ------- FITNESS CALCULATORS -------*/
    double AggregateFitness(const phenotype_t & p, size_t start_idx=0) {
        assert(!p.empty());
        assert(start_idx < p.size());
        double fitness = std::accumulate(p.begin() + start_idx, 
                                  p.end(), 0.0);
        return fitness;
    }

    // Calculate the sum of squared errors across each gene
    double SquaredErrorFitness(const phenotype_t & g, const phenotype_t & target_g) {
        assert(g.size() == target_g.size());

        double error_sum = 0.0;
        for (size_t i = 0; i < g.size(); ++i) {
            double diff = target_g[i] - g[i];
            error_sum += diff * diff;
        }
        return -error_sum;
    }
};