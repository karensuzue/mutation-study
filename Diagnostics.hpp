/*
 * Fitness evaluators taken from 
 *      Hernandez, J. G., Lalejini, A., & Ofria, C. (2022).
 *      A suite of diagnostic metrics for characterizing selection schemes
 */
#pragma once

#include "DiagOrganism.hpp"

struct Diagnostics {
    using phenotype_t = emp::vector<double>;
    using genome_t = emp::vector<double>;

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
};