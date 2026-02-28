/*
 * This class represents an organism with both phenotype and genotype representations,
 * consisting of a vector of genes, where each gene is a number from 0 to 100.
 * This representation is taken from:
 *      Hernandez, J. G., Lalejini, A., & Ofria, C. (2022).
 *      A suite of diagnostic metrics for characterizing selection schemes
 */

# pragma once

#include <cassert>
#include <cmath>
#include <numeric>
#include <iostream>

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"

// Gene bounds from Hernandez et al. suite
static constexpr double GENE_MIN = 0.0;
static constexpr double GENE_MAX = 100.0;

// Minimum and maximum mutation rates (per gene)
static constexpr double MUT_MIN = 1e-12; // if 0 can't be shifted
static constexpr double MUT_MAX = 1.0;

// Reflect a value back into [lo, hi] by "rebounding" off the walls.
// Also corrects cases where mutations overshoot by more than one interval width.
// For example: -0.7 -> 0.7, 100.7 -> 99.3
double ReflectIntoRange(double value, double lo=GENE_MIN, double hi=GENE_MAX) {
    assert(lo < hi && "Lower boundary must be less than the upper boundary!");
    const double width = hi - lo;
    double shifted_value = value - lo; // shift to [0, width]
    shifted_value = std::fmod(shifted_value, 2.0 * width);
    if (shifted_value < 0.0) shifted_value += 2.0 * width;
    if (shifted_value > width) shifted_value = 2.0 * width - shifted_value;
    return lo + shifted_value; // shift back to original range
}

// TODO: polymorphism with the other Organism class :)
class DiagOrganism {
    // using Translator = std::function<void(const emp::vector<double> & /*genome*/,
    //                                       emp::vector<double> & /*phenotype*/)>;
    using phenotype_t = emp::vector<double>;
    using genome_t = emp::vector<double>;
private:
    genome_t genome{};
    phenotype_t phenotype{};
    double fitness = 0.0; // obtained by aggregating phenotype values 
    double mut_rate = MUT_MIN;

    size_t start = 0; // starting index for fitness calculations

public:
    DiagOrganism() = default;
    DiagOrganism(const DiagOrganism &) = default;

    // We let external functions handle the conversion from genome to phenotype
    DiagOrganism(const genome_t & init_genome)
      : genome(init_genome), phenotype(init_genome.size(), 0.0) {}
   
    DiagOrganism(size_t genome_size)
      : genome(genome_size, 0.0), phenotype(genome_size, 0.0) {}

    DiagOrganism(size_t genome_size, emp::Random & random) 
      : genome(genome_size, 0.0), phenotype(genome_size, 0.0) {
        for (double & g : genome) {
            g = random.GetDouble(GENE_MIN, GENE_MAX);
        }
    }

    friend std::ostream & operator<<(std::ostream & os, const DiagOrganism & org) {
        os << "Genome="      << org.genome
           << ", Phenotype=" << org.phenotype
           << ", Fitness="   << org.fitness
           << ", Mutation="  << org.mut_rate;
        return os;
    }

    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }
    void SetGenome(const genome_t & g) { genome = g; }

    phenotype_t & GetPhenotype() { return phenotype; }
    const phenotype_t & GetPhenotype() const { return phenotype; }
    void SetPhenotype(const phenotype_t & p) { phenotype = p; }

    size_t GetGenomeSize() const { return genome.size(); }

    double GetFitness() const { return fitness; }
    void SetFitness(double f) { fitness = f; }

    double GetMutationRate() const { return mut_rate; }
    void SetMutationRate(double m) { mut_rate = std::clamp(m, MUT_MIN, MUT_MAX); }

    // // Translate genome into phenotype
    // void UpdatePhenotype() {
    //     if (!translator) phenotype = genome;
    //     else {
    //         // Let translator handle phenotype and genotype size differences
    //         translator(genome, phenotype); 
    //     }
    // }

    size_t GetStartIndex() const { return start; }
    void SetStartIndex(size_t s) {
        assert(s < phenotype.size());
        start = s;
    }

    // Aggregate phenotype by summing traits from start index
    void UpdateFitnessFromPhenotype() {
        assert(!phenotype.empty());
        assert(start < phenotype.size());
        fitness = std::accumulate(phenotype.begin() + start, 
                                  phenotype.end(), 0.0);
    }

    // Ensure all genes lie in [0, 100]
    void RepairGenome() {
        for (double & g : genome) g = ReflectIntoRange(g, GENE_MIN, GENE_MAX);
    }

    // Point mutations: for each gene, with prob mut_rate add N(0, 1) and rebound into [0,100].
    // Mutator mutations (optional): with prob pi, multiply mutation rate by 2^x where x ~ N(mean, sigma^2),
    // mean = b*sigma^2 (bias upward if b>0). Clamp to [MUT_MIN, MUT_MAX].
    DiagOrganism Mutate(emp::Random & random, 
                        bool const_mutation,
                        double sigma=0.138, 
                        double b=0.0, 
                        double pi=0.8) const {
        // Mutate genome using current mutation rate
        genome_t new_genome = genome;
        for (double & g : new_genome) {
            if (random.P(mut_rate)) {
                const double step = random.GetNormal(0.0, 1.0);
                g = ReflectIntoRange(g + step, GENE_MIN, GENE_MAX);
            }
        }

        // Inherit parent's mutation rate by default
        double mut_child = mut_rate;

        // OPTIONAL MUTATOR MUTATION
        if (!const_mutation && random.P(pi)) { 
            const double mean = b * sigma * sigma; // for upward bias, if b > 0

            // x = log2(mut_child / mut_parent) ~ N(mean, sigma^2)
            // For sigma=0.138, a 1-sigma draw is about a 1.10x increase (or 0.91x decrease)
            const double x = random.GetNormal(mean, sigma);
            mut_child *= std::pow(2.0, x); // multiplicative shift
            mut_child = std::clamp(mut_child, MUT_MIN, MUT_MAX);
        }

        DiagOrganism offspring(new_genome);
        offspring.SetFitness(0.0); // remove parent's fitness
        offspring.SetMutationRate(mut_child);
        offspring.SetStartIndex(start);
        // phenotype left as zeros by constructor, translator/evaluation will fill later
        return offspring;
    }
};