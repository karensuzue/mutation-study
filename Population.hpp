/*
 * This class instantiates evolutionary runs.
 */

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cassert>
#include <iomanip>
#include <functional>

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"
#include "emp/bits/Bits.hpp"

#include "Organism.hpp"
#include "Select.hpp"
#include "Evaluate.hpp"

struct GenerationStats {
  size_t generation;
  double avg_f;
  double best_f;
  // double median_fitness;
  size_t best_id; // fittest

  double avg_mut; // mutation rate
  double highest_mut;
  size_t highest_mut_id;
};

class Population {
    using phenotype_t = emp::vector<double>;
    using genome_t = emp::vector<double>;
    using pop_t = emp::vector<Organism>;

    // Selects ONE parent index
    using selector_fn_t = std::function<size_t(const pop_t &, emp::Random &)>;
    // Translate genome to phenotype 
    using translator_fn_t = std::function<phenotype_t(const genome_t &)>;
    using evaluator_fn_t = std::function<double(const phenotype_t &)>;

private:
    pop_t organisms;

    size_t genome_size = 100;
    size_t max_generations = 50000;
    // size_t max_replicates = 20;
    size_t max_replicates = 100;
    size_t pop_size = 3600;

    double gene_min = -100.0;
    // double gene_min = 0.0;
    double gene_max = 100.0;

    // per-site mutation rate, applied to all orgs
    double init_mut_rate = 0.0;
    bool const_mutation_rate = false; // once toggled, this keeps mutation constant

    /* ------ SINGLE ENVIRONMENT ------ */
    bool valley_crossing = true; // toggles sawtooth transformation

    /* ------ CHANGING ENVIRONMENT ------ */
    // size_t change_env_step = 300; // change target genome every ? generations
    size_t change_env_step = 1000;
    size_t change_per_step = 1; // no. of target genes to modify at each environmental change
    genome_t target_genome;

    // TODO!!!!!!!!!!!!!!!!!!!!!!!!!
    bool rand_phenotype = false; // turns on stochastic phenotype expression

    size_t tour_size = 3;

    size_t generation = 0; // current generation
    size_t print_step = 100;   

    Select selector;
    Evaluate evaluator;

    selector_fn_t selector_fn;
    translator_fn_t translator_fn;
    evaluator_fn_t evaluator_fn;

    emp::String selector_name = "Tournament";
    emp::String translator_name = "ChangingEnv"; // genome to phenotype translator
    emp::String evaluator_name = "SquaredError"; // calculates fitness

    emp::vector<GenerationStats> history;

    // Taken from: https://github.com/jgh9094/ECJ-2023-Suite-Of-Diagnostic-Metrics-For-Characterizing-Selection-Schemes
    // Multi-valley crossing data
    // valley peaks for each floored integer gene value
    const phenotype_t peaks = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  8.0,  9.0,
                                9.0, 11.0, 11.0, 11.0, 14.0, 14.0, 14.0, 14.0, 18.0, 18.0,
                                18.0, 18.0, 18.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 29.0,
                                29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 36.0, 36.0, 36.0, 36.0,
                                36.0, 36.0, 36.0, 36.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0,
                                44.0, 44.0, 44.0, 53.0, 53.0, 53.0, 53.0, 53.0, 53.0, 53.0,
                                53.0, 53.0, 53.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0,
                                63.0, 63.0, 63.0, 63.0, 74.0, 74.0, 74.0, 74.0, 74.0, 74.0,
                                74.0, 74.0, 74.0, 74.0, 74.0, 74.0, 86.0, 86.0, 86.0, 86.0,
                                86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 99.0};
    // unique peaks
    const phenotype_t peaks_set = {8.0, 9.0, 11.0, 14.0, 18.0, 23.0, 29.0, 36.0, 44.0, 53.0, 63.0, 74.0, 86.0, 99.0};
    // where do the dips start?
    const double dips_start = 8.0;
    // where do dips end?
    const double dips_end = 99.9;

    emp::BitVector change_indices; // reusable buffer for ChangeTarget

public:
    Population() = default;
    Population(const Population &) = default;
    Population(Population &&) = default;

    Population & operator=(const Population &) = default;
    Population & operator=(Population &&) = default;

    Organism & operator[](size_t index) {
        return organisms[index];
    }
    const Organism & operator[](size_t index) const {
        return organisms[index];
    }

    size_t GetGenerations() const { return max_generations; }
    void SetGenerations(size_t g) { max_generations = g; }

    size_t GetReplicates() const { return max_replicates; }
    void SetReplicates(size_t r) { max_replicates = r; }

    size_t GetPopSize() const { return pop_size; }
    void SetPopSize(size_t ps) { pop_size = ps; }

    size_t GetGenomeSize() const { return genome_size; }
    void SetGenomeSize(size_t gs) { genome_size = gs; }

    void SetGeneMin(double g_min) { gene_min = g_min; }
    void SetGeneMax(double g_max) { gene_max = g_max; }

    void SetChangeEnvStep(size_t ces) { change_env_step = ces; }
    void SetChangePerStep(size_t cps) { change_per_step = cps; }
    genome_t GetTarget() const { return target_genome; }

    void SetInitMutation(double mu) { init_mut_rate = mu; }
    double GetInitMutation() const { return init_mut_rate; }
    // void ToggleConstantMutation() { const_mutation_rate = !const_mutation_rate; }
    void SetConstantMutation(bool cs) { const_mutation_rate = cs; }
    bool IsConstantMutation() const { return const_mutation_rate; }

    // void ToggleValleyCrossing() { valley_crossing = !valley_crossing; }
    void SetValleyCrossing(bool vc) { valley_crossing = vc; }
    bool IsValleyCrossing() const { return valley_crossing; }

    // TODO
    void SetStochasticPhenotype(bool sp) { rand_phenotype = sp; }
    bool IsStochasticPhenotype() const { return rand_phenotype; }

    void SetSelector(const emp::String & selector) { selector_name = selector; }
    void SetTranslator(const emp::String & translator) { translator_name = translator; }
    void SetEvaluator(const emp::String & evaluator) { evaluator_name = evaluator; } 

    friend std::ostream & operator<<(std::ostream & os, const Population & pop) {
        assert(pop.pop_size == pop.organisms.size() && 
               "pop_size does not match pop.organisms.size().");
        for (size_t i = 0; i < pop.pop_size; ++i) {
            os << pop[i] << '\n';
        }
        os << std::endl;
        return os;
    }

    // Copy the given 'org' a 'pop_size' number of times
    void InitializeUniform(const Organism & org) {
        assert(org.GetGenome().size() == genome_size && 
               "Input organism genome size does not match the configured value.");
        organisms.clear(); // justtttt in case
        organisms.reserve(pop_size);
        for (size_t i = 0; i < pop_size; ++i) {
            organisms.push_back(org);
        }
    }

    // Copy the given genome a 'pop_size' number of times
    void InitializeUniform(const genome_t & genome) {
        assert(genome.size() == genome_size && 
               "Input genome size does not match the configured value.");
        organisms.clear();
        organisms.reserve(pop_size);
        for (size_t i = 0; i < pop_size; ++i) {
            organisms.emplace_back(genome, gene_min, gene_max);
        }
    }

    // Create a random genome, make 'pop_size' copies of it
    void InitializeUniform(emp::Random & random) {
        organisms.clear();
        organisms.reserve(pop_size);
        Organism org(genome_size, random, gene_min, gene_max);
        for (size_t i = 0; i < pop_size; ++i) {
            organisms.push_back(org);
        }
    }

    // Initialize the population with all-zero genomes
    void InitializeUniform() {
        organisms.clear();
        organisms.reserve(pop_size);
        for (size_t i = 0; i < pop_size; ++i) {
            organisms.emplace_back(genome_size, gene_min, gene_max);
        }
    }

    // Add a 'pop_size' number of randomized organisms with the configured genome size
    void InitializeRandom(emp::Random & random) {
        organisms.clear();
        organisms.reserve(pop_size);
        for (size_t i = 0; i < pop_size; ++i) {
            organisms.emplace_back(genome_size, random, gene_min, gene_max);
        }        
    }
    
    // Set all organisms to the same mutation rate
    void SetPopulationMutation(double mutation_rate) {
        assert(!organisms.empty() && "The population is empty.");
        for (Organism & org : organisms) {
            org.SetMutationRate(mutation_rate);
        }
    }

    void ConfigureSelector() {
        if (selector_name == "Tournament") {
            selector_fn = [this](const pop_t & pop,
                        emp::Random& random) {
                            return selector.Tournament(pop, tour_size, random);
                        };
        }
        else {
            throw std::runtime_error("Unknown Selector: " + selector_name);
        }
    }

    // This sets up the translator that maps genomes to phenotypes
    void ConfigureTranslator() {
        if (translator_name == "ExploitationRate") {
            translator_fn = [this](const genome_t & g) {
                return evaluator.ExploitationRate(g);
            };
        }
        else if (translator_name == "ChangingEnv") { // impostor! ;)
            translator_fn = [this](const genome_t & g) {
                return evaluator.ExploitationRate(g);
            };
        }
        else {
            throw std::runtime_error("Unknown Translator: " + translator_name);
        }
    }

    void ConfigureEvaluator() {
        if (evaluator_name == "Aggregate") {
            evaluator_fn = [this](const phenotype_t & p) {
                return evaluator.AggregateFitness(p);
            };
        }
        else if (evaluator_name == "SquaredError") {
            evaluator_fn = [this](const phenotype_t & p) {
                return evaluator.SquaredErrorFitness(p, target_genome);
            };
        }
        else {
            throw std::runtime_error("Unknown Evaluator: " + evaluator_name);
        }
    }
    
    // Initialize target genome sequence if `changing_env` is on
    void InitializeTarget(emp::Random & random) {
        target_genome.resize(genome_size);
        for (double & gene : target_genome) {
            gene = random.GetDouble(gene_min, gene_max);
        }
        change_indices.Resize(genome_size);
    }

    // Changes a random gene in the target genome
    void ChangeTarget(emp::Random & random) {
        assert(target_genome.size() == genome_size && "target_genome size does not match genome_size.");
        change_indices.Clear();
        change_indices.ChooseRandom(random, change_per_step);
        for (size_t idx : change_indices) { // Iterates only over bits set to 1
            std::cout << idx << ", ";
            target_genome[idx] = random.GetDouble(gene_min, gene_max);
        }
    }

    // This function computes the fitness of the whole population
    void EvaluateFitness() {
        for (Organism & org : organisms) {
            const genome_t & g = org.GetGenome();

            phenotype_t p = translator_fn(g);

            if (valley_crossing) {
                p = evaluator.MultiValleyCrossing(p, peaks, dips_start, dips_end);
            }

            org.SetPhenotype(p);
            double fitness = evaluator_fn(p);
            org.SetFitness(fitness);
        }
    }

    // In each generation, produce a 'pop_size' number of offspring
    void RunOneGeneration(emp::Random & random) {
        // In case something goes wrong in the process...
        assert(pop_size == organisms.size() && 
               "pop_size does not match pop.organisms.size().");

        pop_t next_pop;
        next_pop.reserve(pop_size);
        for (size_t i = 0; i < pop_size; ++i) {
            const size_t parent_idx = selector_fn(organisms, random);
            next_pop.push_back(organisms[parent_idx].Mutate(random, const_mutation_rate));
        }
        organisms.swap(next_pop);
    }

    void Run(emp::Random & random) {
        history.clear();
        history.reserve(max_generations + 1);

        if (!selector_fn) ConfigureSelector();
        if (!translator_fn) ConfigureTranslator();
        if (!evaluator_fn) ConfigureEvaluator();

        // for now
        if (translator_name == "ChangingEnv") valley_crossing = false; 

        InitializeUniform(random);
        InitializeTarget(random);
        SetPopulationMutation(init_mut_rate);

        for (generation = 0; generation < max_generations; ++generation) {
            if (generation > 0 && generation % change_env_step == 0 && translator_name == "ChangingEnv") {
                ChangeTarget(random);
            }
            EvaluateFitness();
            RecordGeneration(generation);
            // if (generation % print_step == 0) PrintStats(generation);
            RunOneGeneration(random);
        }
        EvaluateFitness();
        RecordGeneration(max_generations);
        // PrintStats(max_generations);
    }

    void MultiRun(const std::string & prefix, size_t num_replicates = 0) {
        if (num_replicates == 0) num_replicates = max_replicates;
        for (size_t replicate = 0; replicate < num_replicates; ++replicate) {
            emp::Random random(replicate + 1);
            Run(random);
            ExportHistory("history_" + prefix + "_" + std::to_string(replicate));
            std::cout << "Replicate " << replicate << " completed!\n";
        }
    }

    void RecordGeneration(size_t gen) {
        assert(!organisms.empty());

        double avg_f = 0.0;
        // double median_f = 0.0;
        double best_f = organisms[0].GetFitness();
        size_t best_id = 0;

        double avg_mut = 0.0;
        double highest_mut = organisms[0].GetMutationRate();
        size_t highest_mut_id = 0; // they start out the same

        assert(pop_size == organisms.size() && 
               "pop_size does not match pop.organisms.size().");
        for (size_t i = 0; i < pop_size; ++i) {
            const Organism & org = organisms[i];
            
            double org_f = org.GetFitness();
            avg_f += org_f;
            if (org_f > best_f) {
                best_f = org_f;
                best_id = i;
            }

            double org_mut = org.GetMutationRate();
            avg_mut += org_mut;
            if (org_mut > highest_mut) {
                highest_mut = org_mut;
                highest_mut_id = i;
            }
        }
        avg_f /= pop_size;
        avg_mut /= pop_size;
        // assert(best_id > -1 && "Best ID could not be logged.");

        history.emplace_back(gen, avg_f, best_f, best_id, avg_mut, highest_mut, highest_mut_id);
    }

    void PrintStats(size_t gen) const {
        std::cout << "Generation=" << gen << ", "
                  << "AvgFit=" << history[gen].avg_f << ", "
                  << "BestFit=" << history[gen].best_f << ", "
                  << "BestID=" << history[gen].best_id << ", "
                  << "AvgMut=" << history[gen].avg_mut << ", "
                  << "HighMut=" << history[gen].highest_mut << ", "
                  << "HighMutID=" << history[gen].highest_mut_id << "\n";
    }


    void ExportHistory(const std::string & filename="history") const {
        assert(history.size() > 0 && "Cannot export history because it is empty.");

        std::ofstream fitness_file(filename + "_fitness.csv");
        std::ofstream mutation_file(filename + "_mutation.csv");

        if (fitness_file.is_open()) {
            fitness_file << "Generation,Best_F,Mean_F,Fittest_ID\n";
            fitness_file << std::fixed << std::setprecision(3);

            for (const GenerationStats & gen : history) {
                fitness_file << gen.generation << ","
                             << gen.best_f << ","
                             << gen.avg_f << ","
                             << gen.best_id << "\n";
            }
        }

        if (mutation_file.is_open()) {
            // mu: per-site mutation rate
            // U: expected number of mutations per organism (genome-wide)
            mutation_file << "Generation,Highest_mu,Mean_mu,Highest_U,Mean_U,Highest_ID\n";
            mutation_file << std::fixed << std::setprecision(6);

            for (const GenerationStats & gen : history) {
                mutation_file << gen.generation << ","
                              << gen.highest_mut << ","
                              << gen.avg_mut << ","
                              << gen.highest_mut * genome_size << ","
                              << gen.avg_mut * genome_size << ","
                              << gen.highest_mut_id << "\n";
            }
        }
    }
};
