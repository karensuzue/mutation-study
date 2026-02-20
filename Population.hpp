#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <stdexcept>
#include <cassert>
#include <iomanip>

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"

#include "Organism.hpp"
#include "Select.hpp"
#include "Fitness.hpp"

struct GenerationStats {
  int generation;
  double avg_f;
  double best_f;
  // double median_fitness;
  size_t best_id; // fittest

  double avg_mut; // mutation rate
  double highest_mut;
  size_t highest_mut_id;
};


class Population {
private:
    emp::vector<Organism> organisms;

    size_t genome_size = 100;
    size_t max_generations = 1000;
    size_t max_replicates = 20;
    size_t pop_size = 3600;

    std::string selector_name = "Tournament";
    // std::string fitness_name = "OneMax";
    // std::string fitness_name = "AggregateMatch";
    // std::string fitness_name = "BitwiseMatch";
    std::string fitness_name = "KDeceptiveTrap";
    // std::string fitness_name = "NKLandscape";

    size_t tour_size = 3;
    size_t intmatch_target = 784; // for IntMatch problem
    size_t num_blocks = 20; // for KDeceptiveTrap problem

    size_t generation = 0; // current generation
    size_t print_step = 100;   

    // int seed = 11;
    // emp::Random random{seed};

    std::shared_ptr<Selector> selector;
    std::shared_ptr<Fitness> fitness_function;

    emp::vector<GenerationStats> history;

    bool const_mutation_rate = false; // once toggled, this keeps mutation constant

public:
    Population() = default;
    Population(const Population &) = default;
    Population(Population &&) = default;

    Population & operator=(const Population &) = default;
    Population & operator=(Population &&) = default;

    Organism & operator[](unsigned int index) {
        return organisms[index];
    }
    const Organism & operator[](unsigned int index) const {
        return organisms[index];
    }

    size_t GetSize() const { return organisms.size(); }

    void ToggleConstantMutation() { const_mutation_rate = !const_mutation_rate; }
    bool IsConstantMutation() { return const_mutation_rate; }

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
        for (size_t i = 0; i < pop_size; ++i) {
            organisms.push_back(org);
        }
    }

    // Copy the given genome a 'pop_size' number of times
    void InitializeUniform(const emp::BitVector & genome) {
        assert(genome.size() == genome_size && 
               "Input genome size does not match the configured value.");
        organisms.clear();
        for (size_t i = 0; i < pop_size; ++i) {
            organisms.emplace_back(genome);
        }
    }

    // Add a 'pop_size' number of isogenic organisms with the configured genome size
    void InitializeUniform(emp::Random & random) {
        organisms.clear(); 
        emp::BitVector genome(genome_size, random);
        for (size_t i = 0; i < pop_size; ++i) {
            organisms.emplace_back(genome);
        }
    }

    // Add a 'pop_size' number of randomized organisms with the configured genome size
    void InitializeRandom(emp::Random & random) {
        organisms.clear();
        for (size_t i = 0; i < pop_size; ++i) {
            emp::BitVector genome(genome_size, random);
            organisms.emplace_back(genome);
        }        
    }
    
    // Set all organisms to the same mutation rate
    void SetPopulationMutation(double mutation_rate) {
        assert(organisms.size() > 0 && "The population is empty.");
        for (Organism & org : organisms) {
            org.SetMutationRate(mutation_rate);
        }
    }

    void ConfigureSelector() {
        if (selector_name == "Tournament") {
            selector = std::make_shared<TournamentSelector>(tour_size);
        }
        else {
            throw std::runtime_error("Unknown Selector: " + selector_name);
        }
    }

    void ConfigureFitnessFunction() {
        if (fitness_name == "OneMax") {
            fitness_function = std::make_shared<OneMaxFitness>();
        }
        else if (fitness_name == "AggregateMatch") {
            fitness_function = std::make_shared<AggregateMatchFitness>(intmatch_target, genome_size);
        }
        else if (fitness_name == "BitwiseMatch") {
            fitness_function = std::make_shared<BitwiseMatchFitness>(intmatch_target, genome_size);
        }
        else if (fitness_name == "KDeceptiveTrap") {
            fitness_function = std::make_shared<KDeceptiveTrapFitness>(genome_size, num_blocks);
        }
        else {
            throw std::runtime_error("Unknown fitness function: " + fitness_name);
        }
    }

    void EvaluateFitness() {
        if (!fitness_function) ConfigureFitnessFunction();
        for (Organism & org : organisms) {
            double fitness = fitness_function->Evaluate(org);
            org.SetFitness(fitness);
        }
    }

    // In each generation, produce a 'pop_size' number of offspring
    void RunOneGeneration(emp::Random & random) {
        assert(pop_size == organisms.size() && 
               "pop_size does not match pop.organisms.size().");

        if (!selector) ConfigureSelector();
        if (!fitness_function) ConfigureFitnessFunction();

        emp::vector<Organism> next_pop(pop_size);

        for (size_t i = 0; i < pop_size; ++i) {
            const size_t parent_idx = selector->Select(organisms, random);
            next_pop[i] = organisms[parent_idx].Mutate(random, const_mutation_rate);
        }

        organisms.swap(next_pop);
    }

    void Run(emp::Random & random) {
        InitializeUniform(random);
        SetPopulationMutation(0.05);

        for (generation = 0; generation < max_generations; ++generation) {
            EvaluateFitness();
            RecordGeneration(generation);
            if (generation % print_step == 0) { 
                PrintStats(generation);
            }
            RunOneGeneration(random);
        }
    }

    void MultiRun(size_t num_replicates = 0) {
        if (num_replicates == 0) num_replicates = max_replicates;
        for (size_t replicate = 0; replicate < num_replicates; ++replicate) {
            history.clear();
            emp::Random random(replicate + 1);
            Run(random);
            ExportHistory("history_" + std::to_string(replicate));
        }
    }

    void RecordGeneration(size_t gen) {
        double avg_f = 0.0;
        // double median_f = 0.0;
        double best_f = 0.0;
        int best_id = -1;

        double avg_mut = 0.0;
        double highest_mut = 0.0;
        int highest_mut_id = 0; // they start out the same

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
        assert(best_id > -1 && "Best ID could not be logged.");

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
            mutation_file << "Generation,Best_mu,Mean_mu,Best_U,Mean_U,Highest_ID\n";
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
