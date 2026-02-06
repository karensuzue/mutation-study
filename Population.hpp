#pragma once

#include <iostream>
#include <string>
#include <memory>
#include <stdexcept>

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"
#include "emp/config/SettingsManager.hpp"

#include "Organism.hpp"
#include "Select.hpp"
#include "Fitness.hpp"


class Population {
private:
    emp::vector<Organism> organisms;

    size_t genome_size = 10;
    size_t max_generations = 15000;
    size_t max_replicates = 20;
    size_t pop_size = 3600;
    
    double mutation_rate = 0.0; // per gene

    std::string selector_name = "Tournament";
    std::string fitness_name = "OneMax";

    size_t tour_size = 2;
    size_t intmatch_target = 5;

    size_t generation = 0;
    size_t print_step = 100;   

    emp::Random random{5};

    std::unique_ptr<Selector> selector;
    std::unique_ptr<Fitness> fitness_function;

public:
    // void SetupConfig(emp::SettingsManager & settings) {
    //     settings.AddSetting("genome_size", genome_size, "How long is genome?", "org_size");
    //     settings.AddSetting("max_generations", max_generations, "How many generations to run?", "max_gen");
    //     settings.AddSetting("max_replicates", max_replicates, "How many replicates should be performed?", "rep");
    //     settings.AddSetting("pop_size", pop_size, "Size of the population", "pop_size");
    //     settings.AddSetting("mutation_rates", mutation_rates, "List of mutation rates to test", "mut_rate");
    //     settings.AddSetting("selector", selector_name, "Which selection scheme to choose? [Tournament, Truncation (not available)]", "select");
    //     settings.AddSetting("tournament_size", tour_size, "Tournament size (if selector=tournament)", "tour_size");
    //     settings.AddSetting("fitness", fitness_name, "Which fitness function to choose? [OneMax, FlipOnes]", "fitness");
    //     settings.AddSetting("print_step", print_step, "How many generations between printing outputs?", "print");

    //     // if (selector_name == "Tournament") {
    //     //     selector = ;

    //     // }
    //     // if (fitness_name == "OneMax") {
    //     //     fitness_fn = ;
    //     // }
    // }


    Organism & operator[](unsigned int index) {
        return organisms[index];
    }

    const Organism & operator[](unsigned int index) const {
        return organisms[index];
    }

    size_t GetSize() const {
        return organisms.size(); // could return pop_size here but...
    }

    // for debug mostly
    void PrintPopulation(std::ostream & os) const {
        for (size_t i = 0; i < organisms.size(); ++i) {
            os << organisms[i] << ", ";
        }
        os << std::endl;
    }

    // Add a 'pop_size' number of isogenic organisms
    void InitializeUniform() {
        organisms.clear(); // justtttt in case
        for (size_t i = 0; i < pop_size; ++i) {
            organisms.emplace_back(genome_size);
        }
    }

    // Add a 'pop_size' number of randomized organisms
    void InitializeRandom() {
        organisms.clear();
        for (size_t i = 0; i < pop_size; ++i) {
            emp::BitVector genome(genome_size, random);
            organisms.emplace_back(genome);
        }        
    }   

    void ConfigureSelector() {
        if (selector_name == "Tournament") {
            selector = std::make_unique<TournamentSelector>(tour_size);
        }
        else {
            throw std::runtime_error("Unknown Selector: " + selector_name);
        }
    }

    void ConfigureFitnessFunction() {
        if (fitness_name == "OneMax") {
            fitness_function = std::make_unique<OneMaxFitness>();
        }
        else if (fitness_name == "IntMatch") {
            fitness_function = std::make_unique<IntMatchFitness>(intmatch_target, genome_size);
        }
        else {
            throw std::runtime_error("Unknown fitness function: " + fitness_name);
        }
    }

    void CalculateFitness() {
        for (size_t i = 0; i < pop_size; ++i) {
            double fitness = fitness_function->Evaluate(organisms[i]);
            organisms[i].SetFitness(fitness);
        }
    }

    // In each generation, produce a 'pop_size' number of offspring
    void RunOneGeneration() {
        if (!selector) ConfigureSelector();
        if (!fitness_function) ConfigureFitnessFunction();

        emp::vector<Organism> next_pop;

        for (size_t i = 0; i < pop_size; ++i) {
            const size_t parent_idx = selector->Select(organisms, random);
            Organism child = organisms[parent_idx].Mutate(random, mutation_rate);
            next_pop.push_back(child);
        }

        organisms = next_pop;
    }

    void Run() {
        for (generation = 0; generation < max_generations; ++generation) {
            CalculateFitness();
            RunOneGeneration();
            if (generation % print_step == 0) { 
                // PrintStats();
            }
        }
    }

    // void PrintStats() {

    // }
};
