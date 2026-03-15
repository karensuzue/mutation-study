#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "DiagWorld.hpp"

struct RunConfig {
    size_t seed = 1;
    double start_U = 1e-3; // GENOME-WIDE mutation rate
    bool const_mut = true;
    bool valley_cross = true;
    bool rand_phenotype = false;
    size_t gens = 1000;
    // size_t reps = 1;
    size_t pop_size = 100;
    size_t genome_size = 10;
};


void PrintUsage() {
    std::cout << "Options:\n"
        << "  --seed <size_t>           Seed\n"
        << "  --start_U <double>        Starting genome-wide mutation rate\n"
        << "  --const_mut <0|1>         1 = constant mutation, 0 = evolving mutation\n"
        << "  --valley_cross <0|1>      1 = apply sawtooth transformation, 0 = don't apply\n"
        << "  --rand_phenotype <0|1>    1 = stochastic phenotype, 0 = non-stochastic phenotype\n"
        << "  --gens <size_t>           Number of generations\n"
        // << "  --reps <size_t>           Number of replicates\n"
        << "  --pop_size <size_t>       Population size\n"
        << "  --genome_size <size_t>    Number of genes per organism\n"
        << "  --help                    Show this message\n";
}


RunConfig ParseArgs(int argc, char * argv[]) {
    RunConfig cfg;

    for (int i = 1; i < argc; ++i) {
        emp::String arg = argv[i];

        // Checks for a value after flag
        auto require_value = [&](const std::string & name) {
            if (i + 1 >= argc) { 
                emp::notify::Error("Missing value for ", name, "."); 
                std::exit(EXIT_FAILURE);
            }
            emp::String next = argv[i + 1];
            if (next[0] == '-') { 
                emp::notify::Error("Missing value for ", name, ", got option '", next, "' instead."); 
                std::exit(EXIT_FAILURE);
            }
        };

        if (arg == "--help") {
            PrintUsage();
            std::exit(EXIT_SUCCESS);
        }
        else if (arg == "--seed") {
            require_value(arg);
            auto val = std::stoull(argv[++i]);
            cfg.seed = static_cast<size_t>(val);
        }
        else if (arg == "--start_U") {
            require_value(arg);
            cfg.start_U = std::stod(argv[++i]);
        }
        else if (arg == "--const_mut") {
            require_value(arg);
            int val = std::stoi(argv[++i]);
            if (val != 0 && val != 1) {
                emp::notify::Error("--const_mut must be 0 or 1.");
                std::exit(EXIT_FAILURE);
            }
            cfg.const_mut = static_cast<bool>(val);
        }
        else if (arg == "--valley_cross") {
            require_value(arg);
            int val = std::stoi(argv[++i]);
            if (val != 0 && val != 1) {
                emp::notify::Error("--valley_cross must be 0 or 1.");
                std::exit(EXIT_FAILURE);
            }
            cfg.valley_cross = static_cast<bool>(val);
        }
        else if (arg == "--rand_phenotype") {
            require_value(arg);
            int val = std::stoi(argv[++i]);
            if (val != 0 && val != 1) {
                emp::notify::Error("--rand_phenotype must be 0 or 1.");
                std::exit(EXIT_FAILURE);
            }
            cfg.rand_phenotype = static_cast<bool>(val);
        }
        else if (arg == "--gens") {
            require_value(arg);
            auto val = std::stoull(argv[++i]);
            if (val <= 0) {
                emp::notify::Error("--gens must be positive.");
                std::exit(EXIT_FAILURE);
            }
            cfg.gens = static_cast<size_t>(val);
        }
        // else if (arg == "--reps") {
        //     require_value(arg);
        //     auto val = std::stoull(argv[++i]);
        //     if (val <= 0) {
        //         emp::notify::Error("--reps must be positive.");
        //         std::exit(EXIT_FAILURE);
        //     }
        //     cfg.reps = static_cast<size_t>(val);
        // }
        else if (arg == "--pop_size") {
            require_value(arg);
            auto val = std::stoull(argv[++i]);
            if (val <= 0) {
                emp::notify::Error("--pop_size must be positive.");
                std::exit(EXIT_FAILURE);
            }
            cfg.pop_size = static_cast<size_t>(val);
        }
        else if (arg == "--genome_size") {
            require_value(arg);
            auto val = std::stoull(argv[++i]);
            if (val <= 0) {
                emp::notify::Error("--genome_size must be positive.");
                std::exit(EXIT_FAILURE);
            }
            cfg.genome_size = static_cast<size_t>(val);
        }
        else {
            emp::notify::Error("Unknown argument: ", arg);
            std::exit(EXIT_FAILURE);
        }
    }
    return cfg;
}

int main(int argc, char * argv[]) { 
    RunConfig cfg = ParseArgs(argc, argv);

    std::cout << "seed = " << cfg.seed << "\n";
    std::cout << "start_U = " << cfg.start_U << "\n";
    std::cout << "const_mut = " << cfg.const_mut << "\n";
    std::cout << "valley_cross = " << cfg.valley_cross << "\n";
    std::cout << "rand_phenotype = " << cfg.rand_phenotype << "\n";
    std::cout << "gens = " << cfg.gens << "\n";
    // std::cout << "reps = " << cfg.reps << "\n";
    std::cout << "pop_size = " << cfg.pop_size << "\n";
    std::cout << "genome_size = " << cfg.genome_size << "\n";
    
    DiagWorld pop;
    
    pop.SetConstantMutation(cfg.const_mut);
    pop.SetValleyCrossing(cfg.valley_cross);
    pop.SetStochasticPhenotype(cfg.rand_phenotype);
    pop.SetGenerations(cfg.gens);
    // pop.SetReplicates(cfg.reps);
    pop.SetPopSize(cfg.pop_size);
    pop.SetGenomeSize(cfg.genome_size);

    const double mu = cfg.start_U / static_cast<double>(cfg.genome_size);
    pop.SetInitMutation(mu);

    // Formatting
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(4) << cfg.start_U;
    const std::string tag = oss.str();

    std::cout << "Running U=" << tag << " (mu=" 
            << std::scientific << std::setprecision(4) << mu << ")\n";

    // pop.MultiRun(tag);
    emp::Random random(cfg.seed);
    pop.Run(random);
    pop.ExportHistory("history_" + tag + "_" + std::to_string(cfg.seed));
}