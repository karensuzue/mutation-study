#pragma once

#include <cassert>
#include <cstdint>

#include "Organism.hpp"

struct Fitness {
    virtual ~Fitness() = default;
    virtual double Evaluate(const Organism &) const = 0; // higher is better!
};

// Try to set every gene to 1
struct OneMaxFitness : Fitness {
    double Evaluate(const Organism & org) const override {
        return static_cast<double>(org.GetGenome().CountOnes());
    }
};

// Given a target unsigned integer, make sure genes sum up to it
struct AggregateMatchFitness : Fitness {
    size_t target;
    size_t genome_len;

    AggregateMatchFitness(size_t target, size_t genome_len) 
        : target(target), genome_len(genome_len) {
            assert(genome_len > 0 && "Genome length can not be 0!");
            assert(target <= genome_len && "Target cannot exceed genome length.");
    }

    double Evaluate(const Organism & org) const override {
        const emp::BitVector & genome = org.GetGenome();
        assert(genome.size() == genome_len && "Genome length mismatch!");
        const size_t ones = org.GetGenome().CountOnes();
        const size_t dist = (ones > target) ? (ones - target) : (target - ones);
        // normalize (best = 1, worst = 0)
        return 1.0 - static_cast<double>(dist) / static_cast<double>(genome_len);
    }
};

// Given a target unsigned integer, ensure that the binary genome as a whole matches it
struct BitwiseMatchFitness : Fitness {
    uint64_t target;
    size_t genome_len;

    BitwiseMatchFitness(uint64_t target, size_t genome_len)
        : target(target), genome_len(genome_len) {
            assert(genome_len > 0 && "Genome length can not be 0!");
            assert (genome_len <= 63 && "BitwiseMatchFitness only supports genome_len <= 63.");
            assert(target < (1ULL << genome_len) && // target < 2^(genome_len)
               "Target value exceeds maximum representable integer for this genome length.");
    }
    
    double Evaluate(const Organism & org) const override {
        const emp::BitVector & genome = org.GetGenome();
        assert(genome.size() == genome_len && "Genome length mismatch!");

        // Convert bit genome into integer
        uint64_t value = 0;
        for (size_t i = 0; i < genome_len; ++i) {
            value |= (uint64_t(genome[i]) << i);
        }

        const uint64_t dist = (value > target) ? (value - target) : (target - value);
        const uint64_t max_value = (1ULL << genome_len) - 1ULL;
        return 1.0 - (double(dist) / double(max_value));
    }
};

// k-bit deceptive trap fitness
// Genome is split into blocks of size 'k'
// For each block with 'u' 1s:
// score = k                     if u == k (global optimum for that block is all ones)
//       = (k - 1) - (a * u)     otherwise (deceptive slope toward all 0s)
// where 'a' is a smoothing parameter to reduce deceptiveness
// So each block has:
//   - local optimum at u=0 with score k-1 (all zeros)
//   - global optimum at u=k with score k (all ones)
// Total score is the sum of scores obtained in each block
struct KDeceptiveTrapFitness : Fitness {
    size_t genome_len;
    size_t num_blocks;
    size_t k;
    double smooth;

    KDeceptiveTrapFitness(size_t genome_len, size_t num_blocks=1, double smooth=0.25)
        : genome_len(genome_len), num_blocks(num_blocks), k(genome_len / num_blocks), smooth(smooth) {
        assert(genome_len > 0 && "Genome length can not be 0!");
        assert(num_blocks > 0 && "num_blocks can not be 0!");
        assert(genome_len % num_blocks == 0 && "Genome length must be divisible by the number of blocks");
        assert(k >= 2 && "k must be at least 2.");
        assert(smooth >= 0.0 && smooth <= 1.0 && "smooth should be in [0, 1].");
    }

    double Evaluate(const Organism & org) const override {
        const emp::BitVector & genome = org.GetGenome();
        assert(genome.size() == genome_len && "Genome length mismatch!");

        double total = 0.0;
        for (size_t block = 0; block < num_blocks; ++block) {
            const size_t start = block * k;
            const size_t count = genome.CountOnes(start, start + k); // [start, end) ?

            const double block_score = (count == k) ? 
                        static_cast<double>(k) : static_cast<double>((k - 1) - (count * smooth));
            total += block_score;
        }

        // Normalize total score [0, 1]
        const double max_total = static_cast<double>(num_blocks * k);
        return total / max_total;
    }
};


