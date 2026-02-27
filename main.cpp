#include <iostream>

#include "DiagWorld.hpp"

int main() {
    emp::Random rand {11};
    DiagWorld world;
    world.Run(rand);

}