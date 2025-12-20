#include <iomanip>
#include <iostream>
#include <string>

#include "header/QueueSimulator.hpp"

int main(int argc, char* argv[]) {
    double lambda = 0.8;
    double mu = 1.0;
    double simTime = 100000.0;
    int maxStates = 50;
    unsigned int seed = 42;

    // parameter parsing
    if (argc > 1)
        lambda = std::stod(argv[1]);
    if (argc > 2)
        mu = std::stod(argv[2]);
    if (argc > 3)
        simTime = std::stod(argv[3]);
    if (argc > 4)
        maxStates = std::stoi(argv[4]);
    if (argc > 5)
        seed = std::stoi(argv[5]);

    // params check
    std::cout << "Parameters: λ=" << lambda << ", μ=" << mu << ", simulation_time=" << simTime
              << std::endl;

    QueueSimulator simulator(lambda, mu, maxStates, seed);
}
