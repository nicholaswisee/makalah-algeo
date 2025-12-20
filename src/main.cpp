#include <iomanip>
#include <iostream>
#include <string>

#include "ExperimentRunner.hpp"
#include "MarkovChainAnalysis.hpp"
#include "QueueSimulator.hpp"

void printUsage(const char* programName) {
    std::cout << "\nUsage: " << programName << " [lambda] [mu] [sim_time] [max_states] [seed]\n"
              << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  lambda     : Arrival rate (default: 0.8)" << std::endl;
    std::cout << "  mu         : Service rate (default: 1.0)" << std::endl;
    std::cout << "  sim_time   : Simulation duration (default: 100000.0)" << std::endl;
    std::cout << "  max_states : Maximum states to track (default: 50)" << std::endl;
    std::cout << "  seed       : Random seed (default: 42)" << std::endl;
    std::cout << "\nExample: " << programName << " 0.8 1.0 100000 50 42" << std::endl;
}

int main(int argc, char* argv[]) {
    // Default parameters
    double lambda = 0.8;
    double mu = 1.0;
    double simTime = 100000.0;
    int maxStates = 50;
    unsigned int seed = 42;

    // Parse command-line arguments
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

    std::cout << "=== M/M/1 Queue Analysis System ===" << std::endl;
    std::cout << "Parameters: λ=" << lambda << ", μ=" << mu << ", simulation_time=" << simTime
              << std::endl;

    // Run discrete-event simulation
    std::cout << "\n[1] Running discrete-event simulation..." << std::endl;
    QueueSimulator simulator(lambda, mu, maxStates, seed);
    simulator.simulate(simTime);
    simulator.printSummary();

    // Export simulation results
    simulator.exportToCSV("mm1_simulation_results.csv");

    // Perform Markov Chain CTMC analysis (only for stable systems)
    if (lambda < mu) {
        std::cout << "\n[2] Performing Continuous-Time Markov Chain (CTMC) analysis..."
                  << std::endl;
        MarkovChainAnalysis analyzer(maxStates, lambda, mu);

        // Print generator matrix for verification
        analyzer.printGeneratorMatrix();

        std::cout << "\nComputing steady-state distribution by solving πQ = 0..." << std::endl;
        auto ctmcSteadyState = analyzer.computeSteadyState();

        std::cout << "\n--- Comparison: Simulation vs CTMC Generator Matrix Method ---"
                  << std::endl;
        std::cout << "State | Simulation | CTMC (πQ=0) | Theoretical" << std::endl;
        std::cout << "------|------------|-------------|------------" << std::endl;

        auto simProbs = simulator.getSteadyStateProbabilities();
        auto theoProbs = simulator.getTheoreticalProbabilities();

        for (int i = 0; i < std::min(10, maxStates); ++i) {
            std::cout << std::setw(5) << i << " | " << std::setw(10) << std::fixed
                      << std::setprecision(6) << simProbs[i] << " | " << std::setw(11)
                      << ctmcSteadyState[i] << " | " << std::setw(10) << theoProbs[i] << std::endl;
        }
    }

    // Run multiple experiments with varying arrival rates
    std::cout << "\n[3] Running experiments with varying arrival rates..." << std::endl;
    ExperimentRunner experiments(mu, maxStates, simTime, seed);
    experiments.runVaryingLambdaExperiments();
    experiments.exportResults("mm1_varying_lambda.csv");

    std::cout << "\n=== Simulation Complete ===" << std::endl;

    return 0;
}
