#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

#include "ExperimentRunner.hpp"
#include "MarkovChainAnalysis.hpp"
#include "QueueSimulator.hpp"

const double RELATIVE_ERROR_TOLERANCE = 0.05;
const double ABSOLUTE_ERROR_TOLERANCE = 0.01;

struct ValidationResult {
    std::string testName;
    double expected;
    double actual;
    double relativeError;
    bool passed;
};

double relativeError(double expected, double actual) {
    if (std::abs(expected) < 1e-10)
        return std::abs(actual - expected);
    return std::abs((actual - expected) / expected);
}

void printValidationResult(const ValidationResult& result) {
    std::cout << std::setw(40) << std::left << result.testName << " | " << std::setw(12)
              << std::fixed << std::setprecision(6) << result.expected << " | " << std::setw(12)
              << result.actual << " | " << std::setw(10) << std::setprecision(4)
              << (result.relativeError * 100) << "% | " << (result.passed ? "✓ PASS" : "✗ FAIL")
              << std::endl;
}

ValidationResult validateValue(const std::string& name, double expected, double actual) {
    ValidationResult result;
    result.testName = name;
    result.expected = expected;
    result.actual = actual;

    if (std::abs(expected) < 1e-10) {
        result.relativeError = std::abs(actual - expected);
        result.passed = result.relativeError < ABSOLUTE_ERROR_TOLERANCE;
    } else {
        result.relativeError = std::abs((actual - expected) / expected);
        result.passed = result.relativeError < RELATIVE_ERROR_TOLERANCE;
    }
    return result;
}

void runStableSystemAnalysis(double lambda, double mu, double simTime, int maxStates,
                             unsigned int seed) {
    double rho = lambda / mu;

    std::cout << "\n[Stable M/M/1 Simulation]\n";
    std::cout << "λ=" << lambda << ", μ=" << mu << ", ρ=" << rho << ", T=" << simTime << "\n";

    QueueSimulator sim(lambda, mu, maxStates, seed);
    sim.simulate(simTime);

    double simL = sim.getAverageQueueLength();
    double theoryL = rho / (1.0 - rho);

    std::cout << "Average queue length:\n";
    std::cout << "  Simulation : " << simL << "\n";
    std::cout << "  Theory     : " << theoryL << "\n";
    std::cout << "  Rel. Error : " << relativeError(theoryL, simL) * 100 << "%\n";

    auto simProbs = sim.getSteadyStateProbabilities();
    auto theoProbs = sim.getTheoreticalProbabilities();

    std::cout << "\nFirst 5 state probabilities:\n";
    std::cout << "n   Simulated    Theoretical\n";
    for (int i = 0; i < 5; ++i) {
        std::cout << i << "   " << simProbs[i] << "   " << theoProbs[i] << "\n";
    }

    sim.exportToCSV("data/mm1_stable_results.csv");
}

void runUnstableSystemAnalysis(double lambda, double mu, double simTime, int maxStates,
                               unsigned int seed) {
    double rho = lambda / mu;

    std::cout << "\n[Unstable M/M/1 Simulation]\n";
    std::cout << "λ=" << lambda << ", μ=" << mu << ", ρ=" << rho << ", T=" << simTime << "\n";

    QueueSimulator sim(lambda, mu, maxStates, seed);
    sim.simulate(std::min(simTime, 10000.0));

    std::cout << "Average queue length (finite horizon): " << sim.getAverageQueueLength() << "\n";

    sim.exportToCSV("data/mm1_unstable_results.csv");
}

void runSystemComparison(double mu, double simTime, int maxStates, unsigned int seed) {
    std::vector<double> rhos = {0.5, 0.7, 0.9, 1.0, 1.1};

    std::cout << "\n[Stability Comparison]\n";
    std::cout << "ρ    Simulated L    Theoretical L\n";

    for (double rho : rhos) {
        double lambda = rho * mu;

        QueueSimulator sim(lambda, mu, maxStates, seed);
        sim.simulate((rho < 1.0) ? simTime : 5000.0);

        std::cout << rho << "    " << sim.getAverageQueueLength() << "    ";

        if (rho < 1.0)
            std::cout << rho / (1.0 - rho) << "\n";
        else
            std::cout << "N/A\n";
    }
}

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
    double lambda = 0.8;
    double mu = 1.0;
    double simTime = 100000.0;
    int maxStates = 50;
    unsigned int seed = 42;

    if (argc > 1)
        lambda = std::stod(argv[1]);
    if (argc > 2)
        mu = std::stod(argv[2]);
    if (argc > 3)
        simTime = std::stod(argv[3]);

    // Create data directory for CSV output
    std::filesystem::create_directories("data");

    double rho = lambda / mu;

    if (rho < 1.0)
        runStableSystemAnalysis(lambda, mu, simTime, maxStates, seed);
    else
        runUnstableSystemAnalysis(lambda, mu, simTime, maxStates, seed);

    runSystemComparison(mu, simTime, maxStates, seed);

    ExperimentRunner experiments(mu, maxStates, simTime, seed);
    experiments.runVaryingLambdaExperiments();
    experiments.exportResults("data/mm1_varying_lambda.csv");

    return 0;
}
