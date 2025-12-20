#include "../header/MarkovChainAnalysis.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

// constructor
MarkovChainAnalyzer::MarkovChainAnalyzer(int states, double arrivalRate, double serviceRate)
    : numStates(states), lambda(arrivalRate), mu(serviceRate) {
    buildTransitionMatrix();
}

void MarkovChainAnalyzer::buildTransitionMatrix() {
    transitionMatrix.resize(numStates, std::vector<double>(numStates, 0.0));

    double totalRate = lambda + mu;
    double pArrival = lambda / totalRate;
    double pService = mu / totalRate;

    // State 0: can only have arrivals
    transitionMatrix[0][0] = 0.0;
    transitionMatrix[1][0] = 1.0;

    // States 1 to numStates-2: can have arrivals or departures
    for (int i = 1; i < numStates - 1; ++i) {
        transitionMatrix[i - 1][i] = pService;
        transitionMatrix[i + 1][i] = pArrival;
    }

    // State numStates-1: truncated (reflecting boundary)
    transitionMatrix[numStates - 2][numStates - 1] = pService;
    transitionMatrix[numStates - 1][numStates - 1] = pArrival;
}

std::vector<double> MarkovChainAnalyzer::computeSteadyState(int maxIterations, double tolerance) {
    std::vector<double> x(numStates, 1.0 / numStates);
    std::vector<double> xNew(numStates);

    for (int iter = 0; iter < maxIterations; ++iter) {
        for (int i = 0; i < numStates; ++i) {
            xNew[i] = 0.0;
            for (int j = 0; j < numStates; ++j) {
                xNew[i] += transitionMatrix[i][j] * x[j];
            }
        }

        double maxDiff = 0.0;
        for (int i = 0; i < numStates; ++i) {
            maxDiff = std::max(maxDiff, std::fabs(xNew[i] - x[i]));
        }

        x = xNew;

        if (maxDiff < tolerance) {
            std::cout << "Eigenvalue method converged in " << iter + 1 << " iterations."
                      << std::endl;
            break;
        }
    }

    // Normalize to ensure sum approaches 1
    double sum = 0.0;
    for (double val : x)
        sum += val;
    for (double& val : x)
        val /= sum;

    return x;
}

void MarkovChainAnalyzer::printTransitionMatrix() const {
    std::cout << "\n=== Transition Matrix (first 10x10 block) ===" << std::endl;
    int printSize = std::min(10, numStates);

    for (int i = 0; i < printSize; ++i) {
        for (int j = 0; j < printSize; ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(4)
                      << transitionMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
