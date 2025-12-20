#include "../header/MarkovChainAnalysis.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

// constructor
MarkovChainAnalysis::MarkovChainAnalysis(int states, double arrivalRate, double serviceRate)
    : numStates(states), lambda(arrivalRate), miu(serviceRate) {
    buildGeneratorMatrix();
}

void MarkovChainAnalysis::buildGeneratorMatrix() {
    generatorMatrix.resize(numStates, std::vector<double>(numStates, 0.0));

    // State 0: can only have arrivals (no departures possible)
    generatorMatrix[0][0] = -lambda;
    generatorMatrix[0][1] = lambda;

    // States 1 to numStates-2: can have arrivals or departures
    for (int i = 1; i < numStates - 1; ++i) {
        generatorMatrix[i][i - 1] = miu;          // departure (service completion)
        generatorMatrix[i][i] = -(lambda + miu);  // diagonal (negative sum of rates)
        generatorMatrix[i][i + 1] = lambda;       // arrival
    }

    // State numStates-1: truncated boundary (can only have departures)
    generatorMatrix[numStates - 1][numStates - 2] = miu;
    generatorMatrix[numStates - 1][numStates - 1] = -miu;
}

std::vector<double> MarkovChainAnalysis::computeSteadyState(double tolerance) {
    std::vector<std::vector<double>> A(numStates, std::vector<double>(numStates + 1, 0.0));

    // Copy Q^T
    for (int i = 0; i < numStates - 1; ++i) {
        for (int j = 0; j < numStates; ++j) {
            A[i][j] = generatorMatrix[j][i];  // Transpose
        }
        A[i][numStates] = 0.0;  // RHS = 0 for πQ = 0
    }

    // Last row: normalization constraint Σπ_i = 1
    for (int j = 0; j < numStates; ++j) {
        A[numStates - 1][j] = 1.0;
    }
    A[numStates - 1][numStates] = 1.0;  // RHS = 1

    for (int col = 0; col < numStates; ++col) {
        // Find pivot
        int maxRow = col;
        double maxVal = std::abs(A[col][col]);
        for (int row = col + 1; row < numStates; ++row) {
            if (std::abs(A[row][col]) > maxVal) {
                maxVal = std::abs(A[row][col]);
                maxRow = row;
            }
        }

        // Swap rows
        std::swap(A[col], A[maxRow]);

        // Check for singular matrix
        if (std::abs(A[col][col]) < tolerance) {
            std::cerr << "Warning: Near-singular matrix in steady-state computation" << std::endl;
            continue;
        }

        // Eliminate column
        for (int row = col + 1; row < numStates; ++row) {
            double factor = A[row][col] / A[col][col];
            for (int k = col; k <= numStates; ++k) {
                A[row][k] -= factor * A[col][k];
            }
        }
    }

    // Back substitution
    std::vector<double> pi(numStates, 0.0);
    for (int i = numStates - 1; i >= 0; --i) {
        pi[i] = A[i][numStates];
        for (int j = i + 1; j < numStates; ++j) {
            pi[i] -= A[i][j] * pi[j];
        }
        pi[i] /= A[i][i];
    }

    // Ensure non-negative probabilities and normalize
    double sum = 0.0;
    for (int i = 0; i < numStates; ++i) {
        if (pi[i] < 0)
            pi[i] = 0;
        sum += pi[i];
    }
    for (int i = 0; i < numStates; ++i) {
        pi[i] /= sum;
    }

    std::cout << "Steady-state computed via Gaussian elimination on generator matrix" << std::endl;

    return pi;
}

std::vector<double> MarkovChainAnalysis::getTheoreticalSteadyState() const {
    std::vector<double> pi(numStates, 0.0);
    double rho = lambda / miu;

    if (rho >= 1.0) {
        std::cerr << "System is unstable (ρ >= 1), no steady-state exists" << std::endl;
        return pi;
    }

    for (int n = 0; n < numStates; ++n) {
        pi[n] = (1.0 - rho) * std::pow(rho, n);
    }

    return pi;
}

void MarkovChainAnalysis::printGeneratorMatrix() const {
    std::cout << "\n=== Generator Matrix Q (first 10x10 block) ===" << std::endl;
    std::cout << "CTMC rate matrix for M/M/1 queue (λ=" << lambda << ", μ=" << miu << ")"
              << std::endl;
    int printSize = std::min(10, numStates);

    for (int i = 0; i < printSize; ++i) {
        for (int j = 0; j < printSize; ++j) {
            std::cout << std::setw(9) << std::fixed << std::setprecision(4) << generatorMatrix[i][j]
                      << " ";
        }
        std::cout << std::endl;
    }
}
