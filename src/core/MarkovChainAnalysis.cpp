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

void MarkovChainAnalysis::gaussianElimination(std::vector<std::vector<double>>& A,
                                              double tolerance) {
    int n = static_cast<int>(A.size());

    for (int col = 0; col < n; ++col) {
        // Find pivot
        int maxRow = col;
        double maxVal = std::abs(A[col][col]);
        for (int row = col + 1; row < n; ++row) {
            if (std::abs(A[row][col]) > maxVal) {
                maxVal = std::abs(A[row][col]);
                maxRow = row;
            }
        }

        // Swap rows to bring pivot to diagonal
        std::swap(A[col], A[maxRow]);

        // Check for singular or near-singular matrix
        if (std::abs(A[col][col]) < tolerance) {
            std::cerr << "Warning: Near-singular matrix at column " << col << std::endl;
            continue;
        }

        // Eliminate entries below pivot
        for (int row = col + 1; row < n; ++row) {
            double factor = A[row][col] / A[col][col];
            for (int k = col; k <= n; ++k) {  // n+1 columns (augmented)
                A[row][k] -= factor * A[col][k];
            }
        }
    }
}

std::vector<double> MarkovChainAnalysis::backSubstitution(
    const std::vector<std::vector<double>>& A) {
    int n = static_cast<int>(A.size());
    std::vector<double> x(n, 0.0);

    for (int i = n - 1; i >= 0; --i) {
        x[i] = A[i][n];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    return x;
}

void MarkovChainAnalysis::normalizeDistribution(std::vector<double>& pi) {
    double sum = 0.0;

    // Ensure non-negative probabilities
    for (double& p : pi) {
        if (p < 0)
            p = 0;
        sum += p;
    }

    // Normalize to sum to 1
    if (sum > 0) {
        for (double& p : pi) {
            p /= sum;
        }
    }
}

std::vector<double> MarkovChainAnalysis::computeSteadyState(double tolerance) {
    // Build augmented matrix: Q^T with last row replaced by normalization
    std::vector<std::vector<double>> A(numStates, std::vector<double>(numStates + 1, 0.0));

    // Copy Q^T (transpose of generator matrix) for rows 0 to numStates-2
    for (int i = 0; i < numStates - 1; ++i) {
        for (int j = 0; j < numStates; ++j) {
            A[i][j] = generatorMatrix[j][i];  // Transpose
        }
        A[i][numStates] = 0.0;
    }

    // Last row: normalization constraint (sum = 1)
    for (int j = 0; j < numStates; ++j) {
        A[numStates - 1][j] = 1.0;
    }
    A[numStates - 1][numStates] = 1.0;

    // Solve using Gaussian elimination
    gaussianElimination(A, tolerance);
    std::vector<double> pi = backSubstitution(A);
    normalizeDistribution(pi);

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
