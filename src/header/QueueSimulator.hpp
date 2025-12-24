#ifndef MM1QUEUESIMULATOR_H
#define MM1QUEUESIMULATOR_H

#include <queue>
#include <random>
#include <string>
#include <vector>

#include "Event.hpp"

class QueueSimulator {
   private:
    double lambda;
    double mu;
    double rho;

    double currentTime;
    int currentQueueLength;
    double totalSimulationTime;
    int maxStateTracked;
    std::vector<double> timeInState;

    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> eventQueue;

    std::mt19937 rng;
    std::exponential_distribution<double> arrivalDist;
    std::exponential_distribution<double> serviceDist;

    void scheduleArrival();
    void handleArrival();
    void handleDeparture();

   public:
    QueueSimulator(double arrivalRate, double serviceRate, int maxStates = 100,
                   unsigned int seed = 42);

    void simulate(double duration);

    std::vector<double> getSteadyStateProbabilities() const;
    std::vector<double> getTheoreticalProbabilities() const;
    double getAverageQueueLength() const;
    double getTheoreticalAverageQueueLength() const;

    void exportToCSV(const std::string& filename) const;
    void printSummary() const;

    // Getters
    double getRho() const {
        return rho;
    }
    double getLambda() const {
        return lambda;
    }
    double getMu() const {
        return mu;
    }
};

#endif
