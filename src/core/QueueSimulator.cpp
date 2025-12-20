#include "../header/QueueSimulator.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>

QueueSimulator::QueueSimulator(double arrivalRate, double serviceRate, 
                                     int maxStates, unsigned int seed)
    : lambda(arrivalRate), mu(serviceRate), currentTime(0.0), 
      currentQueueLength(0), maxStateTracked(maxStates), 
      rng(seed), arrivalDist(arrivalRate), serviceDist(serviceRate) {
    
    rho = lambda / mu;
    timeInState.resize(maxStateTracked + 1, 0.0);
    
    if (rho >= 1.0) {
        std::cerr << "WARNING: System is unstable (λ/μ = " << rho << " >= 1)" << std::endl;
        std::cerr << "Queue will grow indefinitely over time." << std::endl;
    }
}

void QueueSimulator::simulate(double duration) {
  totalSimulationTime = duration;

  scheduleArrival();

  while (!eventQueue.empty() && currentTime < duration) {
    Event nextEvent = eventQueue.top();
    eventQueue.pop();

    if (nextEvent.time > duration) {
      int state = std::min(currentQueueLength, maxStateTracked);
      timeInState[state] += (duration - currentTime);
      currentTime = duration;
    }

    int state = std::min(currentQueueLength, maxStateTracked);
    timeInState[state] += (nextEvent.time - currentTime);
    
    currentTime = nextEvent.time;
    
    if (nextEvent.type == ARRIVAL) {
      handleArrival();
    } else {
      handleDeparture();
    }
  }
}

void QueueSimulator::scheduleArrival() {
    double interArrivalTime = arrivalDist(rng);
    eventQueue.push({currentTime + interArrivalTime, ARRIVAL});
}

void QueueSimulator::handleArrival() {
  currentQueueLength++;

  scheduleArrival();

  if (currentQueueLength == 1) {
    double serviceTime = serviceDist(rng);
    eventQueue.push({currentTime + serviceTime, DEPARTURE});
  } 
}

void QueueSimulator::handleDeparture() {
    if (currentQueueLength > 0) {
      currentQueueLength--;
      
    if (currentQueueLength > 0) {
      double serviceTime = serviceDist(rng);
      eventQueue.push({currentTime + serviceTime, DEPARTURE});
    }
  }
}

std::vector<double> QueueSimulator::getSteadyStateProbabilities() const {
    std::vector<double> probabilities(maxStateTracked + 1);
    
    for (int i = 0; i <= maxStateTracked; ++i) {
      probabilities[i] = timeInState[i] / totalSimulationTime;
    }
    
    return probabilities;
}

std::vector<double> QueueSimulator::getTheoreticalProbabilities() const {
    std::vector<double> probabilities(maxStateTracked + 1);
    
    if (rho >= 1.0) {
      std::cerr << "Cannot compute theoretical probabilities for unstable system." << std::endl;
      return probabilities;
    }
    
    for (int n = 0; n <= maxStateTracked; ++n) {
      probabilities[n] = (1 - rho) * std::pow(rho, n);
    }
    
    return probabilities;
}

double QueueSimulator::getAverageQueueLength() const {
    double avgLength = 0.0;
    
    for (int i = 0; i <= maxStateTracked; ++i) {
      avgLength += i * (timeInState[i] / totalSimulationTime);
    }
    
    return avgLength;
}

void QueueSimulator::printSummary() const {
    std::cout << "\n=== M/M/1 Queue Simulation Results ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Arrival rate (λ): " << lambda << std::endl;
    std::cout << "Service rate (μ): " << mu << std::endl;
    std::cout << "Utilization (ρ): " << rho << std::endl;
    std::cout << "Simulation time: " << totalSimulationTime << std::endl;
    std::cout << "System status: " << (rho < 1.0 ? "STABLE" : "UNSTABLE") << std::endl;
    
    std::cout << "\n--- Average Queue Length ---" << std::endl;
    std::cout << "Simulated: " << getAverageQueueLength() << std::endl;
    
    if (rho < 1.0) {
      std::cout << "Theoretical: " << getTheoreticalAverageQueueLength() << std::endl;
      std::cout << "Relative error: " 
        << std::abs(getAverageQueueLength() - getTheoreticalAverageQueueLength()) 
            / getTheoreticalAverageQueueLength() * 100
        << "%" << std::endl;
    }
    
    std::cout << "\n--- Steady-State Distribution (first 10 states) ---" << std::endl;
    std::cout << "State | Simulated  | Theoretical | Abs. Error" << std::endl;
    std::cout << "------|------------|-------------|------------" << std::endl;
    
    auto simProbs = getSteadyStateProbabilities();
    auto theoProbs = getTheoreticalProbabilities();
    
    for (int i = 0; i < std::min(10, maxStateTracked + 1); ++i) {
        std::cout << std::setw(5) << i << " | "
                  << std::setw(10) << simProbs[i] << " | "
                  << std::setw(11) << theoProbs[i] << " | "
                  << std::setw(10) << std::abs(simProbs[i] - theoProbs[i])
                  << std::endl;
    }
}