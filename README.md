# M/M/1 Queue Analysis using Markov Chains

Analysis of Steady-State Behavior in Server Queues using Markov Chains and Eigenvalues in the M/M/1 Model.

This project accompanies a research paper for the **IF2123 Linear and Geometrical Algebra** course I took in semester 3.


## Disclaimer

I am by no means an expert, therefore this project is not more than a personal attempt of implementing the research I've done.

## Overview

This project analyzes the steady-state behavior of M/M/1 queues using:
- **Continuous-Time Markov Chain (CTMC)** theory
- **Generator matrix** formulation and balance equations
- **Discrete-event simulation** for validation

The M/M/1 queue models a single-server system where:
- Arrivals follow a Poisson process (rate λ)
- Service times are exponentially distributed (rate μ)
- The system is stable when ρ = λ/μ < 1

## Project Structure

```
makalah-algeo/
├── docs/                           # LaTeX paper and figures
│   ├── main.tex                    # Main paper source
│   ├── main.pdf                    # Compiled paper
│   └── *.png, *.jpg                # Figures and images
├── src/                            # C++ simulation source code
│   ├── main.cpp                    # Main driver program
│   ├── Makefile                    # Build configuration
│   ├── core/                       # Core implementation
│   │   ├── QueueSimulator.cpp      # Discrete-event simulator
│   │   ├── MarkovChainAnalysis.cpp # CTMC analysis
│   │   └── ExperimentRunner.cpp    # Experiment automation
│   ├── header/                     # Header files
│   │   ├── Event.hpp
│   │   ├── QueueSimulator.hpp
│   │   ├── MarkovChainAnalysis.hpp
│   │   └── ExperimentRunner.hpp
│   ├── scripts/                    # Python visualization
│   │   └── plot_mm1_results.py
│   ├── data/                       # Simulation output (CSV)
│   └── plots/                      # Generated plots
└── README.md
```

## Requirements

### C++ Simulation
- C++17 compatible compiler (g++, clang++)
- Make build system

### Python Visualization
- Python 3.8+
- pandas
- matplotlib
- numpy

## Building and Running

### Build the Simulation

```bash
cd src
make build
```

### Run the Simulation

```bash
# Run with default parameters (λ=0.8, μ=1.0, T=100000)
make run

# Or run with custom parameters
./bin/queue_simulator [lambda] [mu] [sim_time] [max_states] [seed]

# Example: High utilization test
./bin/queue_simulator 0.95 1.0 100000 50 42
```

### Generate Plots

```bash
cd src
python scripts/plot_mm1_results.py
```

### Clean Build Artifacts

```bash
make clean        # Clean all build artifacts
make clean-obj    # Clean only object files
make clean-output # Clean only CSV output files
```

## Paper

The full research paper is available in `docs/main.pdf`, covering:
- Theoretical foundations (matrices, eigenvalues, Markov chains)
- M/M/1 queue model and CTMC representation
- Steady-state analysis via generator matrix
- Implementation details
- Experimental results and discussion

## Author

**Nicholas Wise Saragih Sumbayak** (13524037)
