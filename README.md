# Hofstadter Lattice Model

A Julia implementation for studying the Hofstadter model and related quantum many-body physics on lattices. This repository contains tools for exact diagonalization of the Hubbard-Hofstadter model and analysis of Laughlin wavefunctions.

## Overview

The Hofstadter model describes non-interacting electrons on a 2D lattice in the presence of a uniform magnetic field. This implementation extends the model to include:

- **Hubbard interaction**: On-site Coulomb repulsion
- **Hard-core bosons**: Infinite repulsion limit
- **Laughlin wavefunctions**: Fractional quantum Hall states
- **Exact diagonalization**: Full many-body spectrum calculation

## Project Structure

```
├── Main Scripts/           # Core model implementation
│   ├── Lattice.jl         # Lattice geometry and neighbor functions
│   ├── Model.jl           # Hofstadter and Hubbard model Hamiltonians
│   ├── MBBasis.jl         # Many-body basis construction
│   ├── Operators.jl       # Quantum operators (creation/annihilation)
│   ├── ED.jl              # Exact diagonalization routines
│   └── Utilities.jl       # Helper functions and utilities
├── Laughlin Scripts/       # Fractional quantum Hall states
│   ├── GeneralizedLaughlin.jl  # Laughlin wavefunction construction
│   └── JacobiThetaFunction.jl # Theta function implementations
├── SolveModel.jl          # Main script for model solving
└── Overlap.jl             # Wavefunction overlap analysis
```

## Key Features

### 🧲 Hofstadter Model
- Single-particle Hofstadter Hamiltonian with magnetic flux `α`
- Periodic boundary conditions with magnetic phase factors
- Configurable lattice sizes (Nx × Ny)

### ⚛️ Many-Body Physics
- **Hubbard interaction**: On-site repulsion parameter `U`
- **Hard-core constraint**: Infinite repulsion limit for bosons
- **Particle number conservation**: Fixed particle number `pn`
- **Exact diagonalization**: Full spectrum calculation

### 🌊 Laughlin Wavefunctions
- Generalized Laughlin states for fractional quantum Hall effect
- Center-of-mass and relative coordinate separation
- Overlap analysis with exact diagonalization eigenstates
- Jacobi theta function implementations

### 📊 Analysis Tools
- Real-space density visualization
- Wavefunction overlap calculations
- Coefficient optimization routines
- Hilbert-Schmidt norm calculations

## Usage

### Basic Model Solving

```julia
# Include the main scripts
include("SolveModel.jl")

# Set parameters
pn = 2          # Number of particles
Nx = 4          # Lattice size in x-direction
Ny = 4          # Lattice size in y-direction
alpha = 1/4     # Magnetic flux per plaquette
U = 1           # Hubbard interaction strength
periodicity = true
HardCore = true
Nev = 10        # Number of eigenvalues to compute

# Solve the model
E, ψ = Solve(pn, Nx, Ny, alpha, periodicity, HardCore, U, Nev)
```

### Laughlin Wavefunction Analysis

```julia
# Include overlap analysis script
include("Overlap.jl")

# Generate Laughlin states
ψ0, ψ1 = GeneralizedLaughlin(basis, Nx, Ny, UpperLimit, type)

# Calculate overlaps with exact diagonalization states
overlap_01 = Overlap(ψ0, ψ[:,1])
overlap_02 = Overlap(ψ0, ψ[:,2])
```

## Dependencies

- **Julia**: Version 1.6 or higher
- **QuantumOptics.jl**: Quantum state manipulation
- **Plots.jl**: Visualization
- **LaTeXStrings.jl**: LaTeX formatting in plots
- **Combinatorics.jl**: Combinatorial functions

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/hofstadter-lattice-model.git
cd hofstadter-lattice-model
```

2. Install Julia dependencies:
```julia
using Pkg
Pkg.add(["QuantumOptics", "Plots", "LaTeXStrings", "Combinatorics"])
```

## Examples

### Energy Spectrum
The model can compute the full many-body energy spectrum, revealing:
- Hofstadter butterfly structure in single-particle limit
- Interaction-induced gaps and level repulsion
- Ground state properties and phase transitions

### Wavefunction Analysis
- Compare exact diagonalization eigenstates with Laughlin ansätze
- Analyze real-space density distributions
- Study quantum Hall physics in finite systems

## Physics Background

The Hofstadter model describes electrons on a square lattice with a uniform magnetic field. The Hamiltonian is:

```
H = -t ∑ᵢⱼ (cᵢ†cⱼ e^(iφᵢⱼ) + h.c.) + U ∑ᵢ nᵢ↑nᵢ↓
```

where:
- `t` is the hopping amplitude
- `φᵢⱼ` are magnetic phase factors
- `U` is the on-site Hubbard interaction
- `α = p/q` is the magnetic flux per plaquette

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This implementation builds on the rich history of Hofstadter model studies and fractional quantum Hall physics. Special thanks to the Julia community and quantum physics researchers who have inspired this work.
