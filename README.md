# On-the-Fly Crystal

**An automated high-performance framework for constructing potential energy surfaces of transition metal complexes**

[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fjcc.27285-blue)](https://doi.org/10.1002/jcc.27285)
[![Language](https://img.shields.io/badge/Language-Fortran_90-green)]()
[![HPC](https://img.shields.io/badge/HPC-SLURM%20%7C%20Frontera-orange)]()

## Overview

On-the-Fly Crystal is an automated computational framework written in **Fortran 90** that constructs potential energy surfaces (PES) for transition metal complexes by coupling intelligent lattice-growth algorithms with automated DFT calculations on HPC clusters. The code systematically explores multidimensional normal-mode coordinate space, interfaces directly with **Gaussian 16** for electronic structure calculations, and uses a **binary sort tree** for efficient duplicate detection, enabling the evaluation of **80,000+ molecular structures in under 72 hours** with a **20× speedup** over manual approaches.

## The Problem

Constructing accurate potential energy surfaces for transition metal polyhydride complexes is one of the most computationally demanding tasks in quantum chemistry. Each point on the PES requires a full DFT calculation, and the number of required points scales exponentially with dimensionality. Traditional approaches require researchers to manually generate geometries, prepare input files, submit jobs, monitor completion, handle failures, and curate results, a process that takes weeks and is prone to human error.

## How It Works

On-the-Fly Crystal automates this entire workflow through four tightly integrated components:

### 1. Simplicial Growth Algorithm
The code uses a **simplicial growth strategy** to systematically expand a lattice of grid points outward from an equilibrium geometry in N-dimensional normal-mode coordinate space. At each iteration, new candidate points are generated along each coordinate direction (positive and negative), with support for **symmetric coordinates** that skip redundant negative directions.

### 2. Binary Sort Tree for Duplicate Detection
All evaluated lattice points are stored in a custom **binary sort tree** data structure, enabling O(log n) lookup to determine whether a candidate point has already been evaluated. This prevents redundant DFT calculations and ensures efficient memory usage even for datasets exceeding 80,000 structures. A separate tree tracks rejected points above the energy cutoff.

### 3. Automated Gaussian 16 Interface
For each new lattice point, the code:
- Applies a **coordinate transformation matrix** to convert normal-mode displacements into Cartesian atomic coordinates
- Automatically generates a Gaussian 16 input file (`.gjf`) with the correct molecular geometry, basis set (def2-TZVPP), DFT functional (TPSSh), and HPC resource allocation
- Monitors the Gaussian output log in real time using a **sleep-wait polling loop**, detecting SCF convergence, reading final energies, and handling error conditions (e.g., distance matrix failures)
- Parses the SCF energy from the output and applies an **energy cutoff criterion** to accept or reject the lattice point

### 4. Energy-Based Acceptance/Rejection
Each evaluated point is classified based on its DFT energy relative to a configurable cutoff (`Ecut0`). Points below the cutoff are added to the accepted lattice (`gog` tree) and written to the output file. Points above the cutoff are stored in a separate rejected lattice (`nogo` tree) to prevent re-evaluation. This ensures computational resources are focused on the chemically relevant region of the PES.

## Technical Specifications

| Feature | Detail |
|---------|--------|
| Language | Fortran 90 |
| Compiler | Intel Fortran (`ifort`) |
| DFT Engine | Gaussian 16 (TPSSh/def2-TZVPP) |
| HPC Scheduler | SLURM (`sbatch`) |
| Tested On | TACC Frontera (56 cores/node, 64 GB memory) |
| Dimensionality | Configurable N-dimensional (currently 2D–6D) |
| Max Structures | 9,999,999 (configurable allocation) |
| Duplicate Detection | Binary sort tree, O(log n) |
| Lattice Strategy | Simplicial growth with symmetric coordinate support |
| Output | `outputLattices.txt`, coordinates and DFT energies |

## File Structure

```
On-the-Fly-Crystal/
├── crystal.f90            # Main source code
├── submit.sh              # SLURM job submission script
├── gin.gjf                # Gaussian input (auto-generated at runtime)
├── gin.log / gw.log       # Gaussian output (monitored at runtime)
├── radiux.txt             # Growth direction vectors (auto-generated)
├── outputLattices.txt     # Final output: coordinates + energies
└── README.md
```

## Quick Start

### Prerequisites
- Intel Fortran compiler (`ifort`)
- Gaussian 16
- SLURM-managed HPC cluster

### Setup
1. Set the equilibrium energy in `gin.log` (SCF Done format):
   ```
   SCF Done:  E(RTPSSh) =  -2638.09
   ```

2. Configure dimensionality and parameters in `crystal.f90`:
   ```fortran
   d = 2                    ! Number of normal-mode dimensions
   Ecut0 = -2638.09         ! Energy cutoff (Hartrees)
   symmetric(1) = 1         ! Set to 1 for symmetric coordinates
   ```

3. Compile and submit:
   ```bash
   ifort crystal.f90 -o crystal
   sbatch submit.sh
   ```

4. Results are written to `outputLattices.txt` as coordinate-energy pairs.

## Key Algorithms

### Binary Sort Tree (`sort` subroutine)
Custom implementation for multidimensional vector comparison using lexicographic ordering. Returns whether a vector is new (`unique = ±1`) or already exists (`unique = 0`) in the tree, along with the insertion point (`iii`).

### Simplicial Growth (`radiux.txt` generation)
Automatically constructs growth direction vectors based on dimensionality `d`. For each dimension, generates unit vectors in positive and (if not symmetric) negative directions. The origin point is always included as the seed.

### Gaussian Monitor (sleep-wait loop)
Polls the Gaussian log file at 1-second intervals, verifying:
1. Correct geometry was used (coordinate matching within 0.001 tolerance)
2. SCF convergence was achieved
3. No fatal errors (e.g., distance matrix problems)

## Test System

The default configuration models an **iron polyhydride complex** \[FeH₂(PH₃)₄\]⁺ a transition metal system with:
- Iron center with 15 hydrogen atoms and 4 phosphorus ligands
- Charge: +1, Multiplicity: 1 (singlet)
- DFT: TPSSh functional, def2-TZVPP basis set
- 2D PES along selected normal-mode coordinates

This system was chosen for its relevance to hydrogen exchange mechanisms in organometallic chemistry, as published in *J. Phys. Chem. A* (2023).

## Publications

1. **M. Aarabi**, J. Eckert, B. Poirier, "Automated PESs up to Six-Dimensional for Transition Metal Complex Using Parallel On-the-fly Crystal Code," *J. Chem. Theory Comput.* (2026), Submitted.

2. **M. Aarabi**, A. Pandey, B. Poirier, "On-the-fly Crystal: How to reliably and automatically characterize and construct PESs," *J. Comput. Chem.* 45 (2024) 1261–1278. [DOI: 10.1002/jcc.27285](https://doi.org/10.1002/jcc.27285)

3. **M. Aarabi** et al., "Quantum dynamical investigation of dihydrogen–hydride exchange in a transition metal polyhydride complex," *J. Phys. Chem. A* 127 (2023) 6385–6399.

## Author

**Mahdi Aarabi, Ph.D.**  
Computational Scientist 

## Citation

If you use this code in your research, please cite:

```bibtex
@article{aarabi2024crystal,
  title={On-the-fly Crystal: How to reliably and automatically characterize 
         and construct potential energy surfaces},
  author={Aarabi, Mahdi and Pandey, Amit and Poirier, Bill},
  journal={Journal of Computational Chemistry},
  volume={45},
  pages={1261--1278},
  year={2024},
  doi={10.1002/jcc.27285}
}
```

