# multi-scale-GJN

Code for the paper: *Asymptotic product-form steady-state for generalized Jackson networks in multi-scale heavy traffic.*

> **Paper:** [arXiv:2304.01499](https://arxiv.org/abs/2304.01499)

This repository contains:

- **simulation/** – C++ source codes for running simulations.
- **simulated_data/** – Output data generated from the simulations.
- **data_visualization/** – Python scripts for analyzing and visualizing the simulated data.

---

## for simulation running 
g++ -O3 -march=native -std=c++17 -DNDEBUG -o simulate gjn_multi_general.cpp
./simulate

## for indepedence
g++ -O3 -march=native -std=c++17 -DNDEBUG -o gjn_multi_general_independence_representative gjn_multi_general_independence_representative.cpp
./gjn_multi_general_independence_representative
