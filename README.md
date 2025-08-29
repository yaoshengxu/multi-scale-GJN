# multi-scale-GJN

Code for the paper: *Asymptotic product-form steady-state for generalized Jackson networks in multi-scale heavy traffic.*

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
