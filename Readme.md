# Simulation code: Group-sparse Sensor Selection for GEVD problems

This repository contains the simulation code that was used to produce the results presented in [1]. It is a MATLAB implementation of group-sparse sensor selection methods for GEVD problems.

The repository contains an implementation of our proposed method along with an implementation for:
- Exhaustive search
- Forward greedy selection
- Backward greedy selection
- S. A. Hamza, M. G. Amin, Sparse Array Beamforming Design for Wideband Signal Models, IEEE Trans. Aerosp.
Electron. Syst. 57 (2) (2021) 1211-1226. doi:10.1109/TAES.2020.3037409
- F. Qi, W. Wu, Z. L. Yu, Z. Gu, Z. Wen, T. Yu, Y. Li, Spatiotemporal-Filtering-Based Channel Selection for
Single-Trial EEG Classification, IEEE Trans. Cybern. 51 (2) (2021) 558-567. doi:10.1109/TCYB.2019.2963709.

## Requirements
This package requires:
- [CVX](http://cvxr.com/cvx/)
- [minFunc](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)

The software was developed with CVX version 2.2 using the MOSEK optimizer.

## Usage

Two simulations are run in the paper:

- Number of channels (C) = 25; Number of time lags (L) = 3; Number of output filters (K) = 1. This simulation is stored in the folder `rank-deficient-r1-c25-l3-k1`. To generate results for this scenario call `main(seed)` by providing a random seed number.
- Number of channels (C) = 25; Number of time lags (L) = 2; Number of output filters (K) = 2. This simulation is stored in the folder `rank-deficient-r1-c25-l2-k2`. To generate results for this scenario call `main(seed)` by providing a random seed number.

## Reference

[1] J. Dan, S. Geirnaert, A. Bertrand, "Grouped Variable Selection for Generalized Eigenvalue Problems," 2021. arXiv:2105.13667 (https://arxiv.org/abs/2105.13667)
