# Learning Dantzig’s Pivot Rule with a Neural Network

This project is an exercise in building a data pipeline around the Simplex algorithm. I implement Dantzig’s pivot rule, generate supervised data from random linear programs (LPs), train a small neural network to predict the entering variable, and evaluate the learned rule against Dantzig.

## Overview
I solve LPs in the form: maximize cᵀx subject to Ax ≤ b, x ≥ 0. The focus is the pivot rule (entering-variable choice). Dantzig’s rule picks the largest reduced cost; the model learns to imitate this decision from data.

## What’s Included
- Script 1 — Simplex (baseline reference).
- Script 2 — Data generation: run Simplex with Dantzig; log per-pivot (X = reduced costs, y = chosen column) → pivot_data_big.npz.
- Script 3 — Training: PyTorch MLP, feature normalization, cross-entropy on pivot labels → pivot_model.pt, norm_stats.npz.
- Script 4 — Evaluation: compare NN vs Dantzig on fresh LPs; reports mean pivots, 95% CI, min/max, and win/tie counts.
- Script 5 — Visualization: scatter (NN vs Dantzig), histogram of (NN−Dantzig) differences, ECDF; saves PNGs.

## Concepts & Tools
Simplex • Dantzig’s rule • supervised learning • PyTorch (MLP) • NumPy (tableau + data) • Matplotlib (plots) • reproducibility via fixed RNG seeds • cached results (.npz).

## How to Run
1) Generate data: `python script2_generate_data.py`  
2) Train model: `python script3_train_nn.py`  
3) Evaluate: `python script4_evaluate.py`  
4) Visualize: `python script5_plot.py`

## Notes
- LP size is controlled by m (constraints) and n (variables); larger sizes yield longer pivot paths and clearer differences.  
- Random Gaussian LPs often make Dantzig hard to beat; expect the NN to be close, not better.  
- This is an educational project intended to demonstrate a clean data → train → evaluate → visualize pipeline around a classic algorithm.
