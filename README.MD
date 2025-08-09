# Simplex Algorithm (MATLAB and Python)

This project provides implementations of the Simplex algorithm for solving linear programs in standard form:

  maximize cᵗx  
  subject to Ax = b, x ≥ 0

## Input:
- A: m × n matrix
- b: m × 1 vector
- c: n × 1 vector
- BAS: m × 1 vector of indices representing a feasible basis (1-based in MATLAB, 0-based in Python)

## Output:
Both implementations return [stat, sol]:
- stat = 0: optimal solution found → sol is the optimal x
- stat = 1: problem is unbounded → sol is an unbounded direction

## Files:
- `simplex.m`: MATLAB implementation
- `text.m`: MATLAB test cases
- `simplex.py`: Python version, ready for Google Colab

## Example (common to both):
maximize x₁ + 4x₂  
subject to:  
 x₁ + 2x₂ ≤ 2  
 x₁ + x₂ ≤ 4  
→ converted to standard form with slack variables

Inputs:
- c = [1; 4; 0; 0]
- A = [1 2 1 0; 1 1 0 1]
- b = [2; 4]
- BAS = [1; 4] in MATLAB or [0, 3] in Python

## Usage:
MATLAB:
```matlab
[stat, sol] = simplex(c, A, b, BAS);
