# Polytope MATLAB scripts
MATLAB scripts for simulating and calling Polytope epitope barcodes. Companion code for: *Polytope epitope barcoding reveals spatial clonal dynamics of T cell development in vivo* (manuscript submitted to *Nature*).

## Contents
- `src/simulateColorCodes.m` — simulate Polytope color-code generation given input cell numbers and recombination patterns.
- `src/runGMMPolytope.m` — call color codes from MFI data using Gaussian Mixture Models (GMM).
- `example_data/input_for_GMM.csv` — example input file for `runGMMPolytope.m`.

## Requirements
- MATLAB R20XXx or later.
- Statistics and Machine Learning Toolbox (for GMM functions like `fitgmdist`).

## Quick start
```matlab
% From MATLAB, set path to src
addpath('src');

% 1) Simulation
% (Edit the script parameters inside simulateColorCodes.m as needed, then run:)
simulateColorCodes

% 2) GMM-based code calling
data = readmatrix('example_data/input_for_GMM.csv');
% runGMMPolytope.m may expect a filename or a matrix; adjust as needed:
runGMMPolytope('example_data/input_for_GMM.csv');  % or runGMMPolytope(data)


