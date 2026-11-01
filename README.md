# Polytope barcode tools

This repository contains MATLAB scripts and interactive HTML tools for simulating and analyzing Polytope epitope barcodes. Companion code for: *Polytope epitope barcoding reveals spatial clonal dynamics of T cell development in vivo* (manuscript submitted).

## Interactive tools

Open the [Polytope Barcode Tools website](https://dpostrach.github.io/Polytope/).

- [Polytope barcode recombination network](https://dpostrach.github.io/Polytope/interactive_tools/polytope_recombination_network.html) — explore the Polytope barcode recombination network.
- [Random barcode generation simulator](https://dpostrach.github.io/Polytope/interactive_tools/random_Polytope_barcode_generation.html) — simulate random barcode generation and spatial clonal outgrowth.

## MATLAB scripts

- `src/simulateColorCodes.m` — simulate Polytope color code generation given input cell numbers and recombination patterns.
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
