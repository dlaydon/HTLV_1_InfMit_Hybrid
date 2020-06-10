# HTLV-1 Hybrid Model

This repository contains the code for a model of within-host Human T-Lymphotropic Virus Type-1 (HTLV-1) persistence, developed by Daniel J. Laydon, Vikram Sunkara, Lies Boelen, Charles R.M. Bangham and Becca Asquith. Preprint and details available at <https://doi.org/10.1101/799197>

Script contains functions for hybrid model of chronic HTLV-1 infection/within host persistence. Model divides HTLV-1 proviral load (number of infected cells) into clones, where clones are defined as populations of indentically infected cells with common site of proviral integration. Clones proliferate via mitotic spread. Clones are created by infectious spread.

The hybrid model is comprised of two systems; i) deterministic system modelled by series of ODEs (for large clones); ii) stochastic system modelled by multiple birth-death processes (for smaller clones).  

The model is written in R.

## Scripts
- `CloneBirthDeath.R` contains functions for birth-death processes and the change in clone frequency probabilities over time.
- `CloneTrajectories.R` calls functcions from `CloneBirthDeath.R` to store look-up tables for the summary statistics of clone frequency probability distributions for each clone at each age or time. These tables are used when running the model.
- `HybridFunctions.R` contains functions for running the hybrid model.
- `QuickGridSearch.R` runs the hybrid model for a grid of values of the rate of infectious spread, given a fixed rate of mitotic spread.
- `OptimizeFit.R` fits the hybrid model via one-dimensional optimization, using values returned from `QuickGridSearch.R` to narrow the search space.
- `FitHybrid.R` calls i) `CloneTrajectories.R`, ii) `QuickGridSearch.R` and iii) `OptimizeFit.R` to fit the rate of infectious spread for a single patient data set.

## Inputs
The estiamted clone frequency distributions for each patient blood sample are given in directory [EstDists](./Inputs/EstDists). 

## Description
Patient blood sample characteristics are given in `PatientSampleInfo.txt`.