# Guided Random Rearrangement for Store Layout Simulation

## Requirements

- MATLAB R2021a or later
- Required `.mat` files (see below)
- Sufficient memory and computation time for multi-period simulations

## Overview

`guidedRandomRearrangement.m` simulates retail store layout evolution using guided randomization. It evaluates and adjusts department-subcategory assignments over multiple periods using historical data and profitability metrics. Unlike optimization-based approaches, this method uses randomized shuffles with practical constraints to mimic real-world layout experimentation.

## Key Features

- Multi-trial, multi-period layout simulation
- Metrics for profitability, habitual behavior, and cross-selling
- Uses seed layouts and department-specific space constraints
- No solver or external modeling language required

## Function

Function: guidedRandomRearrangement(NumberofTrials, NumberofPeriods)

Inputs:
- NumberofTrials: Number of independent layout simulation runs
- NumberofPeriods: Number of time periods (e.g., months)

Outputs (via global variables):
- Analysis4_MetricsOverTime: Performance metrics for each trial/period
- Analysis4_AislesAssigned: Aisle assignment results per period
- Analysis4_FinalShelfAllocation: Shelf allocation details per period

## Required Data Files

Before running, load the following `.mat` files into the MATLAB workspace:

- AssociationData.mat
- InitialLayouts.mat
- Subcategory_FP_Growth.mat
- Profitability_Matrix.mat
- RawTransactionData.mat
- Department_FP_Growth.mat

## How It Works

1. For each trial:
   - Period 1 uses a seed layout.
   - Subsequent periods apply guided random rearrangement based on:
     - Subcategory space constraints
     - Department size
     - Shelf segment availability

2. Metrics are calculated per layout:
   - Metric 1: Visibility-driven expected profit
   - Metric 2: Historical impulse (habitual aisle continuity)
   - Metric 3: Cross-selling potential based on intra-shelf proximity

## Usage Example

From the MATLAB command window:

guidedRandomRearrangement(1000, 20)

## Output

Save outputs manually after execution:

- Analysis4_MetricsOverTime.mat
- Analysis4_AislesAssigned.mat
- Analysis4_FinalShelfAllocation.mat

## Notes

- No AMPL or solver setup required.
- Randomization is controlled for feasibility and diversity.
- Suitable for large-scale retail layout simulations with limited optimization infrastructure.
