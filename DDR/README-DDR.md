# AMPL Models for Store Layout Optimization

## Requirements

- AMPL installed and licensed
- Solvers: BARON (for Step 1) and CPLEX (for Steps 2 and 3)
- MATLAB AMPL API bindings configured (run `setupOnce` in MATLAB)
- Ensure paths to .mod files are correct when running the MATLAB function

## Overview

This repository contains AMPL model files and input files (.mat) used in the `dataDrivenRearrangement.m` simulation for retail store layout optimization. The layout evolves across multiple periods and is optimized through three sequential steps.

Each step calls a separate AMPL model via MATLAB. All input data is prepared and passed in-code; no .dat files are needed.

## Steps and Models

Step 1: Virtual Category Grouping  
Purpose: Groups subcategories into virtual categories within each department to maximize intra-group profitability and association strength.  
AMPL model: Step1_Grouping.mod  
Solver: BARON  
Inputs: Subcategory list, profitability vector, association matrix  
Output: Binary assignment of subcategories to virtual categories

Step 2: Shelf Space Allocation  
Purpose: Allocates shelf space to each virtual category based on capacity and layout constraints.  
AMPL models:  
- Step2_SingleShelf.mod (for aisles with 1 shelf)  
- Step2_TwoShelf.mod (for aisles with 2 shelves)  
Solver: CPLEX  
Inputs: Shelf segment list, subcategory space bounds, profitability and support per virtual category  
Output: Allocation of shelf segments to virtual categories

Step 3: Aisle Assignment  
Purpose: Assigns departments to aisles based on feasibility and layout history  
AMPL model: Step3_AisleAssign.mod  
Solver: CPLEX  
Inputs: Feasibility matrix, past aisle assignments  
Output: Department-to-aisle assignment

## How to Use

These models are invoked automatically from MATLAB. To use them in the full simulation, provide the paths to each .mod file when calling the main function:

Example call from MATLAB:

dataDrivenRearrangement(1000, 20, 'BARON', 'CPLEX', 'CPLEX', ...
    'C:\\Models\\Step1_Grouping.mod', ...
    'C:\\Models\\Step2_SingleShelf.mod', ...
    'C:\\Models\\Step2_TwoShelf.mod', ...
    'C:\\Models\\Step3_AisleAssign.mod');

Ensure the AMPL API is initialized before running:

>> setupOnce

## Notes

- All model input is passed dynamically via MATLAB; do not use .dat files.
- Each model is solver-agnostic but tuned for BARON (Step 1) and CPLEX (Steps 2 and 3).
- Outputs from AMPL are processed in MATLAB and stored in global structures for further analysis.
