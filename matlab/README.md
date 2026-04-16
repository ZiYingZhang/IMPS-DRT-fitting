# MATLAB workflow for IMPS fitting

This folder contains a MATLAB-based workflow for IMPS fitting and DRT-related analysis.

## Current contents

### MATLAB script
- `IMPS_elaticnet_beta_en.m`  
  Main MATLAB script for IMPS fitting and regularized analysis.

### Demo input data
- `0.1 V.txt`
- `0.2 V.txt`
- `0.3 V.txt`
- `0.4 V.txt`

These files are example IMPS datasets measured at different applied potentials.

### Output files
For each input dataset, the script generates three output files:

- `*_DRT.txt`  
  DRT-related output data

- `*_H.txt`  
  Intermediate or processed fitting-related output

- `*_JV.txt`  
  Additional output data exported by the script

For example, for `0.1 V.txt`, the corresponding outputs are:

- `0.1 V.txt_DRT.txt`
- `0.1 V.txt_H.txt`
- `0.1 V.txt_JV.txt`

## Purpose

This folder is intended to provide:

- a MATLAB implementation for IMPS fitting
- demo datasets for testing
- exported result files for comparison and verification
- a starting point for further development and repository organization

## How to use

1. Open MATLAB
2. Set the current working directory to this `matlab/` folder
3. Open and run `IMPS_elaticnet_beta_en.m`
4. Use one of the provided demo datasets as input
5. Check the generated output files after the script finishes

## Notes

At the current stage, all demo data, script files, and output files are stored directly in this folder for convenience.


## Requirements

Before running the MATLAB script, please make sure that **CVX** has been installed.

CVX can be downloaded from the official website:

- https://cvxr.com/cvx/download/

After downloading and extracting the package, open MATLAB, switch to the CVX folder, and run:

```matlab
cvx_setup
