# Construction of ellipsoid convex model with outlier detection and its application in NRBTO

Thank you for visiting the repo for the current research!

## Prerequisites

1. MATLAB(R) of version 2017b or higher is required.

2. The MATLAB(R) version of the CVX toolbox ought to be obtained via the link:
<http://cvxr.com/cvx/download/>.

3. The MATLAB(R) code of the MMA optimizer should be obtained following the instructions published on:
<https://people.kth.se/~krille/>.

After the MMA codes are obtained, please modify Line 195 in "det_lof_2d.m" in folder "NRBTO" with the main function name of the MMA optimizer.

Also, please modify Line 222 in "nrbto_lof_2d.m" in folder "NRBTO" with the main function name of the MMA optimizer.

## Construction of the ellipsoid convex model

### The codes for this module are in the folder "construct_ECM"

The input parameters are included in the file "find_outlier_settings.json", where

1. "sample_filename": The Microsoft(R) Excel(R) file that contains the samples data.

2. "normalization": The flag of the normalization of the samples, "true" for normalize, "false" for not normalize.

Execute the main function "find_outlier.m" to construct the ellipsoid convex model with outliers detection. The results are printed to the command window, as well as recorded in the MAT-file.

## NRBTO

### The codes for this module are in the folder "nrbto_lof"

In this folder,

1. "deterministic_cantilever_beam.jason" contains the parameters for the deterministic topology optimization. The deterministic optimization is executed using the function "det_lof_2d.m".

2. The parameters for the batch jobs of the NRBTO are stored in the sequential files "NRBTO_cantilever_beam_##.m", where "##" is the sequence number starting from "01". The NRBTO batch jobs are executed using the function "NRBTO_batch_job.m".
