# Construction of ellipsoid convex model with outlier detection and its application in Non-Probabilistic Topology Optimization

Thank you for visiting the repo for the current research!

## Prerequisites

1. MATLAB(R) of version 2017b or higher is required.

2. The MATLAB(R) version of the CVX toolbox ought to be obtained via the link:
<https://cvxr.com/cvx/download/>.

3. The MATLAB(R) code of the MMA optimizer should be obtained following the instructions published on:
<https://people.kth.se/~krille/>.

Please install the CVX toolbox in MATLAB(R) following the directions provided therein.

After the MMA codes are obtained, please modify Line 195 in "det_lof_2d.m" in folder "NRTO" with the main function name of the MMA optimizer (usually left unchanged as "mmasub" therein).

Also, please modify Line 222 in "nrto_lof_2d.m" in folder "NRTO" with the main function name of the MMA optimizer (usually left unchanged as "mmasub" therein).

## Construction of the ellipsoid convex model

### The codes for this module are in the folder "construct_ECM"

The input parameters are included in the file "find_outlier_settings.json", where

1. "sample_filename": The Microsoft(R) Excel(R) file that contains the samples data.

2. "normalization": The flag of the normalization of the samples, "true" for normalize, "false" for not normalize.

Execute the main function "find_outlier.m" to construct the ellipsoid convex model with outliers detection. The results are printed to the command window, as well as recorded in the MAT-file.

## NRTO

### The codes for this module are in the folder "NRTO"

In this folder,

1. "deterministic_cantilever_beam.jason" contains the parameters for the deterministic topology optimization. The deterministic optimization is executed using the function "det_lof_2d.m".

2. The parameters for the batch jobs of the NRTO are stored in the sequential files "NRTO_cantilever_beam_##.m", where "##" is the sequence number starting from "01". The NRTO batch jobs are executed using the function "NRTO_batch_job.m".

## Detection of outliers and construction of ellipsoid convex model

### The 3D example

![Three dimentional example](./image/outlier_ECM_3D.jpg)

## Application in Non-probabilistic Robust Topology optimization

### Deterministic solution's evolution

![Deterministic solution](./image/deterministic.gif)

### Evolution of robust solution based on the normal samples

![Robust solution based on normal samples](./image/NRTO_normal.gif)

### Evolution of robust solution based on all samples

![Robust solution based on all samples](./image/NRTO_all.gif)
