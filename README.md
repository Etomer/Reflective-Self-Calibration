# Reflective-Self-Calibration

## Introduction
This repository contains code for the paper "Sensor node calibration in presence of a dominant reflective plane" at EUSIPCO 2022

## Abstract
Recent advances in simultaneous estimation of both receiver and sender positions in ad-hoc sensor networks have made it possible to automatically calibrate node positions – a prerequisite for many applications. In man-made environments there are often large planar reflective surfaces that give significant reverberations. In this paper, we study geometric problems of receiver-sender node calibration in the presence of such reflective planes. We establish a rank-1 factorization problem that can be used to simplify the estimation. We also show how to estimate offsets, in the Time difference of arrival case, using only the rank constraint. Finally, we present a new solver for the minimal cases of sender-receiver position estimation. These contributions result in a powerful stratified approach for the node calibration problem, given a reflective plane. The methods are verified with both synthetic and real data.

## Reproducting Results
To run the solvers in *src/solvers/* you might have to (depending on your system) mex the ".cpp" files. To do this, you need to install [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page). You can then run the following in MATLAB
```
mex(['-I' cg_eigen_dir],'solver_toa_34_mirror_r_dist.cpp')
```
where `cg_eigen_dir` is the path to where Eigen is installed.
