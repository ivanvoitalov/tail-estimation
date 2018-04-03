# Tail Index Estimation for Degree Sequences of Complex Networks

## Script Summary

This script provides a set of tools to estimate tail index of degree distribution of complex networks. Calculated estimates can be then plotted along with PDF/CCDF of degree distribution. 

## Dependencies

This script requires several packages to be installed:
* NumPy
* Matplotlib

## Installation

This script does not require installation. Simply download 'tail_estimation.py' for your Python version (2 or 3) from the corresponding directory at this repository.

## Simple Usage Example

We provide the simplest usage example.

## Command Line Options

## Implemented Estimators

Currently several classical estimators for tail index estimation are implemented:
* Hill estimator, including smoothed Hill estimator and adjusted Hill estimator;
* moments estimator;
* Pickands estimator;
* kernel-type estimator.

## Double-bootstrap for Optimal Threshold Estimation

The package implement double-bootstrap estimation of the optimal order statistic for several estimators:
* Hill;
* Moments;
* Kernel-type.

Hill double-boostrap: [Danielsson et al. (2001)](https://www.riskresearch.org/papers/DanielssonHaanPengVries2001/)
Moments double-bootstrap: [Draisma et al. (1999)](https://link.springer.com/article/10.1023/A:1009900215680)
Kernel-type double-bootstrap: [Groeneboom et al. (2003)](https://www.jstor.org/stable/3448443)

## Usage Example

## Additional Info
