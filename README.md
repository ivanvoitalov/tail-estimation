# Tail Index Estimation for Degree Sequences of Complex Networks

[summary](#script-summary)

This script is intended to be a single-file simple solution to the tail index estimation problem in complex networks. It consists of several well-established estimators combined into one toolbox along with some useful plotting routines that usually help to analyze a given degree distribution.

## Dependencies

The script requires the following packages to be installed in addition to the standard libraries shipped with Python:
* NumPy
* Matplotlib

## Installation

No installation is required. Simply download _tail-estimation.py_ file from the directory corresponding to your Python version (either 2 or 3) and use it as demonstrated in the Examples section below.

## Required Input Format

The script processes a degree sequence in the form of:
```
k n(k)
```
where **k** is a node's degree and **n(k)** is the number of nodes with such degree in the network.

## Simple Usage Example

Here we provide the simplest usage example of the _tail-estimation_ script. 
Suppose we want to compute the tail index of a degree distribution of the [CAIDA network](http://konect.uni-koblenz.de/networks/as-caida20071105) provided by [KONECT database](http://konect.uni-koblenz.de/). We first convert the network in the form of an edge list to the format of _degree counts_ as indicated above. The converted data can be found under the _Examples_ directory under the **CAIDA_KONECT.dat** name. Then we run the _tail-estimation.py_ script as follows:
```
python tail-estimation.py /path/to/degree/sequence /path/to/output/plots/file
```

This will produce an image file with plots as well as some STDOUT messages reporting estimated tail indices. An example of such image for the CAIDA network is given below:
![CAIDA Output](https://raw.githubusercontent.com/ivanvoitalov/tail-estimation/master/Figures/CAIDA_output.png)

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
