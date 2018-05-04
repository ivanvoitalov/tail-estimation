# Tail Index Estimation for Degree Sequences of Complex Networks

## Summary

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
where **k** is a node's degree and **n(k)** is the number of nodes with such degree in the network. Note that these two number are whitespace-separated.

## Simple Usage Example

Here we provide the simplest usage example of the _tail-estimation_ script. 
Suppose we want to compute the tail index of a degree distribution of the [CAIDA network](http://konect.uni-koblenz.de/networks/as-caida20071105) provided by [KONECT database](http://konect.uni-koblenz.de/). We first convert the network from the format of an edge list to the format of _degree counts_ as indicated above. The converted data can be found under the _Examples_ directory under the **CAIDA_KONECT.dat** name. Then we run the _tail-estimation.py_ script as follows:
```
python tail-estimation.py /path/to/degree/sequence /path/to/output/plots/file
```

This will produce a collection of plots as well as some STDOUT messages reporting estimated tail indices. Most users would be interested in just one-number tail index estimates according to three estimators we have implemented so far. They are reported to the STDOUT in the following form (for the CAIDA network example):
```
**********
Adjusted Hill estimated gamma: 2.09313899261
**********
Moments estimated gamma: 2.11325330189
**********
Kernel-type estimated gamma: 2.13032486828
**********
```
An example of plots generated for the CAIDA network is given below:
![CAIDA Output](https://raw.githubusercontent.com/ivanvoitalov/tail-estimation/master/Figures/CAIDA_output.png)

### What does it mean?!

Although at the first glance the plots produced by the script may seem to be complicated, it is very easy to interpret them for your network! The main thing to notice is that all tail index estimates are plotted in terms of parameter <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi$" title="$\xi$" /></a> that is related to the tail index of the PDF of degree distribution <a href="https://www.codecogs.com/eqnedit.php?latex=$\gamma$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\gamma$" title="$\gamma$" /></a> as follows: <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi=&space;1/(\gamma&space;-&space;1)$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi=&space;1/(\gamma&space;-&space;1)$" title="$\xi= 1/(\gamma - 1)$" /></a>. Here we list the description of what is exactly shown on each subfigure, starting from the top left one:
1. Log-binned probability density function (PDF) of a given degree sequence.
2. Complementary cumulative distribution function (CCDF) of a given degree sequence.
3. Smooth and adjusted Hill estimates of <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi=&space;1/(\gamma&space;-&space;1)$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi=&space;1/(\gamma&space;-&space;1)$" title="$\xi= 1/(\gamma - 1)$" /></a> parameter on a linear scale as a function of the number of included degree sequence order statistic <a href="https://www.codecogs.com/eqnedit.php?latex=$\kappa$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\kappa$" title="$\kappa$" /></a>. Star marker shows best estimate of the tail index <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi$" title="$\xi$" /></a> based on the double-bootstrap procedure for Hill estimator.
4. Smooth and adjusted Hill estimates of <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi$" title="$\xi$" /></a> parameter on a semilog scale as a function of the number of included degree sequence order statistic <a href="https://www.codecogs.com/eqnedit.php?latex=$\kappa$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\kappa$" title="$\kappa$" /></a>. Star marker shows best estimate of the tail index <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi$" title="$\xi$" /></a> based on the double-bootstrap procedure for Hill estimator.
5. Pickands, moments and kernel-type estimates of <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi$" title="$\xi$" /></a> parameter on a linear scale as a function of the number of included degree sequence order statistic <a href="https://www.codecogs.com/eqnedit.php?latex=$\kappa$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\kappa$" title="$\kappa$" /></a>. Star markers show best estimate of the tail index <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi$" title="$\xi$" /></a> based on the double-bootstrap procedure for moments and kernel-type estimators.
6. Pickands, moments and kernel-type estimates of <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi$" title="$\xi$" /></a> parameter on a semilog scale as a function of the number of included degree sequence order statistic <a href="https://www.codecogs.com/eqnedit.php?latex=$\kappa$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\kappa$" title="$\kappa$" /></a>. Star markers show best estimate of the tail index <a href="https://www.codecogs.com/eqnedit.php?latex=$\xi$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\xi$" title="$\xi$" /></a> based on the double-bootstrap procedure for moments and kernel-type estimators.

## Advanced Examples

### Generating diagnostic plots

### Changing _epsstop_ parameter

### Analyzing non-network data

### Noise and why is it added

## Implemented Estimators

Currently several classical estimators for tail index estimation are implemented:
* Hill estimator, including smoothed Hill estimator;
* moments estimator;
* Pickands estimator;
* kernel-type estimator.

## Double-bootstrap for Optimal Threshold Estimation

The package implements double-bootstrap estimation of the optimal order statistic for several estimators:
* Hill;
* Moments;
* Kernel-type.

Hill double-boostrap: [Danielsson et al. (2001)](https://www.riskresearch.org/papers/DanielssonHaanPengVries2001/)
Moments double-bootstrap: [Draisma et al. (1999)](https://link.springer.com/article/10.1023/A:1009900215680)
Kernel-type double-bootstrap: [Groeneboom et al. (2003)](https://www.jstor.org/stable/3448443)

## Command Line Options
