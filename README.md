# spatialSONET
Generates Second Order Networks (SONETs) with prescribed second order motif frequencies and double-sided spatial factor.  Simulates the spiking neuron networks with Brian2 simulator.  Analyzes simulation results with synchronous event detection.

## Overview
### Network Generation
Python code that generates random directed graphs with prescribed second order statistcs and one-dimensional spatial structure.  Additional structure (i.e., higher-order statistics) are approximately minimized.

Second order statistics are prescribed based on frequencies of four types of two-edge motifs (subnetworks containing two edges):
* reciprocal motifs (two nodes reciprocally connected:  o <-> o)
* convergent motifs (two nodes each connecting onto a third node:  o -> o <- o)
* divergent motifs (two nodes each receving a connection from a third node: o <- o -> o)
* chain motifs (three nodes connected in sequence: o -> o -> o)

One-dimensional spatial structure is controlled by left and right parameters L\_left and L\_right where the probability of node _j_ connecting to node _i_ is approximately Gaussian with mean standard deviation L\_left if _j_>_i_ and L\_right if _j_<_i_ .
<!--((_i_-_j_)^2)/2(L\_left)^2 or ((_i_-_j_)^2)/2(L\_right)^2, whichever is positive-->

For homogenous spatial structure as in the [SONET project](https://github.com/dqnykamp/sonets "SONET Github") by Zhao et al., let L\_left = L\_right = inf.
Reference: L. Zhao, B. Beverlin II, T Netoff, and D. Q. Nykamp.  Synchronization from second order network connectivity statistics.
*Frontiers in Computational Neuroscience*, 5:28, 2011.  [Publisher's web site](http://dx.doi.org/10.1007/s10827-011-0373-5).

### Simulation and Event Analysis

## Requirements
* Python 3.6 or later
* [Brian 2](https://brian2.readthedocs.io/en/stable/index.html) simulator for spiking neural networks (not needed for network generation)

## Usage Notes
