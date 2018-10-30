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

One-dimensional spatial structure is controlled by left and right parameters L\_left and L\_right where the probability of node _j_ connecting to node _i_ is approximately Gaussian with standard deviation L\_left if _j_>_i_ and L\_right if _j_<_i_ .
<!--((_i_-_j_)^2)/2(L\_left)^2 or ((_i_-_j_)^2)/2(L\_right)^2, whichever is positive-->

For homogenous spatial structure as in the [SONET project](https://github.com/dqnykamp/sonets "SONET Github") by Zhao et al., let L\_left = L\_right = inf.

Reference: L. Zhao, B. Beverlin II, T Netoff, and D. Q. Nykamp.  Synchronization from second order network connectivity statistics.
*Frontiers in Computational Neuroscience*, 5:28, 2011.  [Publisher's web site](http://dx.doi.org/10.1007/s10827-011-0373-5).

### Simulation and Event Analysis
Python code that utilizes Brian to simulate on networks of excitatory neurons with instantaneous synapses using a leaky integrate\-and\-fire model.  Analyzes simulation results for traveling waves of neuron firings using an original synchronous event detection algorithm.

## Requirements
* Python 3.6 or later
* cython (Python package that utilizes C)
* C compiler
* [Brian](https://brian2.readthedocs.io/en/stable/index.html) simulator for spiking neural networks (not needed for network generation)

## Usage Notes
### Generating Networks
Scripts needed:
* create_P.py
* useP_to_createW_withC.py
* produceW.pyx \- must be compiled into produceW.c (use the command below)
* generate_Ws.py \- use this to create multiple adjacency matrices with same parameters (optional)

Usage details:  
To generate an adjacency matrix for a single network in a Python program, first process the produceW.pyx cython code to create the produceW.c file (this file usually needs to be generated on each machine where it will be used).  To do this, try 
```python setup.py build_ext --inplace```
in the command line.

Next, use the `create_P` function in create_P.py with
```python
P = create_P(N, L_left, L_right, p_AVG)
```
where N is the number of nodes in the network, L\_left and L\_right are as described in the Network Generation Overview (above), and p\_AVG is the average number of connections per neuron in the network. P is an NxN matrix representing probabilities of directed connections between neurons so P_ij = P(W_ij = 1) for any nodes i != j.

Then, to generate a single adjacency matrix use the `create_W` function in useP_to_createW_withC.py 
```python
W = create_W(N, P, alpha_recip, alpha_conv, alpha_div, alpha_chain)
```
W is the adjacency matrix for the network where W_ij=1 if there is a connection from node j to node i and is zero otherwise.  The parameter P is the probability matrix from create_P (above).  This create_W function operates under the assumption that the sum of each column of P is the same (i.e., each neuron has the same out-degree distribution), which our create_P function satisfies.

The alphas specify the probability of a motif W_ij and W_kl according to P(W_ij=1 and W_kl =1) = P_ij\*P_kl(1+alpha) where alpha reperesent the appropriate motif parameter, as follows.  If i, j, and k represent distinct nodes, then alpha=alpha_recip for the motifs W_ij and W_ji, alpha=alpha_conv for the motifs are W_ij and W_ik, alpha=alpha_div for the motifs W_ij and W_kj, and alpha=alpha_chain for the motifs W_ij and W_jk, as described below.  The alphas can be interpreted as specifying the motif frequency or the covariance between edges in the motif.  

Only a small range of parameters alpha_recip, alpha_conv, alpha_div, and alpha_chain represent actual networks.  The function create_W raise an exception if it was not able to generate a network.  Given their definition, the range of the alphas might appear that they could range approximately between -1 and (1/p_AVG)-1.  However, in reality, their range is much smaller.  Further, if individual entries in matrix P are too large, then there will be less flexibility among the alpha values.

The parameters alpha_conv and alpha_div must be non-negative, as they determine the variance of the in- and out-degree distributions.  However, when alpha_conv and alpha_div are close to their maximum values, possible combinations of alpha_recip and alpha_chain are highly restricted.  

The parameter alpha_recip can be negative.  If alpha_recip is close to its maximum value, then network is nearly undirected, which requires alpha_conv, alpha_div, and alpha_chain to be approximately equivalent.  Keeping alpha_recip much smaller allows more flexibility in the other parameters. In fact, the algorithm is designed to ignore the alpha_recip parameter to allow for more control and attention to the other alpha parameters.

The parameter alpha_chain determines the covariance between the in- and out-degrees.  The variance of the in- and out-degrees are determined by alpha_conv and alpha_div, respectively, so the limits of alpha_chain depend on alpha_conv and alpha_div.

Multiple networks with the same parameters can be generated using the generate_Ws.py script, which calls the other scripts.

### Simulations and Event Analysis
