Author: Zhiqiang Xu, 10/2014
Email: zxu1@e.ntu.edu.sg

This package contains the MATLAB implementation of the NCG-VB algorithm for SBM inference, proposed in the paper "A Fast Inference Algorithm for Stochastic Blockmodel" (Xu et al., 2012).

Note that this code is only applicable to undirected and unweighted graphs. 

The main files include:

1. 5 datasets: ca-GrQc, ca-HepTh, PGP, twitter, wordnet
2. initial values for natural parameters: they are used in the experiments reported in the paper (note that only
those for the ca-GrQr dataset are included. contact with me for the others)
3. main.m: demonstrate how to use our proposed algorithm
4. init.m: random initialization of natural parameters
5. ncg_sbm.m: the proposed algorithm
6. cg_sbm.m: CG-VB for SBM inference as a baseline
7. evaluate_quality.m: calculate modularity and conductance
8. lower_bound.m: calculate the lower bound function
9. transform.m: transform natural parameters to expectation parameters