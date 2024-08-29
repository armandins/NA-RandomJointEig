# RandomJointEig
Numerical examples for Randomized Joint Eigenvalue Approximation

This repository contains numerical examples for the paper:  Haoze He, Daniel Kressner, 
  Bor Plestenjak, "Randomized methods for computing joint eigenvalues, 
  with applications to multiparameter eigenvalue problems and root finding"
  
Matlab examples run in Matlab 2023a or higher. 

## Numerical experiments for synthetic data

To reproduce numerical results in this section, set variables example, deltadif (only for example=3), Nsamples and run ExampleMain

Possible values of example are:

    - example=1.1 : Left Figure 2 from Example 5.1 (default example)
    - example=1.2 : Right Figure 2 from Example 5.1 
    - example=1.3 : Figure 3 from Example 5.1 (set deltadif to 1e-2, 1e-4, ..., 1e-16)
    - example=1.4 : Right Figure 4 from Example 5.2 (empirical failure of 1-S RQ)
    - example=2.1 : Left Figure 5 from Example 5.3
    - example=2.2 : Right Figure 5 from Example 5.3
    - example=3   : Example 5.4, set also deltadif
    - example=4.1 : Left Figure 8 from Example 5.6
    - example=4.2 : Right Figure 8 from Example 5.6
 
To change the sample size (default is 10000), set variable Nsamples

Run ExampleLemma44 to get Figure 3 (left)

Setting example=1.2 requires Advanpix Multiprecision Computing Toolbox. If not available, change line 45 in testRayleigh.m to use_advapix=0  

## Multiparameter eigenvalue problems

Run ExampleEisenmann to reproduce Example 6.1

Run ExampleBrassTeflon to reproduce Example 6.2

## Roots of polynomial systems

Run ExampleGrafTownsend to reproduce Example 7.2

