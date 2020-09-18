# AutoBQP
A generic heuristic solver for binary optimization problems!

This solver implements a set of heuristic strategies for the unconstrained Binary Quadratic Programming (BQP), and applies [irace](http://iridia.ulb.ac.be/irace) to search for the best combination of components and parameter values for solving the problem at hand.

**Maintainer:** [Marcelo de Souza](https://souzamarcelo.github.io).

**Contributors:** [Marcus Ritt](https://www.inf.ufrgs.br/~mrpritt).

If you have any difficult or want to collaborate with us, please write to me: marcelo.desouza@udesc.br.

## Relevant Literature

The following article describes in detail the AutoBQP solver and presents applications on solving BQP and MaxCut problems:

+ Marcelo de Souza, Marcus Ritt. **Automatic Grammar-Based Design of Heuristic Algorithms for Unconstrained Binary Quadratic Programming**. In Liefooghe A., López-Ibáñez M. (Org.), Lecture Notes in Computer Science, 1ed.: Springer International Publishing, v. 10782, p. 67-84, 2018.

The following article applies the AutoBQP solver to the test assignment problem:

+ Marcelo de Souza, Marcus Ritt. **An Automatically Designed Recombination Heuristic for the Test-Assignment Problem**. 2018 IEEE Congress on Evolutionary Computation (CEC), Rio de Janeiro, 2018.

## Examples

AutoBQP was previously used to solve the following problems:

* BQP (Binary Quadratic Programming)
* MaxCut (Maximum Cut in graphs)
* Test Assignment Problem

File *instance.hpp* provides the functions to read the instances of such problems, as well as the reduction techniques for the test assignment problem.

AutoBQP are being used to explore other problems from literature, including:

* [Maximum Weight Vertex Clique](https://link.springer.com/article/10.1007/s10878-016-9990-2)
* [Maximum Diversity Problem](https://www.sciencedirect.com/science/article/abs/pii/S0377221706000634)
