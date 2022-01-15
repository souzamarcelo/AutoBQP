# AutoBQP

**A Component-Wise Solver to Binary Optimization!**

AutoBQP is a generic black-box solver for binary problems based on the automatic design of heuristic algorithms. This solver implements a set of heuristic strategies for the unconstrained Binary Quadratic Programming (BQP), and applies [irace](http://iridia.ulb.ac.be/irace) to search for the best combination of components and parameter values for solving the problem at hand.

**Maintainer:** [Marcelo de Souza](https://souzamarcelo.github.io).

**Contributors:** [Marcus Ritt](https://www.inf.ufrgs.br/~mrpritt).

If you have any difficult or want to collaborate with us, please write to me: marcelo.desouza@udesc.br.

***

## Relevant Literature

The following article describes in detail the AutoBQP solver and presents applications on solving **BQP** and **MaxCut** problems:

+ Marcelo de Souza, Marcus Ritt. **Automatic Grammar-Based Design of Heuristic Algorithms for Unconstrained Binary Quadratic Programming**. In Liefooghe A., López-Ibáñez M. (Org.), Lecture Notes in Computer Science, 1ed.: Springer International Publishing, v. 10782, p. 67-84, 2018.

The following article applies the AutoBQP solver to the **test assignment** problem:

+ Marcelo de Souza, Marcus Ritt. **An Automatically Designed Recombination Heuristic for the Test-Assignment Problem**. 2018 IEEE Congress on Evolutionary Computation (CEC), Rio de Janeiro, 2018.

Please, make sure to reference us if you use AutoBQP in your research.

```bibtex
@inproceedings{SouzaEtAl2018autobqp,
  title        = {Automatic Grammar-Based Design of Heuristic Algorithms for Unconstrained Binary Quadratic Programming},
  author       = {Souza, Marcelo and Ritt, Marcus},
  booktitle    = {Evolutionary Computation in Combinatorial Optimization},
  year         = {2018},
  publisher    = {Springer},
  pages        = {67--84}
}
```

***

## Requirements

AutoBQP is implemented in C++/Python and requires:
+ [GCC/G++](https://gcc.gnu.org) compiler
+ [CMake](https://cmake.org) and [GNU Make](https://www.gnu.org/software/make)
+ [Boost C++ libraries](https://www.boost.org)
+ [Python 3](https://www.python.org)
+ [R](https://www.r-project.org) and [irace](https://iridia.ulb.ac.be/irace)

***

## Usage

To use AutoBQP, follow the steps below.

1. [Download](https://github.com/souzamarcelo/AutoBQP/archive/master.zip) and unzip AutoBQP.
2. Create a directory with the instances of the problem being solved: `instances`.
   + See this example: [aac/instances/](aac/instances)
3. Create a text file with each instance (file name) and arguments: `instances.txt`.
   + See this example: [aac/instances.txt](aac/instances.txt)
   + For example, in line `p3000-1.def --maximize 1 --format pgen --runtime 10`:
     + `p3000-1.def` is the instance file name;
     + `--maximize 1` indicates a maximization problem;
     + `--format pgen` indicates the instance format;
     + `--runtime 10` defines the running time limit in seconds.
4. If needed, provide a function to read the instance in the [src/instance.hpp](src/instance.hpp) file.
   + Observe that the code calls the corresponding function according to the instance format.
5. You can also make problem-specific modifications in files [src/solution.hpp](src/solution.hpp) and [src/runner.cpp](src/runner.cpp), if needed.
6. Run [AutoBQP.py](AutoBQP.py).
   + `./AutoBQP.py --insdir /path/to/instances --insfile /path/to/instances.txt --budget 2000 --parallel 4`.
     + Option `--budget` defines the maximum number of executions for the automatic design process.
     + Option `--parallel` defines the number of cores to be used for parallelization.

You can also just execute an algorithm to solve an instance. For example, given the algorithm described by the corresponding parameter values:

```
--cAlg 2 --cSearch 3 --cSearch 2 --cModification 3 --cModification 6 --cTs 1 --cPert 1 --cStep 3 --pTenureType r --pTenureVariation 71 --pCm 359 --pMaxStagnateType a --pMum 17 --pIc 73 --pMaxStepsType n --pMaxStepsValue 44423 --pP 0.06 --pD1 29 --pD2 72 --pB 13 --pGamma 0.38 --pGammam 15 --pBsize 3 --pSiPartialSize 33
```

To run that algorithm, first compile the grammar implementation:

```bash
cd src
cmake .
make grammar
```

Then run the grammar to generate the source code of the algorithm:

```bash
./grammar --cAlg 2 --cSearch 3 --cSearch 2 --cModification 3 --cModification 6 --cTs 1 --cPert 1 --cStep 3 --pTenureType r --pTenureVariation 71 --pCm 359 --pMaxStagnateType a --pMum 17 --pIc 73 --pMaxStepsType n --pMaxStepsValue 44423 --pP 0.06 --pD1 29 --pD2 72 --pB 13 --pGamma 0.38 --pGammam 15 --pBsize 3 --pSiPartialSize 33 > algorithm.cpp
```

Finally, compile the runner and execute the algorithm (in this case, to solve instance `p3000-1.def`):

```bash
make runner
./runner --ins /path/to/p3000-1.def --maximize 1 --format pgen --runtime 10
```
