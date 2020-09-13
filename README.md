# Multi Choice Knapsack problem 
The Multi Choice Knapsack problem (referred to as MCKP) is a generalisation of the ordinaty knapsack problem, where the set of items in partionned into different classes. The problem consists of finding exactly one item per class in order to maximize the total gain while respecting a total capacity constraint.

Multiple approaches are possible, both to the exact and relaxed problem. This repository contains an implementation of different approches, and compares performance.

This repository contains a multichoice knapsack problem solver, using various techniques. We also use a concrete example of multi-channel communication to illustrate the problem and the solver.

For a full description of the problem : 

## Main class
This class uses the solver/MCKP_Solver class in order to solve the multi-channel communication problem. It offers a way to test the solver and to compare performance of the different approaches. We give 5 test files in the folder testfiles. These are the files we used to test and compare the solver.
The class can be used in two ways :

#### javac Main {filepath}
Solves the problem in the file using all approaches and display the results

Example : javac Main testfiles/test1.txt
```
Maximum rate found by LP_solver is : 365.0
Maximum rate found by greedy is : 365.0
Maximum rate found by DP1 is : 365
Maximum rate found by DP2 is : 365
Maximum rate found by BB is : 365
```

#### javac Main -r {filepath}
Compares runtime of the different approaches solving the problem in file.

Example : javac Main -r testfiles/test1.txt
```
Algorithm | Runtime (ms) | Result
dp1 | 0.1165425 | 365
dp2 | 0.19613430999999998 | 365
bb | 0.04078498 | 365
```


## MCKP Solver
The main class is solver/MCKP_Solver. Each MCKP_Solver object contains problem data, and implements different methods to preprocess the data and solve the problem.

### Preprocessing methods :

#### remove_impossible_terms() : 
Removes impossible terms (terms that can be part of any solution, otherwise the capacity constraint is violated). This method throws an exception if the problem is impossible to solve, ie, if even by taking the smallest weights in every class, the capacity constraint can't be respected.

#### remove_IP_dominated():
Removes IP dominated terms (in place). See document below for full explanation

#### remove_LP_dominated():
Removed LP dominated terms. This is done by inserted non-dominated terms in the field  LP_filtered_data. This is not done in place because this prepocessing step is only usefull when solving the relaxed LP problem.

### Solver methods :

#### public Solution greedy_LP() :
Applies the greedy approach to solve the LP-relaxed problem. Returns a Solution object, containing the rate achieved and an array of terms chosen in each class.

#### public double LP_solver() :
Uses scpsolver library to solve the LP-relaxed problem. The library used the simplex method to solve the problem.

#### public int DP_1() :
Uses a first dynammic programming approach to solve the exact MCKP problem.

#### public int DP_2(int U) :
Uses a second dynammic programming approach. The method needs an upper bound for the maximum rate achieved. An upper bound can be computed by upper_bnd_Rate() method.

#### public int Braunch_and_bound() :
Solves the problem using a branch and bound approach, recursively.

#### public int BB_with_queue() :
A branch and bound solver using a queue.

### Visualise data :
You can visualise data after each step of the preprecessing. 

#### public void visualize_data(int _class, String additionnal_title, boolean lp_filtred)
Plots terms in _class (profit, weight). Set lp_filtred to true to visualise data after without the LP-dominated terms, otherwise the method plots whatever is in data field (Depending on when the method is called, in can plot raw or filtered data).




