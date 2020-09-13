# Multi Choice Knapsack problem 
The Multi Choice Knapsack problem (referred to as MCKP) is a generalisation of the ordinaty knapsack problem, where the set of items in partionned into different classes. The problem consists of finding exactly one item per class in order to maximize the total gain while respecting a total capacity constraint.

Multiple approaches are possible, both to the exact and relaxed problem. This repository contains an implementation of different approches, and compares performance.

This repository contains a multichoice knapsack problem solver, using various techniques. We also use a concrete example of multi-channel communication to illustrate the problem and the solver.

## MCKP Solver
The main class is solver/MCKP_Solver. Each MCKP_Solver object contains problem data, and implements different methods to preprocess the data and solve the problem.

### Preprocessing methods :

#### remove_impossible_terms() : 
Removes impossible terms (terms that can be part of any solution, otherwise the capacity constraint is violated). This method throws an exception if the problem is impossible to solve, ie, if even by taking the smallest weights in every channel, the capacity constraint can't be respected.

#### remove_IP_dominated():
Removes IP dominated terms (in place). See document below for full explanation

#### remove_LP_dominated():
Removed LP dominated terms. This is done by inserted non-dominated terms in the field  LP_filtered_data. This is not done in place because this prepocessing step is only usefull when solving the relaxed LP problem.

