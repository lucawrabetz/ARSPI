# To-Do: 

## Refactoring

- [ ] Use copy constructors for gurobi var vectors (already done in line 313 of M2.cpp with x, check it compiles when solving MIP then do with the rest of them).
- [ ] Change the copy constructor for M2ProblemInstance to take a const object, and verify the consistency of that usage throughout libraries and scripts - 
- [ ] Refactor main into scripts that can be called by python or shell scripts for computational experiments

## Optimizing 

- [x] Remove calls to Arc struct in MIP construction 
- [ ] Remove calls to Arc struct in Benders: 
* used to obtain incoming arcs in flow constraints, either needs: 
    * reverse adjacency list
    * adjacency matrix 
    * dijkstra's 
- [ ] Better implementation of Dijkstra (priority queue) if it will be used 
- [ ] return values of solver calls - for dense graphs the interdiction policy is memory intensive, consider revising all of it for a sparse representation


## Low Priority

- [ ] Integer bitwise representation and updates of a graph's available arc set
- [ ] Bitwise enumeration code 

## Locking down naming
* CapitalCase for all types - classes, structs, enums, etc
* lower_case_with_underscores or lowercase for variables
* ALLCAPS for consts
* mixedCaseWithLeadingLowerCase() for functions
* same_as_variables for accessors and mutators
* File naming (to be defined when we start writing and defining the scripts for experiments :) ): 

## Scripts
-


