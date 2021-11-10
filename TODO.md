# To-Do: 

## Refactoring

- [ ] Use copy constructors for gurobi var vectors (already done in line 313 of M2.cpp with x, check it compiles when solving MIP then do with the rest of them).

## Optimizing 

- [x] Remove calls to Arc struct in MIP construction 
- [ ] Remove calls to Arc struct in Benders: 
* used to obtain incoming arcs in flow constraints, either needs: 
    * reverse adjacency list
    * adjacency matrix 
    * dijkstra's 
- [ ] Better implementation of Dijkstra (priority queue) if it will be used 


## Low Priority

- [ ] Integer bitwise representation and updates of a graph's available arc set
- [ ] Bitwise enumeration code 
