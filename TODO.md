# To-Do: 

## Refactoring / Optimizing 

- [x] Remove calls to Arc struct in MIP construction 
- [ ] Remove calls to Arc struct in Benders: 
* used to obtain incoming arcs in flow constraints, either needs: 
    * reverse adjacency list
    * adjacency matrix 
    * dijkstra's 
