# Adaptive Network Interdiction

This repository contains the source code to solve the adaptive network interdiction problem (adaptive shortest-path interdiction). The code to generate test instances and handle experiments can be found in [scripts](https://github.com/lucawrabetz/ARSPI/tree/master/scripts).

## Key Classes

### LayerGraph

* The LayerGraph class, defined in src/M2.cpp, is used to represent a directed network, constructed from an arc-list text file. The constructor creates the following data structures: 
    * A vector of Arc objects - indexed 0-(m-1) in the same order as the lines of the text file - the Arc objects simply two member ints - i, and j (want to remove - not used to build MIP anymore but is used for flow constraints in benders - either use reverse adjacency list or use Dijkstra's as subroutine)
    * A vector of vectors of ints, outer vector indexed 0-(n-1), of adjacency list for each node. Each inner vector is a list of adjacent nodes;

### AdaptiveInstance
