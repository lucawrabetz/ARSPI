# Adaptive Network Interdiction

This repository contains the source code to solve the adaptive network interdiction problem (adaptive shortest-path interdiction). The code to generate test instances and handle experiments can be found in [scripts](https://github.com/lucawrabetz/ARSPI/tree/master/scripts).

## Classes

### LayerGraph

* The LayerGraph class, defined in src/M2.cpp, is used to represent a directed network, constructed from an arc-list text file. The constructor creates the following data structures: 
    * A vector of Arc objects - indexed 0-(m-1) in the same order as the lines of the text file - the Arc objects simply two member ints - i, and j (want to remove - not used to build MIP anymore but is used for flow constraints in benders - either use reverse adjacency list or use Dijkstra's as subroutine)
    * A vector of vectors of ints, outer vector indexed 0-(n-1), of adjacency list for each node. Each inner vector is a list of adjacent nodes;

### AdaptiveInstance

* Instance files:
    * AdaptiveInstance.name is the name of the instance, i.e. a base_name + graph_info(i.e. size, density). You can always construct appropriate filenames as /dat/base_name/AdaptiveInstance.name""
* Costs 
    * see AdaptiveInstance::initCosts() - pass interdiction = -1 to read from existing costs, or a positive interdiction value to generate arc costs and set interdiction costs to interdiction

## Experiment Scripts

* Experimental scripts are written in python files in scripts, found in [scripts](https://github.com/lucawrabetz/ARSPI/tree/master/scripts). These python files will call executeables compiled from C++. 
