# Adaptive Network Interdiction

This repository contains the source code to solve the adaptive network interdiction problem (adaptive shortest-path interdiction). The code to generate test instances and handle experiments can be found in [scripts](https://github.com/lucawrabetz/ARSPI/tree/master/scripts).

## Directory / file naming convention
* All instances in dat - grouped into sets of instances 
    * A set has a setname <base_name>-<date_created(mm_dd_yy)>-<repeat_number>, an instance/graph is called <setname>-<n>_<10*pr>
    * Create sets using graph gen - passing 3 args - base_name, max nodes, max pr (density param of erdos renyi graph)
        * Example call: python scripts/graphgen.py graphs 1000 0.6
        * Graph files (without costs) will be in dat/graphs-<date>-0
        * A single graph in the set, e.g. with n = 50 and pr = 0.5 will be at dat/graphs-<date>-0/graphs-<date>-0-50_5.txt
            * Notice we use 5 for the pr value instead of 0.5 to not use the period (.) in the naming
            * If the set was created and dat/graphs-<date>-0 already existed, it would create it under the name dat/graphs-<date>-1
        * Graph files are simply edge list files
        * Additionally, the set directory contains a creation logfile (dat/<setname>/<setname>.log), with information every graph in order on every line:
            * n
            * m
            * pr
            * density
            * filename (without the .txt extension)
        * Cost files are generated from within the cpp script, and python will pass a parameter to generate them or not depending on whether they already exist
        * Note: cost files depend on p, the number of scenarios, so a graph may have multiple different cost files, with the naming convention:
            * <instance_name>-costs_<p>.csv

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
