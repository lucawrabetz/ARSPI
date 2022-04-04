#include "../inc/M3.h"

// use this file to test and debug
// write good scripts in other file

int main()
{
    int n=100;
    int p=3;
    int k=1;
    int M = 500;
    int interdiction_cost = 5;
    int min = 30;
    int max = 80;
    int mean = 50;
    int stddev = 10;
    int a = 50;
    int b = 200;
    const string setname = "newdist-04_04_22-1";
    const string name = setname + "-" + to_string(n) + "_3";
    const string directory = "dat/" + setname + "/";
    const string filename = directory + name + ".txt";

    const LayerGraph G = LayerGraph(filename, n);
    int r_0=5;
    AdaptiveInstance m3 = AdaptiveInstance(p, k, r_0, G, directory, name); 
    m3.initCosts(interdiction_cost, a, b, 2);

    auto k_solution = enumSolve(m3, G);
    printSolution(k_solution, "k solution");

    auto kprime_solution = extendByOne(k_solution, m3, G);
    printSolution(kprime_solution, "k+1 solution (heuristic)");

    m3.set_policies(k+1);
    auto kk_solution = enumSolve(m3, G);
    printSolution(kk_solution, "k+1 solution (exact)");


    // int n=20;
    // int p=3;
    // int k=1;
    // int M = 500;
    // int interdiction_cost = 10;
    // int min = 30;
    // int max = 80;
    // int mean = 50;
    // int stddev = 10;
    // int a = 50;
    // int b = 100;
    // const string setname = "test-04_04_22-0";
    // const string name = setname + "-" + to_string(n) + "_3";
    // const string directory = "dat/" + setname + "/";
    // const string filename = directory + name + ".txt";

    // const LayerGraph G = LayerGraph(filename, n);
    // int r_0=n*0.3;
    // AdaptiveInstance m3 = AdaptiveInstance(p, k, r_0, G, directory, name); 
    // m3.initCosts();

    // auto k_solution = enumSolve(m3, G);
    // printSolution(k_solution, "k solution");

    // auto kprime_solution = extendByOne(k_solution, m3, G);
    // printSolution(kprime_solution, "k+1 solution (heuristic)");

    // m3.set_policies(k+1);
    // auto kk_solution = enumSolve(m3, G);
    // printSolution(kk_solution, "k+1 solution (exact)");
}
