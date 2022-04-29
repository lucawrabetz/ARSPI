#include "../inc/M3.h"

// use this file to test and debug
// write good scripts in other file

int main()
{
    int n=32;
    int p=13;
    int k=1;
    int M = 500;
    int interdiction_cost = 5;
    int min = 30;
    int max = 80;
    int mean = 50;
    int stddev = 10;
    int a = 50;
    int b = 150;
    float fraction = 0.1;

    const string setname = "tests-04_29_22-0";
    const string name = setname + "-" + to_string(n) + "_1_3";
    const string directory = "dat/" + setname + "/";
    const string filename = directory + name + ".txt";

    vector<int> fake_vec = vector<int>();
    vector<vector<int>> fake_vec2 = vector<vector<int>>();

    const Graph G = Graph(filename, n);
    int r_0=3;
    AdaptiveInstance m3 = AdaptiveInstance(p, k, r_0, G, directory, name); 
    // m3.initCosts();
    m3.initCosts(fraction, a, b, 3, G, true);

    // k = 1
    auto k1_solution = enumSolve(m3, G);
    printSolution(k1_solution, "k=1 solution exact");

    // k = 2
    auto k2_solution = extendByOne(k1_solution, m3, G);
    printSolution(k2_solution, "k=2 solution (heuristic, using exact k=1 as input)");

    m3.set_policies(2);
    auto k2exact_solution = enumSolve(m3, G);
    printSolution(k2exact_solution, "k=2 solution (exact)");

    // k = 3
    auto k3_solution = extendByOne(k2exact_solution, m3, G);
    printSolution(k3_solution, "k=3 solution (heuristic, using exact k=2 as input)");

    m3.set_policies(3);
    auto k3exact_solution = enumSolve(m3, G);
    printSolution(k3exact_solution, "k=3 solution (exact)");

    // k = 4
    // auto k4_solution = extendByOne(k3exact_solution, m3, G);
    // printSolution(k4_solution, "k=4 solution (heuristic, using exact k=3 as input)");

    // m3.set_policies(4);
    // auto k4exact_solution = enumSolve(m3, G);
    // printSolution(k4exact_solution, "k=4 solution (exact)");

    // // k = 5
    // auto k5_solution = extendByOne(k4exact_solution, m3, G);
    // printSolution(k5_solution, "k=5 solution (heuristic, using exact k=4 as input)");

    // m3.set_policies(5);
    // auto k5exact_solution = enumSolve(m3, G);
    // printSolution(k5exact_solution, "k=5 solution (exact)");










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

    // const Graph G = Graph(filename, n);
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
