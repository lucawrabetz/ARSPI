#include "../inc/M3.h"

// solve instance with enumeration algorithm
// example call:
// ./bin/enum set1_09-17-21 set1_09-17-21_26_0.4 26 5 1 6 500 30 80 10 1
// output: <objective> <runtime(ms)> (need to find objective)

int main(int argc, char* argv[]) {
    const string base_name = argv[1];
    const string name = argv[2];

    int n = stoi(argv[3]);
    int p = stoi(argv[4]);
    int k = stoi(argv[5]);
    int budget = stoi(argv[6]);
    int M = stoi(argv[7]);
    int min = stoi(argv[8]);
    int max = stoi(argv[9]);
    int interdiction_cost = stoi(argv[10]);
    int costs = stoi(argv[11]); // 1 if costs must be generated

    const string directory = "dat/" + base_name + "/";
    const string filename = directory + name + ".txt";

    const LayerGraph G = LayerGraph(filename, n);
    AdaptiveInstance m3 = AdaptiveInstance(p, k, budget, G, directory, name);
    if (costs == 1) {
        m3.initCosts(interdiction_cost, min, max);
    }
    else {m3.initCosts();}

    long begin = getCurrentTime();
    auto solution = enumSolve(m3, G);
    long runtime = getCurrentTime() - begin;

    cout << runtime << endl;
}
