#include "../inc/M3.h"

// solve instance with set partitioning model
// example call:
// ./bin/sp graphs 3 27 5 1 2 500

// output: <objective> <runtime(ms)> (<MIPGap> removed for now)

int main(int argc, char* argv[]) {
    const string set_name = argv[1];
    int k_0 = stoi(argv[2]);
    int n = stoi(argv[3]);
    int p = stoi(argv[4]);
    int k = stoi(argv[5]);
    int budget = stoi(argv[6]);
    int M = stoi(argv[7]);

    const string name = set_name + "-" + to_string(n) + "_" + to_string(k_0);
    const string directory = "dat/" + set_name + "/";
    const string filename = directory + name + ".txt";

    const Graph G = Graph(filename, n, k_0);
    cout << "hello i made it" << endl;
    AdaptiveInstance m3 = AdaptiveInstance(p, k, budget, G, directory, name);
    m3.initCosts(1, 1, 1, 1, G, false);

    SetPartitioningModel sp_model = SetPartitioningModel(M, m3);
    sp_model.configureModel(G, m3);

    sp_model.solve();

    AdaptiveSolution sp_solution = sp_model.current_solution;
    sp_solution.computeAllObjectives(G, m3);
    sp_solution.logSolution(G, m3, "solution", true);
    long runtime = sp_solution.most_recent_solution_time;

    cout << sp_solution.worst_case_objective << " " << runtime << endl;
}
