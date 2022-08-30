#include "../inc/M3.h"

// solve instance with set partitioning model
// example call:
// ./bin/sp tests-06_20_22-0 2 3 32 3 1 5 500 50 200 0.1 0 3 

// output: <objective> <runtime(ms)> (<MIPGap> removed for now)

int main(int argc, char* argv[]) {
    const string set_name = argv[1];
    int p_0 = stoi(argv[2]);
    int k_0 = stoi(argv[3]);
    int n = stoi(argv[4]);
    int p = stoi(argv[5]);
    int k = stoi(argv[6]);
    int budget = stoi(argv[7]);
    int M = stoi(argv[8]);
    int min = stoi(argv[9]);
    int max = stoi(argv[10]);
    int fraction = stoi(argv[11]);
    int costs = stoi(argv[12]); // 1 if costs must be generated
    int dist = stoi(argv[13]);

    const string name = set_name + "-" + to_string(n) + "_" + to_string(p_0) + "_" + to_string(k_0);
    const string directory = "dat/" + set_name + "/";
    const string filename = directory + name + ".txt";

    const Graph G = Graph(filename, n);
    AdaptiveInstance m3 = AdaptiveInstance(p, k, budget, G, directory, name);
    if (costs == 1) {
        m3.initCosts(fraction, min, max, dist, G, true);
    }
    else {m3.initCosts(fraction, min, max, dist, G, false);}

    SetPartitioningModel sp_model = SetPartitioningModel(M, m3);
    sp_model.configureModel(G, m3);

    sp_model.solve();

    AdaptiveSolution sp_solution = sp_model.current_solution;
    sp_solution.computeAllObjectives(G, m3);
    long runtime = sp_solution.most_recent_solution_time;

    cout << sp_solution.worst_case_objective << " " << runtime << endl;
}
