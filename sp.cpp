#include "solvers.h"

// Solve an adaptive interdiction instance with the set partitioning model.
// Example call:
// ./bin/sp <set_name> <k_0> <n> <p> <k> <budget> <M>
// ./bin/sp aspi_testbed 3 27 5 1 3 500
// Output: <objective> <runtime(ms)> (<MIPGap> removed for now)
// The value k_0 is ONLY PASSED to form the filename for the graph and costs.
// It will not be used in any of the algorithms in the code.

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

    const Graph G = Graph(filename, n);
    AdaptiveInstance m3 = AdaptiveInstance(p, k, budget, G, directory, name);
    m3.ReadCosts();

    SetPartitioningModel sp_model = SetPartitioningModel(M, m3);
    sp_model.configureModel(G, m3);

    sp_model.solve();
    cout << "gurobi objective: " << sp_model.m3_model->get(GRB_DoubleAttr_ObjVal) << endl; 

    AdaptiveSolution sp_solution = sp_model.current_solution;
    sp_solution.ComputeAllObjectives(G, &m3);
    sp_solution.LogSolution(G, m3, "solution", true);
    long runtime = sp_solution.most_recent_solution_time();
    cout << sp_solution.worst_case_objective() << " " << runtime << endl;
}
