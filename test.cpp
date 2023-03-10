#include "solvers.h"

// Simple simple test coverage - check that the 3 solvers (set partitioning, enumeration, benders)
// produce the same result (objective value) on a few small instances.
// Every test prints the algorithm and instance name1 followed by "PASS" if it passes, "FAIL" if it fails.
// CURRENT: just running sp, need to add enum and benders.

int main(int argc, char*argv[]) {
    const string set_name = "tests";
    const string directory = "dat/" + set_name + "/";
    // First test - 5 node instance with 1 follower:
    int k_0 = 1;
    int n = 5;
    int p = 1;
    int k = 1;
    int budget = 1;
    int M = 50;
    const string name1 = set_name + "-" + to_string(n) + "_" + to_string(k_0);
    const string filename1 = directory + name1 + ".txt";
    const Graph G1 = Graph(filename1, n);
    AdaptiveInstance test1(p, k, budget, G1, directory, name1);
    test1.ReadCosts();
    GRBEnv* env = new GRBEnv(); // Initialize global gurobi environment.
    env->set(GRB_DoubleParam_TimeLimit, 3600); // Set time limit to 1 hour.
    // SetPartitioningBenders benders = SetPartitioningBenders(M, test1);
    // benders.ConfigureSolver(G1, test1);
    // benders.Solve();

    SetPartitioningModel sp = SetPartitioningModel(M, test1, env);
    sp.ConfigureSolver(G1, test1);
    sp.Solve();
    AdaptiveSolution sp_solution = sp.current_solution();
    sp_solution.ComputeAllObjectives(G1, &test1);
    cout << "Instance 1, k = 1: " << sp_solution.worst_case_objective() << endl;

    // Second test - 6 node instance with 3 followers, k = 1,...,3:
    // k_0 = 3;
    // n = 6;
    // p = 3;
    // const string name2 = set_name + "-" + to_string(n) + "_" + to_string(k_0);
    // const string filename2 = directory + name2 + ".txt";
    // cout << filename2 << endl;
    // const Graph G2 = Graph(filename2, n);

    // for (int k=1; k<4; k++) {
    //     AdaptiveInstance test2(p, k, budget, G2, directory, name2);
    //     test2.ReadCosts();
    //     SetPartitioningModel sp = SetPartitioningModel(M, test2);
    //     sp.ConfigureSolver(G2, test2);
    //     sp.Solve();
    //     AdaptiveSolution sp_solution = sp.current_solution();
    //     sp_solution.ComputeAllObjectives(G2, &test2);
    //     cout << "Instance 2, k = " << k << ": " << sp_solution.worst_case_objective() << endl;
    // }
}
