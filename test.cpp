#include "solvers.h"

// Simple simple test coverage - check that the 3 solvers (set partitioning,
// enumeration, benders) produce the same result (objective value) on a few
// small instances. Every test prints the algorithm and instance name1 followed
// by "PASS" if it passes, "FAIL" if it fails.


void SolveAndPrint(const Graph& G, AdaptiveInstance& instance,
                   std::vector<ASPI_Solver>& solvers, GRBEnv* env, int big_m) {
  for (const auto& solver : solvers) {
    if (solver == MIP) {
      SetPartitioningModel sp = SetPartitioningModel(big_m, instance, env);
      sp.ConfigureSolver(G, instance);
      AdaptiveSolution sp_solution = sp.Solve();
      // sp_solution.ComputeAllObjectives(G, &instance);
      sp_solution.LogSolution(G, "MIP - " + instance.name(), false);
    } else if (solver == BENDERS) {
      SetPartitioningBenders benders =
          SetPartitioningBenders(big_m, instance, env);
      benders.ConfigureSolver(G, instance);
      AdaptiveSolution benders_solution = benders.Solve();
      benders_solution.LogSolution(G, "Benders - " + instance.name(), false);
    }
  }
}

int main() {
  GRBEnv* env = new GRBEnv();  // Initialize global gurobi environment.
  env->set(GRB_DoubleParam_TimeLimit, 3600);  // Set time limit to 1 hour.
  const std::string set_name = "tests";
  const std::string directory = "dat/" + set_name + "/";
  // First test - 5 node instance with 1 follower:
  int k_0 = 1;
  int n = 5;
  int p = 1;
  int k = 1;
  int budget = 1;
  int M = 100;
  int interdiction_delta = 10;
  const std::string name1 =
      set_name + "-" + std::to_string(n) + "_" + std::to_string(k_0);
  const std::string filename1 = directory + name1 + ".txt";
  const Graph G1 = Graph(filename1, n);
  AdaptiveInstance test1(p, k, budget, G1, directory, name1);
  test1.ReadCosts(interdiction_delta);
  std::vector<ASPI_Solver> solvers{MIP, BENDERS};
  SolveAndPrint(G1, test1, solvers, env, M);
  // Second test - 6 node instance with 3 followers, k = 1,...,3:
  k_0 = 3;
  n = 6;
  p = 3;
  const std::string name2 =
      set_name + "-" + std::to_string(n) + "_" + std::to_string(k_0);
  const std::string filename2 = directory + name2 + ".txt";
  const Graph G2 = Graph(filename2, n);
  for (int k = 1; k < 4; ++k) {
    AdaptiveInstance test2(p, k, budget, G2, directory, name2);
    test2.ReadCosts(interdiction_delta);
    SolveAndPrint(G2, test2, solvers, env, M);
  }
  const std::string synthetic_set = "aspi_testbed";
  const std::string synthetic_dir = "dat/" + synthetic_set + "/";
  // Synthetic test 1:
  k_0 = 3;
  n = 52;
  p = 5;
  budget = 3;
  M = 1000;
  interdiction_delta = 100;
  const std::string name3 =
      synthetic_set + "-" + std::to_string(n) + "_" + std::to_string(k_0);
  const std::string filename3 = synthetic_dir + name3 + ".txt";
  const Graph G3 = Graph(filename3, n);
  for (int k = 1; k < 4; ++k) {
    AdaptiveInstance test3(p, k, budget, G3, synthetic_dir, name3);
    test3.ReadCosts(interdiction_delta);
    SolveAndPrint(G3, test3, solvers, env, M);
  }
  k_0 = 5;
  const std::string name4 =
      synthetic_set + "-" + std::to_string(n) + "_" + std::to_string(k_0);
  const std::string filename4 = synthetic_dir + name4 + ".txt";
  const Graph G4 = Graph(filename4, n);
  for (int k = 1; k < 4; ++k) {
    AdaptiveInstance test4(p, k, budget, G4, synthetic_dir, name4);
    test4.ReadCosts(interdiction_delta);
    SolveAndPrint(G4, test4, solvers, env, M);
  }
}
