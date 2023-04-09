#include <iomanip>
#include <limits>

#include "solvers.h"

// Simple simple test coverage - check that the 3 solvers (set partitioning,
// enumeration, benders) produce the same result (objective value) on a few
// small instances. Every test prints the algorithm and instance name1 followed
// by "PASS" if it passes, "FAIL" if it fails.

typedef std::numeric_limits<double> dbl;

void SolveAndPrintTest(const Graph& G, AdaptiveInstance& instance,
                       std::vector<ASPI_Solver>& solvers, GRBEnv* env,
                       int debug = 0) {
  std::vector<double> adaptive_objectives;
  double greedy_adaptive_objective;
  for (const auto& solver : solvers) {
    if (solver == MIP) {
      SetPartitioningModel sp = SetPartitioningModel(instance, env);
      sp.ConfigureSolver(G, instance);
      AdaptiveSolution sp_solution = sp.Solve(G, instance, env);
      if (debug == 2)
        sp_solution.LogSolution(G, "MIP - " + instance.name(), true);
      else if (debug == 1)
        sp_solution.LogSolution(G, "MIP - " + instance.name(), false);
      adaptive_objectives.push_back(sp_solution.worst_case_objective());
    } else if (solver == BENDERS) {
      SetPartitioningBenders benders =
          SetPartitioningBenders(instance, env);
      benders.ConfigureSolver(G, instance);
      AdaptiveSolution benders_solution = benders.Solve(G, instance, env);
      if (debug == 2)
        benders_solution.LogSolution(G, "Benders - " + instance.name(), true);
      else if (debug == 1)
        benders_solution.LogSolution(G, "Benders - " + instance.name(), false);
      adaptive_objectives.push_back(benders_solution.worst_case_objective());
    } else if (solver == ENUMERATION) {
      AdaptiveSolution enum_solution = EnumSolve(instance, G, env);
      if (debug == 2)
        enum_solution.LogSolution(G, "Enumeration - " + instance.name(), true);
      else if (debug == 1)
        enum_solution.LogSolution(G, "Enumeration - " + instance.name(), false);
      adaptive_objectives.push_back(enum_solution.worst_case_objective());
    } else if (solver == GREEDY) {
      AdaptiveSolution greedy_solution = GreedyAlgorithm(instance, G, env);
      // enum_solution.ComputeAllObjectives(G, instance);
      if (debug == 2)
        greedy_solution.LogSolution(G, "Enumeration - " + instance.name(),
                                    true);
      else if (debug == 1)
        greedy_solution.LogSolution(G, "Enumeration - " + instance.name(),
                                    false);
      // We won't add the objective to adaptive objectives, as we don't want to
      // use it in the test, which is only for exact algorithms. However, we
      // save it to the stand alone value, greedy_adaptive_objective, so we can
      // print it next to the test outcome.
      greedy_adaptive_objective = greedy_solution.worst_case_objective();
      // adaptive_objectives.push_back(enum_solution.worst_case_objective());
    }
  }
  double prev = adaptive_objectives[0];
  for (double& obj : adaptive_objectives) {
    if (std::abs(obj - prev) >= EPSILON) {
      std::cout << "FAIL: " << instance.name() << std::endl;
      std::cout.precision(dbl::max_digits10 - 1);
      for (double& x : adaptive_objectives) {
        std::cout << std::scientific << x << std::endl;
      }
      return;
    }
    prev = obj;
  }
  std::cout << "PASS: " << instance.name() << ", k = " << instance.policies()
            << ", exact objective: " << adaptive_objectives[0]
            << ", greedy approximation: " << greedy_adaptive_objective
            << std::endl;
}

// void RunKTestsForInstance(const Graph& G, AdaptiveInstance& instance, const std::vector<ASPI_Solver>& solvers, GRBEnv* env, int debug = 0){
//   k_0 = 3;
//   n = 6;
//   p = 3;
//   const std::string name2 =
//       set_name + "-" + std::to_string(n) + "_" + std::to_string(k_0);
//   const std::string filename2 = directory + name2 + ".txt";
//   const Graph G2 = Graph(filename2, n);
//   for (int k = 1; k < 4; ++k) {
//     AdaptiveInstance test2(p, k, budget, G2, directory, name2);
//     test2.ReadCosts(interdiction_delta);
//     SolveAndPrintTest(G2, test2, solvers, env, M, debug_level);
//   }
// }

int main() {
  GRBEnv* env = new GRBEnv();  // Initialize global gurobi environment.
  env->set(GRB_DoubleParam_TimeLimit, 3600);  // Set time limit to 1 hour.
  const int debug_level = 0;
  const std::string set_name = "tests";
  const std::string directory = "dat/" + set_name + "/";
  // First test - 5 node instance with 1 follower:
  int k_0 = 1;
  int n = 5;
  int p = 1;
  int k = 1;
  int budget = 1;
  const std::string name1 =
      set_name + "-" + std::to_string(n) + "_" + std::to_string(k_0);
  const std::string filename1 = directory + name1 + ".txt";
  const Graph G1 = Graph(filename1, n);
  AdaptiveInstance test1(p, k, budget, G1, directory, name1);
  test1.ReadCosts();
  std::vector<ASPI_Solver> solvers{MIP, BENDERS, ENUMERATION, GREEDY};
  SolveAndPrintTest(G1, test1, solvers, env, debug_level);
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
    test2.ReadCosts();
    SolveAndPrintTest(G2, test2, solvers, env, debug_level);
  }
  const std::string synthetic_set = "aspi_testbed";
  const std::string synthetic_dir = "dat/" + synthetic_set + "/";
  // Synthetic test 1:
  k_0 = 3;
  n = 52;
  p = 5;
  budget = 3;
  const std::string name3 =
      synthetic_set + "-" + std::to_string(n) + "_" + std::to_string(k_0);
  const std::string filename3 = synthetic_dir + name3 + ".txt";
  const Graph G3 = Graph(filename3, n);
  for (int k = 1; k < 4; ++k) {
    AdaptiveInstance test3(p, k, budget, G3, synthetic_dir, name3);
    test3.ReadCosts();
    SolveAndPrintTest(G3, test3, solvers, env, debug_level);
  }
  k_0 = 5;
  const std::string name4 =
      synthetic_set + "-" + std::to_string(n) + "_" + std::to_string(k_0);
  const std::string filename4 = synthetic_dir + name4 + ".txt";
  const Graph G4 = Graph(filename4, n);
  for (int k = 1; k < 4; ++k) {
    AdaptiveInstance test4(p, k, budget, G4, synthetic_dir, name4);
    test4.ReadCosts();
    SolveAndPrintTest(G4, test4, solvers, env, debug_level);
  }
}
