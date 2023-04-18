#include <cstring>
#include <dirent.h>
#include <cmath>

#include <iomanip>
#include <limits>
#include <unordered_map>

#include "solvers.h"

typedef std::numeric_limits<double> dbl;

bool IsCostFile(const std::string& name) {
  size_t pos = name.find(".csv");
  if (pos == std::string::npos) return false;
  return true;
}

bool IsRunFile(const std::string& name) {
  size_t pos = name.find("run");
  if (pos == std::string::npos) return false;
  return true;
}

int DigitsToInt(const std::vector<int> num_digits) {
  int i = num_digits.size() - 1;
  int num = 0;
  for (int x : num_digits) {
    num += x * std::pow(10, i);
    i--;
  }
  return num;
}

std::pair<int, int> GetNodesKZero(const std::string& name, int start) {
  int i = start;
  std::vector<int> n_digits;
  std::vector<int> kzero_digits;
  while ('0' <= name[i]  && name[i] <= '9') {
    n_digits.push_back(name[i] - '0');
    i++;
  }
  i++;
  while ('0' <= name[i]  && name[i] <= '9') {
    kzero_digits.push_back(name[i] - '0');
    i++;
  }
  return {DigitsToInt(n_digits), DigitsToInt(kzero_digits)};
}

std::pair<int, int> GetScenariosID(const std::string& name) {
  int i = name.size()-5;
  std::vector<int> p_digits;
  std::vector<int> id_digits;
  while ('0' <= name[i]  && name[i] <= '9') {
    id_digits.insert(id_digits.begin(), name[i] - '0');
    i--;
  }
  i--;
  while ('0' <= name[i]  && name[i] <= '9') {
    p_digits.insert(p_digits.begin(), name[i] - '0');
    i--;
  }
  return {DigitsToInt(p_digits), DigitsToInt(id_digits)};
  return {};
}

const std::string CSV_HEADER = "set_name,instance_name,nodes,arcs,k_zero,density,scenarios,budget,policies,MIP_OPTIMAL,MIP_objective,MIP_gap,MIP_time,BENDERS_OPTIMAL,BENDERS_objective,BENDERS_gap,BENDERS_time,BENDERS_cuts_rounds,ENUMERATION_OPTIMAL,ENUMERATION_objective,ENUMERATION_time,GREEDY_objective,GREEDY_time";

std::string SolveAndPrintTest(const std::string& set_name, const ProblemInput& problem,
                       std::vector<ASPI_Solver>& solvers, int debug = 0) {
  std::vector<double> adaptive_objectives;
  double greedy_adaptive_objective;
  std::string final_csv_string = set_name + "," + problem.instance_.name();
  final_csv_string += "," + std::to_string(problem.G_.nodes());
  final_csv_string += "," + std::to_string(problem.G_.arcs());
  final_csv_string += "," + std::to_string(problem.k_zero_);
  double nodes = problem.G_.nodes();
  double arcs = problem.G_.arcs();
  double density = (arcs) / ((nodes) * (nodes-1));
  final_csv_string += "," + std::to_string(density);
  final_csv_string += "," + std::to_string(problem.instance_.scenarios());
  final_csv_string += "," + std::to_string(problem.budget_);
  final_csv_string += "," + std::to_string(problem.policies_);
  std::string log_line = " ---------- ";
  for (const auto& solver : solvers) {
    if (solver == MIP) {
      SetPartitioningModel sp = SetPartitioningModel(problem);
      sp.ConfigureSolver(problem);
      AdaptiveSolution sp_solution = sp.Solve(problem);
      if (debug == 2)
        sp_solution.LogSolution(problem, true);
      else if (debug == 1)
        sp_solution.LogSolution(problem, false);
      adaptive_objectives.push_back(sp_solution.worst_case_objective());
      std::string optimal;
      if (sp_solution.optimal()) {optimal = "OPTIMAL";}
      else {optimal = "NOT OPTIMAL";}
      final_csv_string += "," + optimal;
      final_csv_string += "," + std::to_string(sp_solution.worst_case_objective());
      final_csv_string += "," + std::to_string(sp_solution.mip_gap());
      final_csv_string += "," + std::to_string(sp_solution.solution_time());
      log_line += "MIP: ";
      log_line += optimal + " - ";
      log_line += std::to_string(sp_solution.worst_case_objective());
      log_line += ", " + std::to_string(sp_solution.solution_time());
      log_line += "ms ----- ";
    } else if (solver == BENDERS) {
      SetPartitioningBenders benders = SetPartitioningBenders(problem);
      benders.ConfigureSolver(problem);
      AdaptiveSolution benders_solution = benders.Solve(problem);
      if (debug == 2)
        benders_solution.LogSolution(problem, true);
      else if (debug == 1)
        benders_solution.LogSolution(problem, false);
      adaptive_objectives.push_back(benders_solution.worst_case_objective());
      std::string optimal;
      if (benders_solution.optimal()) {optimal = "OPTIMAL";}
      else {optimal = "NOT OPTIMAL";}
      final_csv_string += "," + optimal;
      final_csv_string += "," + std::to_string(benders_solution.worst_case_objective());
      final_csv_string += "," + std::to_string(benders_solution.mip_gap());
      final_csv_string += "," + std::to_string(benders_solution.solution_time());
      final_csv_string += "," + std::to_string(benders_solution.lazy_cuts_rounds());
      log_line += "BENDERS: ";
      log_line += optimal + " - ";
      log_line += std::to_string(benders_solution.worst_case_objective());
      log_line += ", " + std::to_string(benders_solution.solution_time());
      log_line += "ms ----- ";
    } else if (solver == ENUMERATION) {
      AdaptiveSolution enum_solution = EnumSolve(problem);
      if (debug == 2)
        enum_solution.LogSolution(problem, true);
      else if (debug == 1)
        enum_solution.LogSolution(problem, false);
      adaptive_objectives.push_back(enum_solution.worst_case_objective());
      std::string optimal;
      if (enum_solution.optimal()) {optimal = "OPTIMAL";}
      else {optimal = "NOT OPTIMAL";}
      final_csv_string += "," + optimal;
      final_csv_string += "," + std::to_string(enum_solution.worst_case_objective());
      final_csv_string += "," + std::to_string(enum_solution.solution_time());
      log_line += "ENUMERATION: ";
      log_line += optimal + " - ";
      log_line += std::to_string(enum_solution.worst_case_objective());
      log_line += ", " + std::to_string(enum_solution.solution_time());
      log_line += "ms ----- ";
    } else if (solver == GREEDY) {
      AdaptiveSolution greedy_solution = GreedyAlgorithm(problem);
      if (debug == 2)
        greedy_solution.LogSolution(problem, true);
      else if (debug == 1)
        greedy_solution.LogSolution(problem, false);
      // We won't add the objective to adaptive objectives, as we don't want to
      // use it in the test, which is only for exact algorithms. However, we
      // save it to the stand alone value, greedy_adaptive_objective, so we can
      // print it next to the test outcome.
      greedy_adaptive_objective = greedy_solution.worst_case_objective();
      final_csv_string += "," + std::to_string(greedy_solution.worst_case_objective());
      final_csv_string += "," + std::to_string(greedy_solution.solution_time());
      log_line += "GREEDY: ";
      log_line += std::to_string(greedy_solution.worst_case_objective());
      log_line += ", " + std::to_string(greedy_solution.solution_time());
      log_line += "ms ---------- ";
      // adaptive_objectives.push_back(enum_solution.worst_case_objective());
    }
  }
  std::cout << log_line << std::endl;
  double prev = adaptive_objectives[0];
  // FOR TESTING
  // for (double& obj : adaptive_objectives) {
  //   if (std::abs(obj - prev) >= EPSILON) {
  //     std::cout << "FAIL: " << problem.instance_.name() << std::endl;
  //     std::cout.precision(dbl::max_digits10 - 1);
  //     for (double& x : adaptive_objectives) {
  //       std::cout << std::scientific << x << std::endl;
  //     }
  //   }
  //   prev = obj;
  // }
  // std::cout << "------------------------  PASS, EXACT OBJECTIVE: " << adaptive_objectives[0] << ", GREEDY APPROXIMATION: " << greedy_adaptive_objective
  //           << "  ------------------------" << std::endl;
  return final_csv_string;
}

void RunAllInstancesInSetDirectory(const int max_policies, const std::string& set_name) {
  GRBEnv* env = new GRBEnv();  // Initialize global gurobi environment.
  // Use seconds, since gurobi takes the parameter value in seconds.
  env->set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_S);  // Set time limit.
  std::cout << std::endl;

  /* Max Policies */
  /* ------------ */
  std::vector<ASPI_Solver> solvers{MIP, BENDERS, ENUMERATION, GREEDY};

  /* Set Name */
  std::string full_path_s = DATA_DIRECTORY + set_name; 
  char* full_path = new char[full_path_s.length() + 1];
  std::strcpy(full_path, full_path_s.c_str());
  DIR* set_directory = opendir(full_path);
  /* -------- */

  struct dirent* entity;
  entity = readdir(set_directory);
  auto t = std::time(nullptr);
  auto tt = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tt, "%d-%m-%Y--%H-%M-%S");
  std::string today = oss.str();
  std::string base_name = DATA_DIRECTORY + set_name + "/" + set_name + "_run_" + today;
  std::string result_file_name =  base_name + ".csv";
  std::ofstream result_file(result_file_name);
  result_file << CSV_HEADER << std::endl;
  while (entity != NULL) {
    std::string name = entity->d_name;
    entity = readdir(set_directory);
    if (name.size() < set_name.size())
      continue;  // This is not a data file, data files start with <set_name> so
                 // are at least as long as set_name.size().
    std::string set_name_in_name(name.begin(), name.begin() + set_name.size());
    if (set_name_in_name != set_name) continue; // For sanity, check that the data file matches the set_name which should always be the case.
    if (IsRunFile(name)) continue; // Cannot run on a output data file (in case we are running again on this directory).
                                   // Since we include full timestamps in the result file names we may do this, e.g. with a different max policies, and won't ovewrite result file names.
    // We will loop through all cost files (which are analogous to all InstanceInputs) parse their names to get params for both GraphInput and InstanceInput, and then run the our single run function for every k from 1 to max k (as long as k <= scenarios).
    if (IsCostFile(name)) {
      std::pair<int, int> n_kzero = GetNodesKZero(name, set_name.size()+1);
      int nodes = n_kzero.first;
      int kzero = n_kzero.second;
      std::pair<int, int> scenarios_id = GetScenariosID(name);
      int scenarios = scenarios_id.first;
      int id = scenarios_id.second;
      GraphInput graph_input(set_name, nodes, kzero);
      InstanceInput instance_input(graph_input, scenarios, id);
      for (int k = 1; k <= max_policies; k++) {
        if (k > instance_input.scenarios_) break;
        const ProblemInput problem(instance_input, k, env);
        std::cout << "RUNNING INSTANCE: " << problem.instance_.name() << ", K = " << std::to_string(k) << std::endl;
        std::string result = SolveAndPrintTest(set_name, problem, solvers);
        std::cout << std::endl;
        result_file << result << std::endl;
      }
    } 
  }
  result_file.close();
  closedir(set_directory);
}

int main(int argc, char *argv[]) {
  const std::string set_name = argv[1];
  const int max_k = std::stoi(argv[2]);
  long before = GetCurrentTime();
  RunAllInstancesInSetDirectory(max_k, set_name);
  long after = GetCurrentTime();
  long dur = after - before;
  std::cout << dur << std::endl;
}
