#include "solvers.h"

void usage(char* name) {
  std::cout << "usage: ";
  std::cout << name << " setname [-m] [-g]" << std::endl;
  std::cout
      << "    [-m]: manual symmetric constraints (pass 0 for none, 1 for "
         "assignment constraints, or 2 for non-decreasing cluster constraints)"
      << std::endl;
  std::cout << "    [-g]: gurobi symmetry detection (pass -1 for automatic, 0 "
               "for off, 1 "
               "for conservative, or 2 for aggressive)"
            << std::endl;
}

void runfile_instructions(bool param_error = false) {
  std::cout << "run parameter file is not set up correctly:" << std::endl;
  if (param_error) {
    std::cout
        << "    check your parameters: (min_policies -- max_policies) and "
           "(min_budget -- max_budget) should represent "
           "non-decreasing ranges"
        << std::endl;
    return;
  }
  std::cout
      << "    please provide a file dat/<setname>_run.txt in the root project "
         "directory with four lines, each containing an integer representing "
         "the following parameters, in the following order:"
      << std::endl;
  std::cout << "    <min_policies>" << std::endl;
  std::cout << "    <max_policies>" << std::endl;
  std::cout << "    <min_budget>" << std::endl;
  std::cout << "    <max_budget>" << std::endl;
}

int main(int argc, char* argv[]) {
  int manual_symmetry_constraints =
      MANUAL_SYMMETRY_NONE;  // (Constants defined in solvers.h) - 0 = none, 1 =
                             // assignment type constraints, 2 = non-decreasing
                             // cluster constraints. Default to 0.
  int gurobi_symmetry_detection =
      -1;  // -1 for automatic, 0 = off, 1 = conservative, 2 = aggressive.
           // Default to -1.
  double greedy_mip_gap_threshold = -1;  // Always -1 (causes it to be ignored
                                         // in the greedy algorithm) for now.
  if (argc < 2) {
    usage(argv[0]);
    return 1;
  }
  const std::string set_name = argv[1];
  int min_policies = -1;
  int max_policies = -1;
  int min_budget = -1;
  int max_budget = -1;
  std::string run_file_path = "dat/" + set_name + "_run.txt";
  std::ifstream run_file(run_file_path);
  std::string param;
  if (run_file.is_open()) {
    int i = 0;
    while (run_file) {
      std::getline(run_file, param);
      if (i == 0) {
        // Min policies.
        if (param.empty()) {
          runfile_instructions();
          return 1;
        }
        min_policies = std::stoi(param);
      } else if (i == 1) {
        // Max policies.
        if (param.empty()) {
          runfile_instructions();
          return 1;
        }
        max_policies = std::stoi(param);
      } else if (i == 2) {
        // Min budget.
        if (param.empty()) {
          runfile_instructions();
          return 1;
        }
        min_budget = std::stoi(param);
      } else if (i == 3) {
        // Max budget.
        if (param.empty()) {
          runfile_instructions();
          return 1;
        }
        max_budget = std::stoi(param);
      } else if (i == 4) {
        // End of file, do nothing.
        continue;
      } else if (i > 4) {
        if (param.empty()) {
          continue;
        } else {
          runfile_instructions();
          return 1;
        }
      }
      i++;
    }
  }
  if (std::min(std::min(min_policies, min_budget),
               std::min(max_policies, max_budget)) < 0) {
    runfile_instructions();
    return 1;
  }
  if (min_policies > max_policies) {
    runfile_instructions(true);
    return 1;
  }
  if (min_budget > max_budget) {
    runfile_instructions(true);
    return 1;
  }
  if (argc > 2) {
    if (argc == 3 || argc == 5) {
      std::cout << "argc: " << argc << std::endl;
      // If optional args -m or -g are passed, they must also have a value.
      // So argc should be 2 (if neither are passed) or 4 or 6.
      usage(argv[0]);
      return 1;
    }
    if (argc > 6) {
      // Too many arguments.
      usage(argv[0]);
      return 1;
    }
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-m") == 0) {
        i++;
        int arg = std::stoi(argv[i]);
        if (arg > 2 || arg < 0) {
          // Manual symmetry constraint parameter should be 0, 1, or 2.
          usage(argv[0]);
          return 1;
        }
        manual_symmetry_constraints = arg;
      }
      if (strcmp(argv[i], "-g") == 0) {
        i++;
        int arg = std::stoi(argv[i]);
        if (arg > 2 || arg < -1) {
          // Gurobi symmetry detection parameter should be -1, 0, 1, or 2.
          usage(argv[0]);
          return 1;
        }
        gurobi_symmetry_detection = arg;
      }
    }
  }
  const std::vector<ASPI_Solver> solvers{MIP};
  RunAllInstancesInSetDirectory(
      min_policies, max_policies, min_budget, max_budget, set_name, solvers,
      manual_symmetry_constraints, gurobi_symmetry_detection,
      greedy_mip_gap_threshold);
}
