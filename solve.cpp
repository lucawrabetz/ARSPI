#include <fstream>
#include <unordered_map>

#include "solvers.h"

std::unordered_map<char, ASPI_Solver> flag_to_solver{
    {'m', MIP},
    {'b', BENDERS},
    {'e', ENUMERATION},
    {'g', GREEDY},
};

std::ofstream InitializeEmptyResultsFile(const std::string& set_name) {
  auto t = std::time(nullptr);
  auto tt = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tt, "%d-%m-%Y--%H-%M-%S");
  std::string today = oss.str();
  std::string base_name = DATA_DIRECTORY;
  base_name.append(set_name);
  base_name.append("/");
  base_name.append(set_name);
  base_name.append("_run_");
  base_name.append(today);
  std::string result_file_path = base_name;
  result_file_path.append(".csv");
  std::ofstream result_file(result_file_path, std::ios::app);
  if (!result_file.is_open()) {
    throw std::runtime_error(
        std::string("Failed to open file: ").append(result_file_path));
  }
  std::string csv_header = INSTANCE_INFO_COLUMN_HEADERS;
  csv_header.append(",").append(OUTPUT_COLUMN_HEADERS);
  result_file << csv_header << std::endl;
  return result_file;
}

std::ofstream InitializeAppendResultsFile(const std::string& append_file) {
  std::string result_file_path =
      std::string(DATA_DIRECTORY).append(append_file).append(".csv");
  std::ofstream result_file(result_file_path, std::ios::app);
  if (!result_file.is_open()) {
    throw std::runtime_error(
        std::string("Failed to open file: ").append(result_file_path));
  }
  return result_file;
}

void usage(char* name) {
  std::cout << "usage: ";
  std::cout << name << " setname [-s] [-m] [-g] [-a] [-t] [-u]" << std::endl;
  std::cout
      << "    [-s]: solver to use -- pass m for mip, b for benders, e for "
         "enumeration, g for greedy (case insensitive) -- default is mip "
      << std::endl
      << "    -- accepts multiple solvers / concatenation, i.e. mb or bm for "
         "mip and benders, ebgm for all solvers, etc."
      << std::endl;
  std::cout
      << "    [-m]: manual symmetric constraints (pass 0 for none, 1 for "
         "assignment constraints, or 2 for non-decreasing cluster constraints)"
      << std::endl;
  std::cout << "    [-g]: gurobi symmetry detection (pass -1 for automatic, 0 "
               "for off, 1 "
               "for conservative, or 2 for aggressive)"
            << std::endl;
  std::cout
      << "    [-a]: append mode, pass a file name in /dat without the "
         ".csv extension, which will be added. for example, if the "
         "script should append to test.csv, put test.csv in dat - "
         "i.e. ARSPI/dat/test.csv should exist, and pass ./bin/solve -a test"
      << std::endl;
  std::cout << "    [-t]: test mode, will run k = 1, ..., k_0 for every graph, "
               "budget = k_0 and fail if objectives of exact algorithms don't "
               "match on optimality (no check for run file)"
            << std::endl;
  std::cout << "    [-u]: uninterdicted mode, will run k = 0, for every "
               "instance using mip (no check for run file)"
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
  std::cout << "    please provide a file dat/<setname>_run.txt "
               "with four lines, each containing an integer representing "
               "the following parameters, in the following order:"
            << std::endl;
  std::cout << "    <min_policies>" << std::endl;
  std::cout << "    <max_policies>" << std::endl;
  std::cout << "    <min_budget>" << std::endl;
  std::cout << "    <max_budget>" << std::endl;
}

char lower_case(char c) {
  // Returns lower case value a-z if it is passed a lower case char or upper
  // case char a-z or A-Z, otherwise returns 'A'.
  if ('a' <= c && c <= 'z') return c;
  if ('A' <= c && c <= 'Z') return c + ('a' - 'A');
  return 'A';
}

int main(int argc, char* argv[]) {
  std::vector<ASPI_Solver> solvers;
  std::string append_file = "";  // Append file defaults to empty string, which
                                 // signals no append file to the library.
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
  bool test = false;
  bool uninterdicted = false;
  if (argc > 2) {
    if (argc > 11) {
      // Too many arguments.
      // Max args:
      // setname -m M -g G -s S -a A -u -t (11)
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
      if (strcmp(argv[i], "-s") == 0) {
        i++;
        std::string solver_string(argv[i]);
        if (solver_string.empty()) {
          usage(argv[0]);
          return 1;
        }
        for (char c : solver_string) {
          auto it = flag_to_solver.find(lower_case(c));
          if (it == flag_to_solver.end()) {
            usage(argv[0]);
            return 1;
          }
          solvers.push_back(it->second);
        }
      }
      if (strcmp(argv[i], "-a") == 0) {
        i++;
        append_file = argv[i];
      }
      if (strcmp(argv[i], "-u") == 0) {
        uninterdicted = true;
      }
      if (strcmp(argv[i], "-t") == 0) {
        test = true;
      }
    }
  }
  if (test) {
    TestOnAllInstancesInSetDirectory(set_name, manual_symmetry_constraints,
                                     gurobi_symmetry_detection,
                                     greedy_mip_gap_threshold);
    return 0;
  }
  std::ofstream result_file;
  if (append_file.empty()) {
    result_file = InitializeEmptyResultsFile(set_name);
  } else {
    result_file = InitializeAppendResultsFile(append_file);
  }
  if (uninterdicted) {
    UninterdictedRunOnAllInstancesInsetDirectory(set_name, result_file);
    return 0;
  }
  std::string run_file_path = "dat/" + set_name + "_run.txt";
  std::ifstream run_file(run_file_path);
  std::string param;
  int min_policies = -1;
  int max_policies = -1;
  int min_budget = -1;
  int max_budget = -1;
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
  if (solvers.empty()) solvers.push_back(MIP);  // Default solver is mip.
  for (ASPI_Solver solver : solvers) {
    SingleRunOnAllInstancesInSetDirectory(
        min_policies, max_policies, min_budget, max_budget, set_name,
        result_file, solver, manual_symmetry_constraints,
        gurobi_symmetry_detection, greedy_mip_gap_threshold);
  }
  result_file.close();
}
