#include "solvers.h"

int main(int argc, char *argv[]) {
  if (argc != 6) {
    std::cout << argc << std::endl;
    std::cout << "Please input set name, min policies, max policies, min "
                 "budget and max budget."
              << std::endl;
    return 1;
  }
  const std::string set_name = argv[1];
  const int min_k = std::stoi(argv[2]);
  const int max_k = std::stoi(argv[3]);
  const int min_budget = std::stoi(argv[4]);
  const int max_budget = std::stoi(argv[5]);
  const std::vector<ASPI_Solver> solvers{ENUMERATION};
  RunAllInstancesInSetDirectory(min_k, max_k, min_budget, max_budget, set_name,
                                solvers);
}
