#include "solvers.h"

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cout << argc << std::endl;
    std::cout << "Please input set name, min policies, and max policies." << std::endl;
    return 1;
  }
  const std::string set_name = argv[1];
  const int min_k = std::stoi(argv[2]);
  const int max_k = std::stoi(argv[3]);
  const std::vector<ASPI_Solver> solvers{
    MIP,
    BENDERS,
    ENUMERATION,
    GREEDY
  };
  RunAllInstancesInSetDirectory(min_k, max_k, set_name, solvers);
}
