#include "solvers.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << argc << std::endl;
    std::cout << "Please input set name." << std::endl;
    return 1;
  }
  const std::string set_name = argv[1];
  UninterdictedObjectiveForAllInstances(set_name);
}
