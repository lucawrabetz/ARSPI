#include "solvers.h"

#include <string>

long GetCurrentTime() {
  // Helper function to get current time in milliseconds.
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

void PrintSeparator() {}

Graph::Graph(const std::string& filename, int nodes) {
  // Graph constructor - read arc list from file.
  std::string line, word;
  std::ifstream myfile(filename);
  if (myfile.is_open()) {
    nodes_ = nodes;
    int arc_index = 0;
    arc_index_hash_ = std::vector<std::vector<int>>(nodes_);
    adjacency_list_ = std::vector<std::vector<int>>(nodes_);
    rev_arc_index_hash_ = std::vector<std::vector<int>>(nodes_);
    rev_adjacency_list_ = std::vector<std::vector<int>>(nodes_);
    while (getline(myfile, line)) {
      int counter = 0;
      int i, j;
      std::stringstream str(line);
      while (getline(str, word, ' ')) {
        if (counter == 0) {
          i = std::stoi(word);
        } else if (counter == 1) {
          j = std::stoi(word);
        }
        ++counter;
      }
      adjacency_list_[i].push_back(j);
      arc_index_hash_[i].push_back(arc_index);
      rev_adjacency_list_[j].push_back(i);
      rev_arc_index_hash_[j].push_back(arc_index);
      ++arc_index;
    }
    arcs_ = arc_index;
  }
}

void Graph::PrintArc(int a, int i, int index, bool rev) const {
  if (rev) {
    int j = rev_adjacency_list_[i][index];
    std::cout << a << ": (" << j << ", " << i << ")";
  } else {
    int j = adjacency_list_[i][index];
    std::cout << a << ": (" << i << ", " << j << ")";
  }
}

void Graph::PrintGraph(bool rev) const {
  // Print arc summary of a graph.
  std::cout << "n: " << nodes_ << ", m: " << arcs_ << std::endl;
  for (int i = 0; i < nodes_; i++) {
    int index = 0;
    for (int a : arc_index_hash_[i]) {
      PrintArc(a, i, index);
      std::cout << std::endl;
      index++;
    }
  }
  // In case you want to test that adjacency_list_ and rev_adjacency_list_ match
  // for sanity.
  if (rev) {
    for (int i = 0; i < nodes_; i++) {
      int index = 0;
      for (int a : rev_arc_index_hash_[i]) {
        PrintArc(a, i, index, true);
        std::cout << std::endl;
        index++;
      }
    }
  }
}

void Graph::PrintGraphWithCosts(const std::vector<std::vector<int>>& costs,
                                int interdiction_delta) const {
  // Print arc summary of a graph with costs.
  std::cout << "n: " << nodes_ << ", m: " << arcs_ << std::endl;
  int scenarios = costs.size();
  for (int i = 0; i < nodes_; i++) {
    int index = 0;
    for (int a : arc_index_hash_[i]) {
      PrintArc(a, i, index);
      std::cout << ", costs: ";
      for (int q = 0; q < scenarios; q++) {
        std::cout << q + 1 << ": " << costs[q][a] << ", ";
      }
      std::cout << "interdiction delta: " << interdiction_delta << std::endl;
    }
  }
}

AdaptiveInstance::AdaptiveInstance(AdaptiveInstance& adaptive_instance,
                                   std::vector<int>& keep_scenarios) {
  // Copy constructor only keeping a subset of the scenarios
  nodes_ = adaptive_instance.nodes_;
  scenarios_ = keep_scenarios.size();
  interdiction_delta_ = adaptive_instance.interdiction_delta_;
  big_m_ = adaptive_instance.big_m_;

  arc_costs_ = std::vector<std::vector<int>>(scenarios_);
  scenario_index_map_ = std::vector<int>(scenarios_);

  for (int q = 0; q < scenarios_; q++) {
    arc_costs_[q] = adaptive_instance.arc_costs_[keep_scenarios[q]];
    scenario_index_map_[q] = keep_scenarios[q];
  }
}

void AdaptiveInstance::ReadCosts() {
  // Read arc costs from a file.
  // Compute interdiction_delta_, big_m_. (Maximum cost over all followers and
  // arcs, times the number of nodes).
  std::string line, word;
  std::ifstream myfile(costs_filename_);
  int q = 0;
  std::vector<int> v;
  int cost;
  int max_cost = INT_MIN;
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      v.clear();
      std::stringstream str(line);
      while (getline(str, word, ',')) {
        cost = std::stoi(word);
        if (cost > max_cost) max_cost = cost;
        v.push_back(cost);
      }
      if (q < scenarios_) {
        arc_costs_.push_back(v);
      }
      ++q;
    }
  }
  interdiction_delta_ = (max_cost)*nodes_;
  big_m_ = ((max_cost)*nodes_) + 1;
}

void AdaptiveInstance::PrintInstance(const Graph& G) const {
  // Print Summary of Problem Instance
  std::cout << "Followers/Scenarios: " << scenarios_ << std::endl;
  G.PrintGraphWithCosts(arc_costs_, interdiction_delta_);
}

double AdaptiveInstance::SPModel(int q, const Graph& G, GRBEnv* env) const {
  // Return shortest path objective function for follower q on instance.
  try {
    // Initialize Model/Environment.
    GRBModel* sp_model = new GRBModel(env);
    sp_model->set(GRB_IntParam_OutputFlag, 0);
    // Add decision variables.
    std::vector<GRBVar> y_var(G.arcs());
    for (int a = 0; a < G.arcs(); ++a) {
      std::string varname = "y_" + std::to_string(a);
      y_var[a] =
          sp_model->addVar(0, 1, arc_costs_[q][a], GRB_CONTINUOUS, varname);
    }
    // Add flow constraints.
    for (int i = 0; i < nodes_; ++i) {
      if (G.arc_index_hash()[i].empty() && G.rev_arc_index_hash()[i].empty())
        continue;
      GRBLinExpr lhs = 0;
      int rhs = 0;
      if (i == 0)
        rhs = 1;
      else if (i == nodes_ - 1)
        rhs = -1;
      for (size_t index = 0; index < G.arc_index_hash()[i].size(); index++) {
        int a = G.arc_index_hash()[i][index];
        lhs += (1) * y_var[a];
      }
      for (size_t index = 0; index < G.rev_arc_index_hash()[i].size();
           index++) {
        int a = G.rev_arc_index_hash()[i][index];
        lhs += (-1) * y_var[a];
      }
      sp_model->addConstr(lhs == rhs);
    }
    sp_model->update();
    // std::string modelname =
    //     "spmodel" + std::to_string(w) + "_" + std::to_string(q) + ".lp";
    // sp_model->write(modelname);
    sp_model->optimize();
    double obj = sp_model->get(GRB_DoubleAttr_ObjVal);
    delete sp_model;
    return obj;
  } catch (GRBException e) {
    std::cout << "Gurobi error number [AdaptiveInstance::SPModel]: "
              << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  } catch (...) {
    std::cout << "Non-gurobi error during optimization "
                 "[AdaptiveInstance::SPModel]"
              << std::endl;
    return -1;
  }
}

double AdaptiveInstance::SolveASPIZeroPolicies(const Graph& G,
                                               GRBEnv* env) const {
  double worst_case_objective = DBL_MAX;
  for (int q = 0; q < scenarios_; q++) {
    double obj = SPModel(q, G, env);
    if (obj < worst_case_objective) worst_case_objective = obj;
  }
  return worst_case_objective;
}

void AdaptiveInstance::ApplyInterdiction(const std::vector<double>& x_bar,
                                         bool rev) {
  // Receive interdiction policy and update costs for the adaptive instance.
  // If rev, we are "removing" the interdiction policy and returning the
  // instance to its original state.
  int arcs = x_bar.size();
  for (int a = 0; a < arcs; ++a) {
    if (x_bar[a] == 1) {
      for (int q = 0; q < scenarios_; ++q) {
        if (rev) {
          arc_costs_[q][a] -= interdiction_delta_;
        } else {
          arc_costs_[q][a] += interdiction_delta_;
        }
      }
    }
  }
}

double ComputeObjectiveOfPolicyForScenario(
    const ProblemInput& problem, const std::vector<double>& binary_policy,
    int q) {
  AdaptiveInstance instance_copy = problem.instance_;
  instance_copy.ApplyInterdiction(binary_policy);
  double objective = instance_copy.SPModel(q, problem.G_, problem.env_);
  return objective;
}

double AdaptiveInstance::ValidatePolicy(std::vector<double>& x_bar,
                                        const Graph& G, GRBEnv* env) {
  // Solve shortest path on interdicted graph - check objectives match.
  double objective = DBL_MAX;
  double sp_result;
  ApplyInterdiction(x_bar);
  for (int q = 0; q < scenarios_; ++q) {
    sp_result = SPModel(q, G, env);
    if (sp_result < objective) {
      objective = sp_result;
    }
  }
  ApplyInterdiction(x_bar, true);
  return objective;
}

void SetPartitioningModel::AddSetPartitioningVariables() {
  // Set partitioning variables.
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; ++q) {
      std::string varname = std::string("H_")
                                .append(std::to_string(w))
                                .append("_")
                                .append(std::to_string(q));
      h_var_[w][q] = sp_model_->addVar(0, 1, 0, GRB_BINARY, varname);
    }
  }
}

void SetPartitioningModel::AddInterdictionPolicyVariables() {
  // Interdiction policy on arcs.
  for (int w = 0; w < policies_; ++w) {
    for (int a = 0; a < arcs_; a++) {
      std::string varname = std::string("x_")
                                .append(std::to_string(w))
                                .append("_")
                                .append(std::to_string(a));
      x_var_[w][a] = sp_model_->addVar(0, 1, 0, GRB_BINARY, varname);
    }
  }
}

void SetPartitioningModel::AddPiShortestPathVariable() {
  // Pi shortest path variable.
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; q++) {
      for (int i = 0; i < nodes_; i++) {
        std::string varname = std::string("pi_")
                                  .append(std::to_string(w))
                                  .append("_")
                                  .append(std::to_string(q))
                                  .append("_")
                                  .append(std::to_string(i));
        pi_var_[w][q][i] = sp_model_->addVar(-GRB_INFINITY, GRB_INFINITY, 0,
                                             GRB_CONTINUOUS, varname);
      }
    }
  }
}

void SetPartitioningModel::AddDualLambdaArcVariable() {
  // arc variable 'lambda'
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; q++) {
      for (int a = 0; a < arcs_; a++) {
        std::string varname = std::string("lambda_")
                                  .append(std::to_string(w))
                                  .append("_")
                                  .append(std::to_string(q))
                                  .append("_")
                                  .append(std::to_string(a));
        lambda_var_[w][q][a] =
            sp_model_->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, varname);
      }
    }
  }
}

void SetPartitioningModel::AddBudgetConstraints() {
  for (int w = 0; w < policies_; ++w) {
    GRBLinExpr linexpr = 0;
    for (int a = 0; a < arcs_; a++) {
      linexpr += x_var_[w][a];
    }
    sp_model_->addConstr(linexpr <= budget_);
  }
}

void SetPartitioningModel::AddObjectiveBoundingConstraint() {
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; q++) {
      GRBLinExpr linexpr = 0;
      linexpr += (pi_var_[w][q][nodes_ - 1] -
                  pi_var_[w][q][0]);  // b^\top pi (our b is simply a
                                      // source-sink single unit of flow)
      for (int a = 0; a < arcs_; ++a) {
        linexpr += -lambda_var_[w][q][a];  // u^\top \cdot lambda (our u is 1)
      }

      linexpr += big_m_ * (1 - h_var_[w][q]);
      sp_model_->addConstr(z_var_ <= linexpr);
    }
  }
}

void SetPartitioningModel::AddArcCostBoundingConstraint(
    const ProblemInput& problem) {
  int j_node, arc;
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; ++q) {
      for (int i = 0; i < nodes_; ++i) {
        for (size_t j = 0; j < problem.G_.adjacency_list()[i].size(); ++j) {
          j_node = problem.G_.adjacency_list()[i][j];
          arc = problem.G_.arc_index_hash()[i][j];
          // add constraint
          sp_model_->addConstr(
              (pi_var_[w][q][j_node] - pi_var_[w][q][i] -
               lambda_var_[w][q][arc]) <=
              problem.instance_.arc_costs()[q][arc] +
                  (problem.instance_.interdiction_delta() * x_var_[w][arc]));
        }
      }
    }
  }
}

void SetPartitioningModel::AddSetPartitioningConstraints() {
  for (int q = 0; q < scenarios_; ++q) {
    GRBLinExpr linexpr = 0;
    for (int w = 0; w < policies_; ++w) {
      linexpr += h_var_[w][q];
    }
    sp_model_->addConstr(linexpr == 1);
  }
}

void SetPartitioningModel::AddRootNodeZeroConstraint() {
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; q++) {
      sp_model_->addConstr(pi_var_[w][q][0] == 0);
    }
  }
}

void SetPartitioningModel::AddAssignmentSymmetryConstraints() {
  for (int q = 0; q < policies_; q++) {
    GRBLinExpr linexpr = 0;
    for (int w = 0; w <= q; w++) {
      linexpr += h_var_[w][q];
    }
    sp_model_->addConstr(linexpr == 1);
  }
}

void SetPartitioningModel::AddNonDecreasingSymmetryConstraints() {
  for (int w = 0; w < policies_ - 1; w++) {
    int v = w + 1;  // "Next" cluster.
    GRBLinExpr lhs = 0;
    GRBLinExpr rhs = 0;
    for (int q = 0; q < scenarios_; q++) {
      lhs += h_var_[w][q];
      rhs += h_var_[v][q];
    }
    sp_model_->addConstr(lhs <= rhs);
  }
}

void SetPartitioningModel::ProcessInputSymmetryParameters(
    const ProblemInput& problem) {
  // Gurobi Symmetry parameters.
  if (problem.gurobi_symmetry_detection_ != GUROBI_SYMMETRY_AUTO) {
    sp_model_->set(GRB_IntParam_Symmetry, problem.gurobi_symmetry_detection_);
  }
  // Manual Symmetry constraints.
  if (problem.manual_symmetry_constraints_ == MANUAL_SYMMETRY_ASSIGNMENT) {
    AddAssignmentSymmetryConstraints();
  }
  if (problem.manual_symmetry_constraints_ == MANUAL_SYMMETRY_NONDECREASING) {
    AddNonDecreasingSymmetryConstraints();
  }
  // if problem.manual_symmetry_constraints_ == MANUAL_SYMMETRY_NONE: do
  // nothing.
}

void SetPartitioningModel::ConfigureSolver(const ProblemInput& problem) {
  try {
    sp_model_->set(GRB_IntParam_OutputFlag, 0);
    // Add model variables.
    AddSetPartitioningVariables();
    AddInterdictionPolicyVariables();
    AddPiShortestPathVariable();
    AddDualLambdaArcVariable();
    // Add model constraints.
    AddBudgetConstraints();
    AddObjectiveBoundingConstraint();
    AddArcCostBoundingConstraint(problem);
    AddSetPartitioningConstraints();
    AddRootNodeZeroConstraint();
    ProcessInputSymmetryParameters(problem);
    sp_model_->update();
  } catch (GRBException e) {
    std::cout << "Gurobi error number [SetPartitioningModel::ConfigureSolver]: "
              << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch (...) {
    std::cout << "Non-gurobi error during optimization "
                 "[SetPartitioningModel::ConfigureSolver]"
              << std::endl;
  }
}

AdaptiveSolution SetPartitioningModel::Solve(const ProblemInput& problem) {
  // Updates current solution to latest solution, or sets unbounded to true.
  // Compute full objective matrix after solving.
  // Returns current solution.
  //
  // Optimize model and measure solution time.
  long begin = GetCurrentTime();
  sp_model_->optimize();
  long time = GetCurrentTime() - begin;
  // Initialize AdaptiveSolution.
  AdaptiveSolution final_solution(false, false, problem);
  final_solution.set_solution_time(time);
  if (sp_model_->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
    // If status is optimal, set unbounded to false and set objective value,
    // policies and partition.
    final_solution.set_unbounded(false);
    final_solution.set_optimal(true);
    final_solution.set_worst_case_objective(
        -sp_model_->get(GRB_DoubleAttr_ObjVal));
    for (int w = 0; w < policies_; ++w) {
      std::vector<double> x_vector;
      for (int a = 0; a < arcs_; a++) {
        x_vector.push_back(x_var_[w][a].get(GRB_DoubleAttr_X));
      }
      final_solution.set_solution_policy(w, x_vector);
      for (int q = 0; q < scenarios_; q++) {
        if (h_var_[w][q].get(GRB_DoubleAttr_X) > 0.5) {
          final_solution.add_to_partition(w, q);
        }
      }
    }
    // Compute Objective Matrix.
    // Skip this as it is time consuming and not needed, only uncomment
    // for debugging.
    // final_solution.ComputeObjectiveMatrix(problem);
  } else if (sp_model_->get(GRB_IntAttr_Status) == GRB_UNBOUNDED) {
    final_solution.set_unbounded(true);
  } else if (sp_model_->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
    // If solution status is hit time limit, set time to TIME_LIMIT_MS.
    // Use milliseconds, as our adaptive solution class stores these values in
    // milliseconds, since GetCurrentTime() returns milliseconds.
    final_solution.set_unbounded(false);
    final_solution.set_solution_time(TIME_LIMIT_MS);
    // Check if there is an incumbent - in this case set solution and objective.
    if (sp_model_->get(GRB_IntAttr_SolCount) > 0) {
      final_solution.set_mip_gap(sp_model_->get(GRB_DoubleAttr_MIPGap));
      final_solution.set_worst_case_objective(
          -sp_model_->get(GRB_DoubleAttr_ObjVal));
      for (int w = 0; w < policies_; ++w) {
        std::vector<double> x_vector;
        for (int a = 0; a < arcs_; a++) {
          x_vector.push_back(x_var_[w][a].get(GRB_DoubleAttr_X));
        }
        final_solution.set_solution_policy(w, x_vector);
        for (int q = 0; q < scenarios_; q++) {
          if (h_var_[w][q].get(GRB_DoubleAttr_X) > 0.5) {
            final_solution.add_to_partition(w, q);
          }
        }
      }
    }
  }
  return final_solution;
}

void BendersCallback::ConfigureIndividualSubModel(const Graph& G, int w,
                                                  int q) {
  submodels_[w][q]->set(GRB_IntParam_OutputFlag, 0);
  // Add decision variables.
  for (int a = 0; a < arcs_; ++a) {
    std::string varname = "y_" + std::to_string(a);
    y_var_[w][q].push_back(submodels_[w][q]->addVar(0, 1, arc_costs_[q][a],
                                                    GRB_CONTINUOUS, varname));
  }
  // Add flow constraints.
  for (int i = 0; i < nodes_; ++i) {
    if (G.arc_index_hash()[i].empty() && G.rev_arc_index_hash()[i].empty())
      continue;
    GRBLinExpr lhs = 0;
    int rhs = 0;
    if (i == 0)
      rhs = 1;
    else if (i == nodes_ - 1)
      rhs = -1;
    for (size_t index = 0; index < G.arc_index_hash()[i].size(); index++) {
      int a = G.arc_index_hash()[i][index];
      lhs += (1) * y_var_[w][q][a];
    }
    for (size_t index = 0; index < G.rev_arc_index_hash()[i].size(); index++) {
      int a = G.rev_arc_index_hash()[i][index];
      lhs += (-1) * y_var_[w][q][a];
    }
    submodels_[w][q]->addConstr(lhs == rhs);
  }
}

void BendersCallback::ConfigureSubModels(const ProblemInput& problem) {
  // Initialize models (k x p).
  for (int w = 0; w < policies_; ++w) {
    submodels_.push_back(std::vector<GRBModel*>());
    for (int q = 0; q < scenarios_; ++q) {
      submodels_[w].push_back(new GRBModel(*env_));
    }
  }
  y_var_ = std::vector<std::vector<std::vector<GRBVar>>>(
      policies_, std::vector<std::vector<GRBVar>>(scenarios_));
  // Configure each individual submodel.
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; ++q) {
      ConfigureIndividualSubModel(problem.G_, w, q);
    }
  }
}

void BendersCallback::UpdateSubModels(bool rev) {
  for (int w = 0; w < policies_; ++w) {
    for (int a = 0; a < arcs_; ++a) {
      if (getSolution(x_var_[w][a]) > 0.5) {
        for (int q = 0; q < scenarios_; ++q) {
          if (rev) {
            int new_cost = arc_costs_[q][a];
            y_var_[w][q][a].set(GRB_DoubleAttr_Obj, new_cost);
          } else {
            int new_cost = arc_costs_[q][a] + interdiction_delta_;
            y_var_[w][q][a].set(GRB_DoubleAttr_Obj, new_cost);
          }
        }
      }
    }
  }
}

void BendersCallback::SolveSubModels() {
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; ++q) {
      submodels_[w][q]->optimize();
    }
  }
}

void BendersCallback::AddLazyCuts() {
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; ++q) {
      GRBLinExpr lhs = z_var_ - big_m_ * (1 - h_var_[w][q]);
      GRBLinExpr rhs = 0;
      for (int a = 0; a < arcs_; ++a) {
        if (y_var_[w][q][a].get(GRB_DoubleAttr_X) > 0.5) {
          rhs += arc_costs_[q][a] + interdiction_delta_ * x_var_[w][a];
        }
      }
      addLazy(lhs <= rhs);
    }
  }
  lazy_cuts_rounds_++;
}

void BendersCallback::callback() {
  if (where == GRB_CB_MIPSOL) {
    iteration_++;
    // Updating upper_bound_ here, as this is equivalent to the "end" of the
    // loop iteration in the paper algorithm (i.e. gurobi just re-solved the
    // master problem, which happens at the end of an iteration in the
    // algorithm, right after we do the subproblems and add cuts, and right
    // before we update the upper bound.
    upper_bound_ = -getDoubleInfo(GRB_CB_MIPSOL_OBJBST);
    if (upper_bound_ - lower_bound_ >= epsilon_) {
      UpdateSubModels();
      SolveSubModels();
      UpdateSubModels(true);
      // temp_bound is equal to:
      // min_q max_w submodel_objective[w][q]
      double temp_bound = DBL_MAX;
      for (int q = 0; q < scenarios_; ++q) {
        double incumbent = DBL_MIN;
        for (int w = 0; w < policies_; ++w) {
          double obj = submodels_[w][q]->get(GRB_DoubleAttr_ObjVal);
          if (obj > incumbent) incumbent = obj;
        }
        if (incumbent < temp_bound) temp_bound = incumbent;
      }
      if (lower_bound_ < temp_bound) {
        // Update lower bound.
        lower_bound_ = temp_bound;
        // Set current_solution_/x_prime (final interdiction solution) to
        // solution x_var_ that was used for the subproblems.
        // for (int w = 0; w < policies_; ++w) {
        //   std::vector<double> binary_policy(arcs_);
        //   for (int a = 0; a < arcs_; ++a) {
        //     if (getSolution(x_var_[w][a]) > 0.5)
        //       binary_policy[a] = 1;
        //     else
        //       binary_policy[a] = 0;
        //   }
        // }
      }
      AddLazyCuts();
    }
  }
}

void SetPartitioningBenders::AddSetPartitioningVariables() {
  // Set partitioning variables.
  for (int w = 0; w < policies_; ++w) {
    for (int q = 0; q < scenarios_; ++q) {
      std::string varname = std::string("H_")
                                .append(std::to_string(w))
                                .append("_")
                                .append(std::to_string(q));
      h_var_[w][q] = benders_model_->addVar(0, 1, 0, GRB_BINARY, varname);
    }
  }
}

void SetPartitioningBenders::AddInterdictionPolicyVariables() {
  // Interdiction policy on arcs.
  for (int w = 0; w < policies_; ++w) {
    for (int a = 0; a < arcs_; a++) {
      std::string varname = std::string("x_")
                                .append(std::to_string(w))
                                .append("_")
                                .append(std::to_string(a));
      x_var_[w][a] = benders_model_->addVar(0, 1, 0, GRB_BINARY, varname);
    }
  }
}

void SetPartitioningBenders::AddBudgetConstraints() {
  for (int w = 0; w < policies_; ++w) {
    GRBLinExpr linexpr = 0;
    for (int a = 0; a < arcs_; a++) {
      linexpr += x_var_[w][a];
    }
    benders_model_->addConstr(linexpr <= budget_);
  }
}

void SetPartitioningBenders::AddSetPartitioningConstraints() {
  for (int q = 0; q < scenarios_; ++q) {
    GRBLinExpr linexpr = 0;
    for (int w = 0; w < policies_; ++w) {
      linexpr += h_var_[w][q];
    }
    benders_model_->addConstr(linexpr == 1);
  }
}

void SetPartitioningBenders::ConfigureSolver(const ProblemInput& problem) {
  try {
    benders_model_->set(GRB_IntParam_OutputFlag, 0);
    benders_model_->getEnv().set(GRB_IntParam_LazyConstraints, 1);
    AddSetPartitioningVariables();
    AddInterdictionPolicyVariables();
    AddBudgetConstraints();
    AddSetPartitioningConstraints();
    callback_ = BendersCallback(problem, z_var_, h_var_, x_var_);
    callback_.ConfigureSubModels(problem);
  } catch (GRBException e) {
    std::cout
        << "Gurobi error number [SetPartitioningBenders::ConfigureSolver]: "
        << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch (...) {
    std::cout << "Non-gurobi error during optimization "
                 "[SetPartitioningBenders::ConfigureSolver]"
              << std::endl;
  }
}

AdaptiveSolution SetPartitioningBenders::Solve(const ProblemInput& problem) {
  // Set callback, optimize and measure solution time.
  benders_model_->setCallback(&callback_);
  long begin = GetCurrentTime();
  benders_model_->optimize();
  long time = GetCurrentTime() - begin;
  AdaptiveSolution final_solution(true, false, problem);
  if (benders_model_->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
    // If optimal, set solution time objective and policy.
    final_solution.set_unbounded(false);
    final_solution.set_optimal(true);
    final_solution.set_solution_time(time);
    final_solution.set_worst_case_objective(
        -benders_model_->get(GRB_DoubleAttr_ObjVal));
    final_solution.set_cuts(callback_.lazy_cuts_rounds());
    for (int w = 0; w < policies_; ++w) {
      std::vector<double> x_vector;
      for (int a = 0; a < arcs_; a++) {
        x_vector.push_back(x_var_[w][a].get(GRB_DoubleAttr_X));
      }
      final_solution.set_solution_policy(w, x_vector);
      for (int q = 0; q < scenarios_; q++) {
        if (h_var_[w][q].get(GRB_DoubleAttr_X) > 0.5) {
          final_solution.add_to_partition(w, q);
        }
      }
    }
    // Skip this as it is time consuming and not needed, only uncomment
    // for debugging.
    // final_solution.ComputeObjectiveMatrix(problem);
  } else if (benders_model_->get(GRB_IntAttr_Status) == GRB_UNBOUNDED) {
    final_solution.set_unbounded(true);
    final_solution.set_solution_time(time);
  } else if (benders_model_->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
    // If we hit time limit, check for incumbent and set policy and objective.
    // If solution status is hit time limit, set time to TIME_LIMIT_MS.
    // Use milliseconds, as our adaptive solution class stores these values in
    // milliseconds, since GetCurrentTime() returns milliseconds.
    final_solution.set_unbounded(false);
    final_solution.set_solution_time(TIME_LIMIT_MS);
    if (benders_model_->get(GRB_IntAttr_SolCount) > 0) {
      final_solution.set_mip_gap(benders_model_->get(GRB_DoubleAttr_MIPGap));
      final_solution.set_cuts(callback_.lazy_cuts_rounds());
      final_solution.set_worst_case_objective(
          -benders_model_->get(GRB_DoubleAttr_ObjVal));
      for (int w = 0; w < policies_; ++w) {
        std::vector<double> x_vector;
        for (int a = 0; a < arcs_; a++) {
          x_vector.push_back(x_var_[w][a].get(GRB_DoubleAttr_X));
        }
        final_solution.set_solution_policy(w, x_vector);
        for (int q = 0; q < scenarios_; q++) {
          if (h_var_[w][q].get(GRB_DoubleAttr_X) > 0.5) {
            final_solution.add_to_partition(w, q);
          }
        }
      }
    }
    // Skip this as it is time consuming and not needed, only uncomment
    // for debugging.
    // final_solution.ComputeObjectiveMatrix(problem);
  }
  return final_solution;
}

std::vector<int> InitKappa(int p, int k) {
  // Initialize partition vector based on total number in set (p) and exact
  // number of partitions required.
  std::vector<int> kappa(p, 0);
  for (int q = (p - k + 1); q < p; ++q) {
    kappa[q] = (q - (p - k));
  }
  return kappa;
}

void RobustAlgoModel::ConfigureModel(const ProblemInput& problem) {
  // Construct baseline model with no constraints.
  algo_model_ = new GRBModel(env_);
  algo_model_->set(GRB_IntParam_OutputFlag, 0);
  // Decision Variables.
  std::string varname = "z";
  z_var_ = algo_model_->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS,
                               varname);  // objective function dummy variable.
  for (int a = 0; a < arcs_; ++a) {
    varname = "x_" + std::to_string(a);
    x_var_.push_back(algo_model_->addVar(0, 1, 0, GRB_BINARY, varname));
  }  // interdiction policy on arcs.
  std::vector<GRBVar> tempvector;
  for (int q = 0; q < scenarios_; ++q) {
    pi_var_.push_back(
        tempvector);  // post interdiction s-i best path (for every q).
    for (int i = 0; i < nodes_; ++i) {
      varname = "pi_" + std::to_string(q) + "_" + std::to_string(i);
      pi_var_[q].push_back(algo_model_->addVar(-GRB_INFINITY, GRB_INFINITY, 0,
                                               GRB_CONTINUOUS, varname));
    }
  }
  for (int q = 0; q < scenarios_; ++q) {
    lambda_var_.push_back(
        tempvector);  // lambda variable on arcs (for every q).
    for (int a = 0; a < arcs_; ++a) {
      varname = "lambda_" + std::to_string(q) + "_" + std::to_string(a);
      lambda_var_[q].push_back(
          algo_model_->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
    }
  }
  // Add Budget Constraint.
  GRBLinExpr linexpr = 0;
  for (int a = 0; a < arcs_; a++) {
    linexpr += x_var_[a];
  }
  algo_model_->addConstr(linexpr <= budget_, "budget");
  algo_model_->update();
  // Arc/Dual Constraints.
  for (int q = 0; q < scenarios_; ++q) {
    dual_constraints_.push_back(std::vector<GRBTempConstr>());
    for (int i = 0; i < nodes_; ++i) {
      for (size_t j = 0; j < problem.G_.adjacency_list()[i].size(); ++j) {
        int next = problem.G_.adjacency_list()[i][j];
        int a = problem.G_.arc_index_hash()[i][j];
        GRBTempConstr constraint =
            pi_var_[q][next] - pi_var_[q][i] - lambda_var_[q][a] <=
            problem.instance_.arc_costs()[q][a] +
                (problem.instance_.interdiction_delta() * x_var_[a]);
        dual_constraints_[q].push_back(constraint);
      }
    }
  }
  // Objective Constraints.
  for (int q = 0; q < scenarios_; ++q) {
    GRBLinExpr linexpr = 0;
    linexpr += (pi_var_[q][nodes_ - 1] - pi_var_[q][0]);  // b^\top pi.
    for (int a = 0; a < arcs_; ++a) {
      linexpr += -lambda_var_[q][a];  // u^\top \cdot lambda.
    }
    z_constraints_.push_back(z_var_ <= linexpr);
  }
}

void RobustAlgoModel::Update(std::vector<int>& subset) {
  // Add constraints to model for subset in partition.
  for (int q : subset) {
    std::string zero_name = "zero_" + std::to_string(q);
    std::string z_name = "z_" + std::to_string(q);
    algo_model_->addConstr(pi_var_[q][0] == 0, zero_name);
    algo_model_->addConstr(z_constraints_[q], z_name);
    int a = 0;
    for (GRBTempConstr constraint : dual_constraints_[q]) {
      std::string dual_name =
          "dual_" + std::to_string(q) + "_" + std::to_string(a);
      algo_model_->addConstr(constraint, dual_name);
      ++a;
    }
  }
  algo_model_->update();
}

void RobustAlgoModel::ReverseUpdate(std::vector<int>& subset) {
  // Remove all constraints from model - i.e. all constraints associated with
  // this subset.
  for (int q : subset) {
    std::string zero_name = "zero_" + std::to_string(q);
    std::string z_name = "z_" + std::to_string(q);
    GRBConstr zero_constraint = algo_model_->getConstrByName(zero_name);
    GRBConstr z_constraint = algo_model_->getConstrByName(z_name);
    algo_model_->remove(zero_constraint);
    algo_model_->remove(z_constraint);
    int a = 0;
    for (GRBTempConstr constraint : dual_constraints_[q]) {
      std::string dual_name =
          "dual_" + std::to_string(q) + "_" + std::to_string(a);
      GRBConstr dual_constraint = algo_model_->getConstrByName(dual_name);
      algo_model_->remove(dual_constraint);
      ++a;
    }
  }
  algo_model_->update();
}

std::pair<double, std::vector<double>> RobustAlgoModel::Solve() {
  // Solve static model, return objective value and policy.
  // If it hits the time limit, return {-1, {-1,....-1}}
  // If there is any status other than time limit or optimal,
  // return {0, {0,....,0}}.
  algo_model_->optimize();
  if (algo_model_->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
    return std::make_pair(-1, std::vector<double>(arcs_, -1));
  } else if (algo_model_->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
    std::vector<double> binary_policy(arcs_, 0);
    double objective = -algo_model_->get(GRB_DoubleAttr_ObjVal);
    for (int a = 0; a < arcs_; ++a) {
      binary_policy[a] = x_var_[a].get(GRB_DoubleAttr_X);
    }
    return std::make_pair(objective, binary_policy);
  } else
    return std::make_pair(0, std::vector<double>(arcs_, 0));
  // This is not necessary and we might even lose speed on warm starts etc.
  // algo_model_->reset();
}

int max_int(int a, int b) {
  if (a <= b)
    return b;
  else
    return a;
}

void AdaptiveSolution::ComputeObjectiveMatrix(const ProblemInput& problem) {
  // Populate / recompute the objectives_ matrix.
  for (int w = 0; w < policies_; w++) {
    for (int q = 0; q < scenarios_; q++) {
      double objective =
          ComputeObjectiveOfPolicyForScenario(problem, solution_[w], q);
      objectives_[w][q] = objective;
    }
  }
}

void AdaptiveSolution::ComputePartition() {
  // Set up the partition when it is uninitialized (only needed for the
  // approximation algorithm). Objective matrix should be populated at this
  // point - i.e. ComputeObjectiveMatrix should be called first. For each
  // scenario, find the policy w that maximizes their objectives (all k will be
  // available in objectives_).
  for (int q = 0; q < scenarios_; ++q) {
    double max_objective = DBL_MIN;
    int partition_index = -1;
    for (int w = 0; w < policies_; ++w) {
      double objective_value = objectives_[w][q];
      if (objective_value > max_objective) {
        max_objective = objective_value;
        partition_index = w;
      }
    }
    partition_[partition_index].push_back(q);
  }
}

void AdaptiveSolution::ComputeAdaptiveObjective() {
  // Minimizing - finding the worst (minimum) objective value, out of all of the
  // follower's objectives (using the best (maximum) policy for each follower).
  // Only to be run once ComputeAllObjectives is run.
  // Only reassigns worst_case_objective_ if it was -1 - otherwise, it does
  // nothing, but will output a warning if the computed objective and previous
  // objective are different by more than epsilon. If unbounded_ was set to
  // true, the function is skipped all together.
  if (unbounded_) return;
  double adaptive_objective = DBL_MAX;
  for (int q = 0; q < scenarios_; ++q) {
    double follower_objective = DBL_MIN;
    for (int w = 0; w < policies_; ++w) {
      if (objectives_[w][q] > follower_objective)
        follower_objective = objectives_[w][q];
    }
    if (follower_objective < adaptive_objective)
      adaptive_objective = follower_objective;
  }
  if (worst_case_objective_ == -1)
    worst_case_objective_ = adaptive_objective;
  else {
    if (std::abs(worst_case_objective_ - adaptive_objective) > EPSILON) {
      std::cout << "WARNING[AdaptiveSolution::ComputeAdaptiveObjective]: "
                << worst_case_objective_ << " (saved), " << adaptive_objective
                << " (computed)" << std::endl;
    }
  }
}

void AdaptiveSolution::LogSolution(const ProblemInput& problem, bool policy) {
  std::string title = problem.instance_.name();
  std::cout << std::endl << std::endl;
  if (title == "") {
    std::cout << "Solution, time (ms): " << solution_time_ << " ms";
  } else {
    std::cout << title << ", time: " << solution_time_ << " ms";
  }
  std::cout << ", k = " << policies_ << ", p = " << scenarios_
            << ", r_0 = " << budget_;
  std::cout << ", nodes = " << nodes_ << ", arcs = " << arcs_ << std::endl;
  std::cout << "Worst Case Objective: " << worst_case_objective_ << std::endl;
  if (policy) {
    std::cout << "Objective Matrix: " << std::endl;
    for (int w = 0; w < policies_; ++w) {
      std::cout << "Policy " << w << ": ";
      for (double objective : objectives_[w]) {
        std::cout << objective << " ";
      }
      std::cout << std::endl;
    }
  }
  for (int w = 0; w < policies_; ++w) {
    if (policy) {
      std::cout << "subset: { ";
      for (int q : partition_[w]) {
        std::cout << q << " ";
      }
      std::cout << "}, interdicted arcs: ";
      std::vector<std::vector<int>> arc_index_hash =
          problem.G_.arc_index_hash();
      std::vector<std::vector<int>> adjacency_list =
          problem.G_.adjacency_list();
      for (int i = 0; i < problem.G_.nodes(); ++i) {
        int index = 0;
        for (int a : arc_index_hash[i]) {
          if (solution_[w][a] > 0.5) {
            problem.G_.PrintArc(a, i, index);
            std::cout << " ";
          }
          index++;
        }
      }
    }
    std::cout << std::endl;
  }
  if (benders_)
    std::cout << "Rounds of lazy cuts: " << lazy_cuts_rounds_ << std::endl;
}

void AdaptiveSolution::MergeEnumSols(AdaptiveSolution sol2,
                                     AdaptiveInstance* instance2,
                                     int split_index) {
  // Merge 2 adaptive solutions
  // The AdaptiveSolution being passed has k=2
  // The AdaptiveSolution being worked on needs its split_index replaced:
  //   - in .partition, by the two subsets in sol2
  //   - in .solution_, by the two policies in sol2
  //   - in both cases, we can just replace split_index by the first, and add
  //   the second to the end
  // Additionally, .policies must be increased by 1
  // Remember - change-of-scenario copy constructor reindexes, but we maintain
  // a map of the original indices (sol2 has this map)
  set_policies(policies_ + 1);
  // reindexed copy of sol2.partition
  std::vector<std::vector<int>> reindexed_partition(sol2.policies());
  int w = 0;
  for (auto& subset : sol2.partition()) {
    for (int i : subset) {
      int q = instance2->scenario_index_map()[i];
      reindexed_partition[w].push_back(q);
    }
    ++w;
  }
  // replace and push_back in partition
  partition_[split_index] = reindexed_partition[0];
  for (int w = 1; w < sol2.policies(); w++) {
    partition_.push_back(reindexed_partition[w]);
  }
  // replace in solution_
  solution_[split_index] = sol2.solution()[0];
  for (int w = 1; w < sol2.policies(); w++) {
    solution_.push_back(sol2.solution()[w]);
  }
}

bool nextKappa(std::vector<int>& kappa, std::vector<int>& max, int k, int p) {
  // update kappa and max place to next partition
  // if this is the last one, return false

  for (int q = (p - 1); q > 0; --q) {
    if ((kappa[q] < k - 1) && (kappa[q] <= max[q - 1])) {
      ++kappa[q];
      max[q] = max_int(max[q], kappa[q]);
      for (int u = q + 1; u <= (p - (k - max[q])); ++u) {
        kappa[u] = 0;
        max[u] = max[q];
      }
      for (int u = (p - (k - max[q])) + 1; u < p; ++u) {
        kappa[u] = (k - (p - u));
        max[u] = kappa[u];
      }
      return true;
    }
  }
  return false;
}

std::vector<std::vector<int>> kappa_to_partition(std::vector<int>& kappa, int k,
                                                 int p) {
  // convert a kappa vector to a partition of subsets
  std::vector<int> subset;
  std::vector<std::vector<int>> partition(k, subset);

  for (int q = 0; q < p; ++q) {
    int subset_index = kappa[q];
    partition[subset_index].push_back(q);
  }

  return partition;
}

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>>
MergeEnumSols(std::pair<std::vector<std::vector<int>>,
                        std::vector<std::vector<double>>>& sol1,
              std::pair<std::vector<std::vector<int>>,
                        std::vector<std::vector<double>>>& sol2,
              int w_index) {
  // Merge 2 solutions e.g. when extending an enum solve, combine the new split
  // worst subset w_index is the position of the split subset in the original
  // solution (sol1)
  // place the first subset/policy of sol2 at position w_index of sol1,
  // replacing the subset that was split
  sol1.first[w_index] = sol2.first[0];
  sol1.second[w_index] = sol2.second[0];
  // place the second subset/policy of sol2 at the end of sol1, extending it to
  // k+1
  sol1.first.push_back(sol2.first[1]);
  sol1.second.push_back(sol2.second[1]);
  return sol1;
}

std::pair<double, std::vector<double>> SolveBendersInEnumSolve(
    ProblemInput& problem, std::vector<int>& update_vector) {
  // Solve the nominal problem for followers in update vector.
  ProblemInput sub_problem(problem, update_vector, 1);
  SetPartitioningBenders sub_benders(sub_problem);
  sub_benders.ConfigureSolver(sub_problem);
  AdaptiveSolution sol = sub_benders.Solve(sub_problem);
  if (sol.optimal())
    return {sol.worst_case_objective(), sol.solution()[0]};
  else
    return {-1, {-1}};
}

AdaptiveSolution EnumSolve(ProblemInput& problem) {
  // Initialize solution object for if we hit time limit.
  AdaptiveSolution time_limit_solution = AdaptiveSolution(false, false, false);
  time_limit_solution.set_solution_time(TIME_LIMIT_MS);
  time_limit_solution.set_worst_case_objective(-1);
  // Pass an AdaptiveInstance and solve using enumeration algorithm.
  // Enumeration maintained as in Orlov paper.
  bool next = true;
  int p = problem.instance_.scenarios();
  int k = problem.policies_;
  int m = problem.G_.arcs();
  // Initialize partitioning 'string' vector as in Orlov Paper.
  // Initialize corresponding max vector.
  std::vector<int> kappa = InitKappa(p, k);
  std::vector<int> max = kappa;
  // Initialize solution vector.
  std::vector<double> arc_vec(m, 0);
  std::vector<std::vector<double>> sol(k, arc_vec);
  try {
    long begin = GetCurrentTime();
    // Objective value maintained here (and optimal partition).
    double best_worstcase_objective = 0;
    std::vector<double> final_objectives(k, 0);
    std::vector<std::vector<int>> best_worstcase_partition;
    // Enumerate while not 'failing' to get next partition.
    while (next) {
      long current_time = GetCurrentTime() - begin;
      // Use milliseconds for time limit to check this, as we measured in
      // milliseconds with GetCurrentTime().
      if (current_time >= TIME_LIMIT_MS) {
        time_limit_solution.set_partition(best_worstcase_partition);
        time_limit_solution.set_worst_case_objective(best_worstcase_objective);
        return time_limit_solution;
      }
      std::vector<std::vector<double>> temp_sol(k, arc_vec);
      double temp_worst_objective = DBL_MAX;
      std::vector<double> temp_objectives(k, 0);
      std::vector<std::vector<int>> partition = kappa_to_partition(kappa, k, p);
      for (int w = 0; w < k; ++w) {
        // For every subset in partition, solve M2 for k=1 using the Benders.
        std::pair<double, std::vector<double>> temp_single_solution =
            SolveBendersInEnumSolve(problem, partition[w]);
        // Check that the static benders model didn't hit the time limit:
        if (temp_single_solution.first == -1) {
          time_limit_solution.set_partition(best_worstcase_partition);
          time_limit_solution.set_worst_case_objective(
              best_worstcase_objective);
          return time_limit_solution;
        }
        temp_objectives[w] = temp_single_solution.first;
        temp_sol[w] = temp_single_solution.second;
        // Update temp_worst_objective and temp_sol if this subset is worse.
        if (temp_objectives[w] < temp_worst_objective) {
          temp_worst_objective = temp_objectives[w];
        }
      }
      if (temp_worst_objective > best_worstcase_objective) {
        best_worstcase_objective = temp_worst_objective;
        sol = temp_sol;
        best_worstcase_partition = partition;
        final_objectives = temp_objectives;
      }
      next = nextKappa(kappa, max, k, p);
    }
    long time = GetCurrentTime() - begin;
    if (time >= TIME_LIMIT_MS) {
      time_limit_solution.set_partition(best_worstcase_partition);
      time_limit_solution.set_worst_case_objective(best_worstcase_objective);
      return time_limit_solution;
    }
    std::vector<std::vector<double>> final_policies(k);
    for (int w = 0; w < k; ++w) {
      final_policies[w] = sol[w];
      // final_policies[w].set_objective(final_objectives[w]);
    }
    // Benders == false, optimal == true.
    AdaptiveSolution final_solution = AdaptiveSolution(
        false, true, problem, best_worstcase_partition, final_policies);
    final_solution.set_worst_case_objective(best_worstcase_objective);
    final_solution.set_solution_time(time);
    // Skip this as it is time consuming and not needed, only uncomment
    // for debugging.
    // final_solution.ComputeObjectiveMatrix(problem);
    return final_solution;
  } catch (GRBException e) {
    std::cout << "Gurobi error number [EnumSolve]: " << e.getErrorCode()
              << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch (...) {
    std::cout << "Non Gurobi error during construction of static robust model "
                 "object [EnumSolve]"
              << std::endl;
  }
  AdaptiveSolution dummy_solution;
  return dummy_solution;
}

std::pair<double, std::vector<double>> SolveBendersInGreedyAlgorithm(
    ProblemInput& problem, int q) {
  // Solve the nominal problem for follower q.
  std::vector<int> Update_vector{q};  // Vector with follower q.
  ProblemInput sub_problem(problem, Update_vector, 1);
  SetPartitioningBenders sub_benders(sub_problem);
  sub_benders.ConfigureSolver(sub_problem);
  if (problem.greedy_mip_gap_threshold_ != -1) {
    sub_benders.SetMIPGap(problem.greedy_mip_gap_threshold_);
  }
  AdaptiveSolution sol = sub_benders.Solve(sub_problem);
  if (sol.optimal())
    return {sol.worst_case_objective(), sol.solution()[0]};
  else
    return {-1, {-1}};
}

int FirstFollowerInGreedyAlgorithm(ProblemInput& problem) {
  std::pair<int, double> min_subset = {-1, DBL_MAX};
  for (int q = 0; q < problem.instance_.scenarios(); q++) {
    double obj = problem.instance_.SPModel(q, problem.G_, problem.env_);
    if (obj < min_subset.second) min_subset = {q, obj};
  }
  return min_subset.first;
}

AdaptiveSolution GreedyAlgorithm(ProblemInput& problem) {
  // Solution in case we hit time limit.
  AdaptiveSolution time_limit_solution = AdaptiveSolution(false, false, false);
  time_limit_solution.set_solution_time(TIME_LIMIT_MS);
  time_limit_solution.set_worst_case_objective(-1);
  long begin = GetCurrentTime();
  std::vector<std::vector<double>> objective_matrix;
  std::unordered_set<int> centers;
  std::vector<std::vector<double>> final_policies;
  // Choose the first scenario to solve SPI for arbitrarily (we'll just use
  // index 0).
  int next_follower = FirstFollowerInGreedyAlgorithm(problem);
  centers.insert(next_follower);
  // Solve problem for next follower, and add a new vector to objective matrix,
  // i.e. the objective value of the new policy interdicting each follower.
  std::pair<double, std::vector<double>> temp_single_solution =
      SolveBendersInGreedyAlgorithm(problem, next_follower);
  if (temp_single_solution.first == -1) return time_limit_solution;
  std::vector<double> new_objectives(problem.instance_.scenarios());
  for (int q = 0; q < problem.instance_.scenarios(); q++) {
    new_objectives[q] = ComputeObjectiveOfPolicyForScenario(
        problem, temp_single_solution.second, q);
  }
  objective_matrix.push_back(new_objectives);
  final_policies.push_back(temp_single_solution.second);
  size_t policies = problem.policies_;
  while (centers.size() < policies) {
    long current_time = GetCurrentTime() - begin;
    // Use milliseconds since we used GetCurrentTime, which is in milliseconds,
    // to measure.
    if (current_time >= TIME_LIMIT_MS) {
      AdaptiveSolution time_limit_solution =
          AdaptiveSolution(false, false, false);
      time_limit_solution.set_solution_time(TIME_LIMIT_MS);
      time_limit_solution.set_worst_case_objective(-1);
      return time_limit_solution;
    }
    // Evaluate all non-center scenarios against the current policies,
    // maintaining the minimum one.
    std::pair<int, double> min_subset = {-1, DBL_MAX};
    for (int q = 0; q < problem.instance_.scenarios(); q++) {
      if (centers.count(q) == 1) continue;
      double max_for_q = DBL_MIN;
      for (const auto& objective_row : objective_matrix) {
        if (objective_row[q] > max_for_q) max_for_q = objective_row[q];
      }
      if (max_for_q < min_subset.second) min_subset = {q, max_for_q};
    }
    next_follower = min_subset.first;
    // Solve SPI for scenario that is currently performing worst (for the
    // interdictor), and add the scenario to the set of centers.
    // std::pair<double, std::vector<double>> temp_single_solution =
    // spi_model.Solve();
    //
    temp_single_solution =
        SolveBendersInGreedyAlgorithm(problem, next_follower);
    if (temp_single_solution.first == -1) return time_limit_solution;
    std::vector<double> new_objectives(problem.instance_.scenarios());
    for (int q = 0; q < problem.instance_.scenarios(); q++) {
      new_objectives[q] = ComputeObjectiveOfPolicyForScenario(
          problem, temp_single_solution.second, q);
    }
    objective_matrix.push_back(new_objectives);
    final_policies.push_back(temp_single_solution.second);
    centers.insert(next_follower);
  }
  std::vector<std::vector<int>> temp_partition(problem.policies_,
                                               std::vector<int>(0));
  // Benders == false, optimal == false.
  AdaptiveSolution final_solution(false, false, problem, temp_partition,
                                  final_policies);
  long solution_time = GetCurrentTime() - begin;
  if (solution_time >= TIME_LIMIT_MS) return time_limit_solution;
  final_solution.SetObjectiveMatrix(objective_matrix);
  final_solution.ComputeAdaptiveObjective();
  final_solution.ComputePartition();
  final_solution.set_solution_time(solution_time);
  return final_solution;
}

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
  while ('0' <= name[i] && name[i] <= '9') {
    n_digits.push_back(name[i] - '0');
    i++;
  }
  i++;
  while ('0' <= name[i] && name[i] <= '9') {
    kzero_digits.push_back(name[i] - '0');
    i++;
  }
  return {DigitsToInt(n_digits), DigitsToInt(kzero_digits)};
}

std::pair<int, int> GetScenariosID(const std::string& name) {
  int i = name.size() - 5;
  std::vector<int> p_digits;
  std::vector<int> id_digits;
  while ('0' <= name[i] && name[i] <= '9') {
    id_digits.insert(id_digits.begin(), name[i] - '0');
    i--;
  }
  i--;
  while ('0' <= name[i] && name[i] <= '9') {
    p_digits.insert(p_digits.begin(), name[i] - '0');
    i--;
  }
  return {DigitsToInt(p_digits), DigitsToInt(id_digits)};
  return {};
}

std::string SolutionPartitionToString(const AdaptiveSolution& solution,
                                      const ProblemInput& problem) {
  // Turn partition into string for column encoding:
  // Partitions will be encoded as follows in the csv column:
  // For every follower (indexed 0-p-1), we write the index (0-k-1) of the
  // policy that they are assigned to, in a string delimeted by '-'. For
  // example, if we have p=3 and k=2, with followers 0,1 in partition 0, and
  // follower 2 in partition 1, we have: "0-0-1".
  std::vector<int> partition_vector(problem.instance_.scenarios());
  std::string partition_string = "";
  std::vector<std::vector<int>> partition = solution.partition();
  int policy = 0;
  for (std::vector<int>& subset : partition) {
    for (int q : subset) {
      partition_vector[q] = policy;
    }
    ++policy;
  }
  int q = 0;
  for (int w : partition_vector) {
    partition_string.append(std::to_string(w));
    if (q < problem.instance_.scenarios() - 1) partition_string.append("-");
    q++;
  }
  return partition_string;
}

std::string SolveAndPrintTest(const std::string& set_name,
                              const ProblemInput& problem,
                              ProblemInput& problem_copyable,
                              const std::vector<ASPI_Solver>& solvers,
                              int debug) {
  std::vector<double> adaptive_objectives;
  std::string final_csv_string = set_name;
  final_csv_string.append(",");
  final_csv_string.append(problem.instance_.name());
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.G_.nodes()));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.G_.arcs()));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.k_zero_));
  double nodes = problem.G_.nodes();
  double arcs = problem.G_.arcs();
  double density = (arcs) / ((nodes) * (nodes - 1));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(density));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.instance_.scenarios()));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.budget_));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.policies_));
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
      if (sp_solution.optimal()) {
        optimal = "OPTIMAL";
      } else {
        optimal = "NOT OPTIMAL";
      }
      final_csv_string.append(",");
      final_csv_string.append(optimal);
      final_csv_string.append(",");
      final_csv_string.append(
          std::to_string(sp_solution.worst_case_objective()));
      final_csv_string.append(",");
      final_csv_string.append(std::to_string(sp_solution.mip_gap()));
      final_csv_string.append(",");
      final_csv_string.append(std::to_string(sp_solution.solution_time()));
      final_csv_string.append(",");
      final_csv_string.append(SolutionPartitionToString(sp_solution, problem));
      log_line.append("MIP: ");
      log_line.append(optimal);
      log_line.append(" - ");
      log_line.append(std::to_string(sp_solution.worst_case_objective()));
      log_line.append(", ");
      log_line.append(std::to_string(sp_solution.solution_time()));
      log_line.append("ms ----- ");
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
      if (benders_solution.optimal()) {
        optimal = "OPTIMAL";
      } else {
        optimal = "NOT OPTIMAL";
      }
      final_csv_string.append(",");
      final_csv_string.append(optimal);
      final_csv_string.append(",");
      final_csv_string.append(
          std::to_string(benders_solution.worst_case_objective()));
      final_csv_string.append(",");
      final_csv_string.append(std::to_string(benders_solution.mip_gap()));
      final_csv_string.append(",");
      final_csv_string.append(std::to_string(benders_solution.solution_time()));
      final_csv_string.append(",");
      final_csv_string.append(
          std::to_string(benders_solution.lazy_cuts_rounds()));
      final_csv_string.append(",");
      final_csv_string.append(
          SolutionPartitionToString(benders_solution, problem));
      log_line.append("BENDERS: ");
      log_line.append(optimal);
      log_line.append(" - ");
      log_line.append(std::to_string(benders_solution.worst_case_objective()));
      log_line.append(",");
      log_line.append(std::to_string(benders_solution.solution_time()));
      log_line.append("ms ----- ");
    } else if (solver == ENUMERATION) {
      AdaptiveSolution enum_solution = EnumSolve(problem_copyable);
      if (debug == 2)
        enum_solution.LogSolution(problem, true);
      else if (debug == 1)
        enum_solution.LogSolution(problem, false);
      adaptive_objectives.push_back(enum_solution.worst_case_objective());
      std::string optimal;
      if (enum_solution.optimal()) {
        optimal = "OPTIMAL";
      } else {
        optimal = "NOT OPTIMAL";
      }
      final_csv_string.append(",");
      final_csv_string.append(optimal);
      final_csv_string.append(",");
      final_csv_string.append(
          std::to_string(enum_solution.worst_case_objective()));
      final_csv_string.append(",");
      final_csv_string.append(std::to_string(enum_solution.solution_time()));
      final_csv_string.append(",");
      final_csv_string.append(
          SolutionPartitionToString(enum_solution, problem));
      log_line.append("ENUMERATION: ");
      log_line.append(optimal);
      log_line.append(" - ");
      log_line.append(std::to_string(enum_solution.worst_case_objective()));
      log_line.append(", ");
      log_line.append(std::to_string(enum_solution.solution_time()));
      log_line.append("ms ----- ");
    } else if (solver == GREEDY) {
      AdaptiveSolution greedy_solution = GreedyAlgorithm(problem_copyable);
      if (debug == 2)
        greedy_solution.LogSolution(problem, true);
      else if (debug == 1)
        greedy_solution.LogSolution(problem, false);
      final_csv_string.append(",");
      final_csv_string.append(
          std::to_string(greedy_solution.worst_case_objective()));
      final_csv_string.append(",");
      final_csv_string.append(std::to_string(greedy_solution.solution_time()));
      final_csv_string.append(",");
      final_csv_string.append(
          SolutionPartitionToString(greedy_solution, problem));
      log_line.append("GREEDY: ");
      log_line.append(std::to_string(greedy_solution.worst_case_objective()));
      log_line.append(", ");
      log_line.append(std::to_string(greedy_solution.solution_time()));
      log_line.append("ms ---------- ");
      // adaptive_objectives.push_back(enum_solution.worst_case_objective());
    }
  }
  std::cout << log_line << std::endl;
  // FOR TESTING
  // double prev = adaptive_objectives[0];
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
  // std::cout << "------------------------  PASS, EXACT OBJECTIVE: " <<
  // adaptive_objectives[0] << ", GREEDY APPROXIMATION: " <<
  // greedy_adaptive_objective
  //           << "  ------------------------" << std::endl;
  return final_csv_string;
}

// TODO: Add consequence of symmetry parameters being passed, manual and gurobi.
void RunAllInstancesInSetDirectory(
    const int min_policies, const int max_policies, const int min_budget,
    const int max_budget, const std::string& set_name,
    const std::vector<ASPI_Solver>& solvers, int manual_symmetry_constraints,
    int gurobi_symmetry_detection, double greedy_mip_gap_threshold) {
  GRBEnv* env = new GRBEnv();  // Initialize global gurobi environment.
  // Use seconds, since gurobi takes the parameter value in seconds.
  env->set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_S);  // Set time limit.
  std::cout << std::endl;
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
  std::string base_name = DATA_DIRECTORY;
  base_name.append(set_name);
  base_name.append("/");
  base_name.append(set_name);
  base_name.append("_run_");
  base_name.append(today);
  std::string result_file_name = base_name;
  result_file_name.append(".csv");
  std::ofstream result_file(result_file_name);

  std::string CSV_HEADER = INSTANCE_INFO_COLUMN_HEADERS;
  for (auto& solver : solvers) {
    if (solver == MIP) {
      CSV_HEADER.append(",");
      CSV_HEADER.append(MIP_COLUMN_HEADERS);
    } else if (solver == BENDERS) {
      CSV_HEADER.append(",");
      CSV_HEADER.append(BENDERS_COLUMN_HEADERS);
    } else if (solver == ENUMERATION) {
      CSV_HEADER.append(",");
      CSV_HEADER.append(ENUMERATION_COLUMN_HEADERS);
    } else if (solver == GREEDY) {
      CSV_HEADER.append(",");
      CSV_HEADER.append(GREEDY_COLUMN_HEADERS);
    }
  }
  result_file << CSV_HEADER << std::endl;
  while (entity != NULL) {
    std::string name = entity->d_name;
    entity = readdir(set_directory);
    if (name.size() < set_name.size())
      continue;  // This is not a data file, data files start with <set_name> so
                 // are at least as long as set_name.size().
    std::string set_name_in_name(name.begin(), name.begin() + set_name.size());
    if (set_name_in_name != set_name)
      continue;  // For sanity, check that the data file matches the set_name
                 // which should always be the case.
    if (IsRunFile(name))
      continue;  // Cannot run on a output data file (in case we are running
                 // again on this directory). Since we include full timestamps
                 // in the result file names we may do this, e.g. with a
                 // different max policies, and won't ovewrite result file
                 // names.
    // We will loop through all cost files (which are analogous to all
    // InstanceInputs) parse their names to get params for both GraphInput and
    // InstanceInput, and then run the our single run function for every k from
    // 1 to max k (as long as k <= scenarios).
    if (IsCostFile(name)) {
      std::pair<int, int> n_kzero = GetNodesKZero(name, set_name.size() + 1);
      int nodes = n_kzero.first;
      int kzero = n_kzero.second;
      std::pair<int, int> scenarios_id = GetScenariosID(name);
      int scenarios = scenarios_id.first;
      int id = scenarios_id.second;
      GraphInput graph_input(set_name, nodes, kzero);
      InstanceInput instance_input(graph_input, scenarios, id);
      for (int k = min_policies; k <= max_policies; k++) {
        for (int budget = min_budget; budget <= max_budget; budget++) {
          if (k > instance_input.scenarios_) break;
          // removing budget > k_zero check for now, because I want to run some
          // instances of that type if (budget > graph_input.k_zero_) break;
          const ProblemInput problem(
              instance_input, k, budget, env, manual_symmetry_constraints,
              gurobi_symmetry_detection, greedy_mip_gap_threshold);
          ProblemInput problem_copyable(
              instance_input, k, budget, env, manual_symmetry_constraints,
              gurobi_symmetry_detection, greedy_mip_gap_threshold);
          std::cout << "RUNNING INSTANCE: " << problem.instance_.name()
                    << ", K = " << std::to_string(k) << std::endl;
          std::string result = SolveAndPrintTest(
              set_name, problem, problem_copyable, solvers, DEBUG);
          std::cout << std::endl;
          result_file << result << std::endl;
        }
      }
    }
  }
  result_file.close();
  closedir(set_directory);
  delete env;
  delete[] full_path;
}

std::string SolveAndPrintUninterdicted(const std::string& set_name,
                                       const ProblemInput& problem) {
  std::string final_csv_string = set_name;
  final_csv_string.append(",");
  final_csv_string.append(problem.instance_.name());
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.G_.nodes()));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.G_.arcs()));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.k_zero_));
  double nodes = problem.G_.nodes();
  double arcs = problem.G_.arcs();
  double density = (arcs) / ((nodes) * (nodes - 1));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(density));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.instance_.scenarios()));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.budget_));
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(problem.policies_));

  // MIP Fill-in. (Put uninterdicted solution here).
  final_csv_string.append(",");
  final_csv_string.append(",");
  final_csv_string.append(std::to_string(
      problem.instance_.SolveASPIZeroPolicies(problem.G_, problem.env_)));
  final_csv_string.append(",");
  final_csv_string.append(",");
  final_csv_string.append(",");
  // Benders Fill-in.
  final_csv_string.append(",");
  final_csv_string.append(",");
  final_csv_string.append(",");
  final_csv_string.append(",");
  final_csv_string.append(",");
  final_csv_string.append(",");
  // Enumeration Fill-in.
  final_csv_string.append(",");
  final_csv_string.append(",");
  final_csv_string.append(",");
  final_csv_string.append(",");
  // Greedy Fill-In.
  final_csv_string.append(",");
  final_csv_string.append(",");
  final_csv_string.append(",");
  std::cout << problem.instance_.name() << " done." << std::endl;
  return final_csv_string;
}

void UninterdictedObjectiveForAllInstances(const std::string& set_name) {
  GRBEnv* env = new GRBEnv();  // Initialize global gurobi environment.
  // Use seconds, since gurobi takes the parameter value in seconds.
  env->set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_S);  // Set time limit.
  std::cout << std::endl;
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
  std::string base_name = DATA_DIRECTORY;
  base_name.append(set_name);
  base_name.append("/");
  base_name.append(set_name);
  base_name.append("_run_");
  base_name.append(today);
  std::string result_file_name = base_name;
  result_file_name.append(".csv");
  std::ofstream result_file(result_file_name);

  std::string CSV_HEADER = INSTANCE_INFO_COLUMN_HEADERS;
  CSV_HEADER.append(",");
  CSV_HEADER.append(MIP_COLUMN_HEADERS);
  CSV_HEADER.append(",");
  CSV_HEADER.append(BENDERS_COLUMN_HEADERS);
  CSV_HEADER.append(",");
  CSV_HEADER.append(ENUMERATION_COLUMN_HEADERS);
  CSV_HEADER.append(",");
  CSV_HEADER.append(GREEDY_COLUMN_HEADERS);
  result_file << CSV_HEADER << std::endl;
  while (entity != NULL) {
    std::string name = entity->d_name;
    entity = readdir(set_directory);
    if (name.size() < set_name.size())
      continue;  // This is not a data file, data files start with <set_name> so
                 // are at least as long as set_name.size().
    std::string set_name_in_name(name.begin(), name.begin() + set_name.size());
    if (set_name_in_name != set_name)
      continue;  // For sanity, check that the data file matches the set_name
                 // which should always be the case.
    if (IsRunFile(name))
      continue;  // Cannot run on a output data file (in case we are running
                 // again on this directory). Since we include full timestamps
                 // in the result file names we may do this, e.g. with a
                 // different max policies, and won't ovewrite result file
                 // names.
    if (IsCostFile(name)) {
      std::pair<int, int> n_kzero = GetNodesKZero(name, set_name.size() + 1);
      int nodes = n_kzero.first;
      int kzero = n_kzero.second;
      std::pair<int, int> scenarios_id = GetScenariosID(name);
      int scenarios = scenarios_id.first;
      int id = scenarios_id.second;
      GraphInput graph_input(set_name, nodes, kzero);
      InstanceInput instance_input(graph_input, scenarios, id);
      const ProblemInput problem(
          instance_input, /*policies=*/0, /*budget=*/0, env,
          /*manual_symmetry_constraints=*/
          MANUAL_SYMMETRY_NONE,
          /*gurobi_symmetry_detection=*/GUROBI_SYMMETRY_AUTO,
          /*greedy_mip_gap_threshold=*/0);
      std::string result = SolveAndPrintUninterdicted(set_name, problem);
      std::cout << std::endl;
      result_file << result << std::endl;
    }
  }
  result_file.close();
  closedir(set_directory);
  delete env;
  delete[] full_path;
}
