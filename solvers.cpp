#include "solvers.h"

long GetCurrentTime() {
  // Helper function to get current time in milliseconds.
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

void PrintSeparator() {
  std::cout << "-------------------------------" << std::endl;
}

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

AdaptiveInstance::AdaptiveInstance(AdaptiveInstance* adaptive_instance,
                                   std::vector<int>& keep_scenarios) {
  // Copy constructor only keeping a subset of the scenarios
  scenarios_ = keep_scenarios.size();
  nodes_ = adaptive_instance->nodes_;
  interdiction_delta_ = adaptive_instance->interdiction_delta_;
  arc_costs_ = std::vector<std::vector<int>>(scenarios_);
  scenario_index_map_ = std::vector<int>(scenarios_);
  for (int q = 0; q < scenarios_; q++) {
    arc_costs_[q] = adaptive_instance->arc_costs_[keep_scenarios[q]];
    scenario_index_map_[q] = keep_scenarios[q];
  }
}

void AdaptiveInstance::ReadCosts() {
  // Read arc costs from a file.
  // Compute interdiction_delta_, big_m_:
  // interdiction_delta_ = (max_q,a c_a * (n)).
  // big_m_ = (max_q,a c_a * (n)) + 1.
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

double AdaptiveInstance::SPModel(int q, const Graph& G, GRBEnv* env) {
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
    return sp_model->get(GRB_DoubleAttr_ObjVal);
  } catch (GRBException e) {
    std::cout
        << "Gurobi error number [SetPartitioningBenders::ConfigureSolver]: "
        << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  } catch (...) {
    std::cout << "Non-gurobi error during optimization "
                 "[SetPartitioningBenders::ConfigureSolver]"
              << std::endl;
    return -1;
  }
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

void SetPartitioningModel::ConfigureSolver(const ProblemInput& problem) {
  try {
    sp_model_ = new GRBModel(*env_);
    sp_model_->set(GRB_IntParam_OutputFlag, 0);
    // Set partitioning variables.
    std::vector<GRBVar> new_vector;
    for (int w = 0; w < policies_; ++w) {
      h_var_.push_back(new_vector);
      for (int q = 0; q < scenarios_; ++q) {
        std::string varname =
            "H_" + std::to_string(w) + "_" + std::to_string(q);
        h_var_[w].push_back(sp_model_->addVar(0, 1, 0, GRB_BINARY, varname));
      }
    }
    // Interdiction policy on arcs.
    for (int w = 0; w < policies_; ++w) {
      x_var_.push_back(new_vector);
      for (int a = 0; a < arcs_; a++) {
        std::string varname =
            "x_" + std::to_string(w) + "_" + std::to_string(a);
        x_var_[w].push_back(sp_model_->addVar(0, 1, 0, GRB_BINARY, varname));
      }
    }
    // Objective value variable, piecewise linearization.
    z_var_ = sp_model_->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS);
    std::vector<std::vector<GRBVar>> newnew_vector;
    // post interdiction flow 'pi'
    for (int w = 0; w < policies_; ++w) {
      pi_var_.push_back(newnew_vector);
      for (int q = 0; q < scenarios_; q++) {
        pi_var_[w].push_back(new_vector);
        for (int i = 0; i < nodes_; i++) {
          std::string varname = "pi_" + std::to_string(w) + "_" +
                                std::to_string(q) + "_" + std::to_string(i);
          pi_var_[w][q].push_back(sp_model_->addVar(
              -GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
        }
      }
    }
    // arc variable 'lambda'
    for (int w = 0; w < policies_; ++w) {
      lambda_var_.push_back(newnew_vector);
      for (int q = 0; q < scenarios_; q++) {
        lambda_var_[w].push_back(new_vector);
        for (int a = 0; a < arcs_; a++) {
          std::string varname = "lambda_" + std::to_string(w) + "_" +
                                std::to_string(q) + "_" + std::to_string(a);
          lambda_var_[w][q].push_back(
              sp_model_->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
        }
      }
    }
    // ------ Constraints ------
    // budget constraint
    for (int w = 0; w < policies_; ++w) {
      GRBLinExpr linexpr = 0;
      for (int a = 0; a < arcs_; a++) {
        linexpr += x_var_[w][a];
      }
      sp_model_->addConstr(linexpr <= budget_);
    }
    // z constraints
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
    // main constraint for each arc
    int i;
    int jn;  // node j
    int a;
    for (int w = 0; w < policies_; ++w) {
      for (int q = 0; q < scenarios_; ++q) {
        for (i = 0; i < nodes_; ++i) {
          for (size_t j = 0; j < problem.G_.adjacency_list()[i].size(); ++j) {
            jn = problem.G_.adjacency_list()[i][j];
            a = problem.G_.arc_index_hash()[i][j];
            // add constraint
            sp_model_->addConstr(
                (pi_var_[w][q][jn] - pi_var_[w][q][i] - lambda_var_[w][q][a]) <=
                problem.instance_.arc_costs()[q][a] +
                    (problem.instance_.interdiction_delta() * x_var_[w][a]));
          }
        }
      }
    }
    // set-partitioning constraint
    for (int q = 0; q < scenarios_; ++q) {
      GRBLinExpr linexpr = 0;
      for (int w = 0; w < policies_; ++w) {
        linexpr += h_var_[w][q];
      }
      sp_model_->addConstr(linexpr == 1);
    }
    // pi[0] = 0
    for (int w = 0; w < policies_; ++w) {
      for (int q = 0; q < scenarios_; q++) {
        sp_model_->addConstr(pi_var_[w][q][0] == 0);
      }
    }
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
  }
  else if (sp_model_->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
    // If solution status is hit time limit, set time to TIME_LIMIT_MS.
    // Use milliseconds, as our adaptive solution class stores these values in milliseconds,
    // since GetCurrentTime() returns milliseconds.
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

void SetPartitioningBenders::ConfigureSolver(const ProblemInput& problem) {
  try {
    // Initialize Model/Environment.
    benders_model_ = new GRBModel(env_);
    benders_model_->set(GRB_IntParam_OutputFlag, 0);
    benders_model_->getEnv().set(GRB_IntParam_LazyConstraints, 1);
    // Add Decision variables - z, h(w)(q) and x(w).
    z_var_ = benders_model_->addVar(0, big_m_, -1, GRB_CONTINUOUS, "z");
    h_var_ = std::vector<std::vector<GRBVar>>(policies_,
                                              std::vector<GRBVar>(scenarios_));
    for (int w = 0; w < policies_; ++w) {
      for (int q = 0; q < scenarios_; ++q) {
        std::string varname =
            "H_" + std::to_string(w) + "_" + std::to_string(q);
        h_var_[w][q] = benders_model_->addVar(0, 1, 0, GRB_BINARY, varname);
      }
    }
    x_var_ =
        std::vector<std::vector<GRBVar>>(policies_, std::vector<GRBVar>(arcs_));
    for (int w = 0; w < policies_; ++w) {
      for (int a = 0; a < arcs_; a++) {
        std::string varname =
            "x_" + std::to_string(w) + "_" + std::to_string(a);
        x_var_[w][a] = benders_model_->addVar(0, 1, 0, GRB_BINARY, varname);
      }
    }
    // Add Budget Constraints.
    for (int w = 0; w < policies_; ++w) {
      GRBLinExpr linexpr = 0;
      for (int a = 0; a < arcs_; a++) {
        linexpr += x_var_[w][a];
      }
      benders_model_->addConstr(linexpr <= budget_);
    }
    // Add Set Partitioning Constraint.
    for (int q = 0; q < scenarios_; ++q) {
      GRBLinExpr linexpr = 0;
      for (int w = 0; w < policies_; ++w) {
        linexpr += h_var_[w][q];
      }
      benders_model_->addConstr(linexpr == 1);
    }
    // Initialize Separation Object, passing decision variables.
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
  }
  else if (benders_model_->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
    // If we hit time limit, check for incumbent and set policy and objective.
    // If solution status is hit time limit, set time to TIME_LIMIT_MS.
    // Use milliseconds, as our adaptive solution class stores these values in milliseconds,
    // since GetCurrentTime() returns milliseconds.
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
  }
  else if (algo_model_->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
    std::vector<double> binary_policy(arcs_, 0);
    double objective = -algo_model_->get(GRB_DoubleAttr_ObjVal);
    for (int a = 0; a < arcs_; ++a) {
      binary_policy[a] = x_var_[a].get(GRB_DoubleAttr_X);
    }
    return std::make_pair(objective, binary_policy);
  }
  else return std::make_pair(0, std::vector<double>(arcs_, 0));
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
  // available in objectives_). The uninitialized partition for the
  // approximation algorithm solution is a matrix of 0s, so we set
  // partition_[w][q] to 1, once we have identified the correct w.
  std::cout << "AdaptiveSolution::ComputePartition" << std::endl;
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
    partition_[partition_index][q] = 1;
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
  }
  if (benders_)
    std::cout << std::endl
              << "Rounds of lazy cuts: " << lazy_cuts_rounds_ << std::endl;
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

AdaptiveSolution EnumSolve(const ProblemInput& problem) {
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
  // Initialize static robust model.
  try {
    long begin = GetCurrentTime();
    RobustAlgoModel static_robust = RobustAlgoModel(problem);
    static_robust.ConfigureModel(problem);
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
        AdaptiveSolution time_limit_solution = AdaptiveSolution(false, false, false);
        time_limit_solution.set_solution_time(TIME_LIMIT_MS);
        time_limit_solution.set_worst_case_objective(-1);
        return time_limit_solution;
      } 
      std::vector<std::vector<double>> temp_sol(k, arc_vec);
      double temp_worst_objective = DBL_MAX;
      std::vector<double> temp_objectives(k, 0);
      std::vector<std::vector<int>> partition = kappa_to_partition(kappa, k, p);
      for (int w = 0; w < k; ++w) {
        // For every subset in partition, solve M2 for k=1.
        static_robust.Update(partition[w]);
        std::pair<double, std::vector<double>> temp_single_solution =
            static_robust.Solve();
        // Check that the static model didn't hit the time limit:
        if (temp_single_solution.first == -1) {
          AdaptiveSolution time_limit_solution = AdaptiveSolution(false, false, false);
          time_limit_solution.set_solution_time(TIME_LIMIT_MS);
          time_limit_solution.set_worst_case_objective(-1);
          return time_limit_solution;
        }
        temp_objectives[w] = temp_single_solution.first;
        temp_sol[w] = temp_single_solution.second;
        static_robust.ReverseUpdate(partition[w]);
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
      AdaptiveSolution time_limit_solution = AdaptiveSolution(false, false, false);
      time_limit_solution.set_solution_time(TIME_LIMIT_MS);
      time_limit_solution.set_worst_case_objective(-1);
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

// void AdaptiveSolution::ExtendByOne(AdaptiveInstance& instance, const Graph&
// G,
//                                    GRBEnv* env, bool mip_subroutine) {
//   // Use an optimal solution found for k, to find a good solution for k+1
//   // Take the worst subset in the optimal partition and "split it in 2"
//   // "Split it in two": solve that subset for k = 2
//   // Return format will always be exactly the same as enum solve, just for k
//   // greater than one Naming - instance_prime is the copy of instance, and
//   // anything_prime is associated with instance_prime int values
//   std::cout << "heuristic - extend by one policy" << std::endl;
//   long begin = GetCurrentTime();
//   // find worst subset in optimal partition
//   double min_subset_obj = DBL_MAX;
//   int min_subset_windex;
//
//   for (int w = 0; w < policies_; ++w) {
//     if (solution_[w].objective() < min_subset_obj) {
//       min_subset_obj = solution_[w].objective();
//       min_subset_windex = w;
//     }
//   }
//
//   int p_prime = partition_[min_subset_windex]
//                     .size();  // the new p - number of scenarios in the
//                     subset
//                               // that we will work on
//
//   if (p_prime > 1) {
//     std::cout << "Subset to split: { ";
//     for (int q : partition_[min_subset_windex]) {
//       std::cout << q << " ";
//     }
//     std::cout << "}" << std::endl;
//     AdaptiveInstance instance_prime =
//         AdaptiveInstance(&instance, partition_[min_subset_windex]);
//     // careful - instance_prime always has k=2 because we are just extending
//     by
//     // ONE
//     instance_prime.set_policies(2);
//     // std::cout << std::endl << std::endl;
//     // std::cout << "instance prime: " << std::endl;
//     // instance_prime.printInstance(G);
//
//     AdaptiveSolution k_prime_solution;
//     if (mip_subroutine) {
//       SetPartitioningModel instanceprime_model =
//           SetPartitioningModel(500, instance_prime, env);
//       instanceprime_model.ConfigureSolver(G, instance_prime);
//       instanceprime_model.Solve();
//       k_prime_solution = instanceprime_model.current_solution();
//       k_prime_solution.ComputeAllObjectives(G, instance_prime);
//     } else {
//       k_prime_solution = EnumSolve(instance_prime, G, env);
//     }
//     long time = GetCurrentTime() - begin;
//     MergeEnumSols(k_prime_solution, &instance_prime, min_subset_windex);
//     solution_time_ = time;
//   }
// }

AdaptiveSolution GreedyAlgorithm(const ProblemInput& problem) {
  long begin = GetCurrentTime();
  RobustAlgoModel static_robust = RobustAlgoModel(problem);
  std::unordered_set<int> centers;
  RobustAlgoModel spi_model = RobustAlgoModel(problem);
  std::vector<std::vector<double>> final_policies;
  spi_model.ConfigureModel(problem);
  // Choose the first scenario to solve SPI for arbitrarily (we'll just use
  // index 0).
  centers.insert(0);
  std::vector<int> Update_vector(1);
  spi_model.Update(Update_vector);
  std::vector<double> single_policy = spi_model.Solve().second;
  spi_model.ReverseUpdate(Update_vector);
  final_policies.push_back(single_policy);
  size_t policies = problem.policies_;
  while (centers.size() < policies) {
    long current_time = GetCurrentTime() - begin;
    // Use milliseconds since we used GetCurrentTime, which is in milliseconds, to measure.
    if (current_time >= TIME_LIMIT_MS) {
      AdaptiveSolution time_limit_solution = AdaptiveSolution(false, false, false);
      time_limit_solution.set_solution_time(TIME_LIMIT_MS);
      time_limit_solution.set_worst_case_objective(-1);
      return time_limit_solution;
    } 
    // Evaluate all non center scenarios against the current policies,
    // maintaining the minimum one.
    std::pair<int, double> min_subset = {-1, DBL_MAX};
    for (int q = 0; q < problem.instance_.scenarios(); q++) {
      if (centers.count(q) == 1) continue;
      for (const std::vector<double>& policy : final_policies) {
        double objective =
            ComputeObjectiveOfPolicyForScenario(problem, policy, q);
        if (objective < min_subset.second) min_subset = {q, objective};
      }
    }
    // Solve SPI for scenario that is currently performing worst (for the
    // interdictor), and add the scenario to the set of centers.
    Update_vector[0] = min_subset.first;
    spi_model.Update(Update_vector);
    std::pair<double, std::vector<double>> temp_single_solution = spi_model.Solve();
    if (temp_single_solution.first == -1) {
      AdaptiveSolution time_limit_solution = AdaptiveSolution(false, false, false);
      time_limit_solution.set_solution_time(TIME_LIMIT_MS);
      time_limit_solution.set_worst_case_objective(-1);
      return time_limit_solution;
    }
    std::vector<double> single_policy = temp_single_solution.second;
    spi_model.ReverseUpdate(Update_vector);
    final_policies.push_back(single_policy);
    centers.insert(min_subset.first);
  }
  std::vector<std::vector<int>> temp_partition(
      problem.policies_, std::vector<int>(problem.instance_.scenarios(), 0));
  // Benders == false, optimal == false.
  AdaptiveSolution final_solution(false, false, problem, temp_partition,
                                  final_policies);
  long solution_time = GetCurrentTime() - begin;
  if (solution_time >= TIME_LIMIT_MS) {
    AdaptiveSolution time_limit_solution = AdaptiveSolution(false, false, false);
    time_limit_solution.set_solution_time(TIME_LIMIT_MS);
    time_limit_solution.set_worst_case_objective(-1);
    return time_limit_solution;
  } 
  final_solution.ComputeObjectiveMatrix(problem);
  final_solution.ComputeAdaptiveObjective();
  final_solution.set_solution_time(solution_time);
  return final_solution;
}

// AdaptiveSolution KMeansHeuristic(AdaptiveInstance& instance, const Graph& G,
//                                  GRBEnv* env) {
//   // Use the KMeans - Style heuristic to solve an Adaptive Instance.
//   long begin = GetCurrentTime();
//   double z_previous = DBL_MIN;
//   // Initialize the first solution - solving the non-robust model for k
//   // followers chosen at random out of p.
//   AdaptiveSolution first_solution = GreedyAlgorithm(instance, G, env);
//   std::vector<std::vector<int>> partitions(instance.policies());
//   // AdaptiveSolution current_solution =
//   //     AdaptiveSolution(false, instance.policies(), instance.scenarios(),
//   //                      instance, partitions, first_solution);
//   // Update the objective (and initialize the current objective).
//   first_solution.ComputeAllObjectives(G, instance);
//   double z_current = first_solution.worst_case_objective();
//   RobustAlgoModel static_model = RobustAlgoModel(instance, env);
//   static_model.ConfigureModel(G, instance);
//   int counter = 0;
//   while (z_previous < z_current) {
//     std::cout << "objective, iteration " << counter << ": " << z_current
//               << std::endl;
//     // ONLY UPDATE MIN OBJECTIVE POLICY
//     double min_objective = DBL_MAX;
//     int subset = -1;
//     for (int w = 0; w < instance.policies(); w++) {
//       if (current_solution.solution()[w].objective() < min_objective) {
//         min_objective = current_solution.solution()[w].objective();
//         subset = w;
//       }
//     }
//     static_model.Update(current_solution.partition()[subset]);
//     Policy new_policy = static_model.Solve();
//     current_solution.set_solution_policy(subset, new_policy);
//     static_model.ReverseUpdate(current_solution.partition()[subset]);
//     // UPDATE All POLICIES
//     // for (int subset=0; subset<instance->policies(); subset++) {
//     //     static_model.Update(current_solution.partition()[subset]);
//     //     Policy new_policy = static_model.Solve();
//     //     current_solution.set_solution_policy(subset, new_policy);
//     //     static_model.ReverseUpdate(current_solution.partition()[subset]);
//     // }
//     // REASSIGNMENT
//     current_solution.ComputeAllObjectives(G, instance);
//     z_previous = z_current;
//     z_current = current_solution.worst_case_objective();
//     counter++;
//   }
//   std::cout << "objective, iteration " << counter << ": " << z_current
//             << std::endl;
//   long runtime = GetCurrentTime() - begin;
//   current_solution.set_solution_time(runtime);
//   return current_solution;
// }
