#include "solvers.h"

long GetCurrentTime() {
  // Helper function to get current time.
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000 + tv.tv_usec / 1000;
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
  policies_ = adaptive_instance->policies_;
  budget_ = adaptive_instance->budget_;
  nodes_ = adaptive_instance->nodes_;
  arcs_ = adaptive_instance->arcs_;
  interdiction_delta_ = adaptive_instance->interdiction_delta_;
  arc_costs_ = std::vector<std::vector<int>>(scenarios_);
  scenario_index_map_ = std::vector<int>(scenarios_);
  for (int q = 0; q < scenarios_; q++) {
    arc_costs_[q] = adaptive_instance->arc_costs_[keep_scenarios[q]];
    scenario_index_map_[q] = keep_scenarios[q];
  }
}

void AdaptiveInstance::ReadCosts(int interdiction_delta) {
  // Read arc costs from a file.
  interdiction_delta_ = interdiction_delta;
  std::string line, word;
  std::string filename =
      directory_ + name_ + "-costs_" + std::to_string(scenarios_) + ".csv";
  std::ifstream myfile(filename);
  int q = 0;
  std::vector<int> v;
  int cost;
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      v.clear();
      std::stringstream str(line);
      while (getline(str, word, ',')) {
        cost = std::stoi(word);
        v.push_back(cost);
      }
      if (q < scenarios_) {
        arc_costs_.push_back(v);
      }
      ++q;
    }
  }
}

void AdaptiveInstance::PrintInstance(const Graph& G) const {
  // Print Summary of Problem Instance
  std::cout << "k: " << policies_ << std::endl;
  std::cout << "p: " << scenarios_ << std::endl;
  G.PrintGraphWithCosts(arc_costs_, interdiction_delta_);
}

int AdaptiveInstance::Dijkstra(int q, const Graph& G) {
  // Compute shortest path 0-n-1 and just return its objective (we never
  // actually need the path).
  std::vector<int> pred(nodes_), distance(nodes_);
  // First item in the pair is the weight/cost, second is the vertex.
  // std::greater allows the smallest distance to appear at top.
  std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
                      std::greater<std::pair<int, int>>>
      pq;
  std::unordered_set<int> visited;
  visited.insert(0);
  pq.push({0, 0});
  pred[0] = 0;
  distance[0] = 0;
  for (int i = 1; i < nodes_; ++i) {
    pred[i] = -1;
    distance[i] = INT_MAX;
  }
  while (!pq.empty()) {
    int node = pq.top().second;
    int node_distance = pq.top().first;
    pq.pop();
    visited.insert(node);
    for (size_t i = 0; i < G.arc_index_hash()[node].size(); ++i) {
      int arc = G.arc_index_hash()[node][i];
      int v = G.adjacency_list()[node][i];
      if (node_distance + arc_costs_[q][arc] < distance[v]) {
        pq.push(std::make_pair(node_distance + arc_costs_[q][arc], v));
        pred[v] = node;
        distance[v] = node_distance + arc_costs_[q][arc];
      }
    }
  }
  return distance[nodes_ - 1];
}

void AdaptiveInstance::ApplyInterdiction(const std::vector<double>& x_bar,
                                         bool rev) {
  // Receive interdiction policy and update costs for the M3 instance
  // If rev, we are "removing" the interdiction policy and returning the
  // instance to its original state
  for (int a = 0; a < arcs_; ++a) {
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

double AdaptiveInstance::ComputeObjectiveOfPolicyForScenario(
    const std::vector<double>& binary_policy, const Graph& G, int q) {
  ApplyInterdiction(binary_policy);
  int objective = Dijkstra(q, G);
  ApplyInterdiction(binary_policy, true);
  return objective;
}

double AdaptiveInstance::ValidatePolicy(std::vector<double>& x_bar,
                                        const Graph& G) {
  // Solve shortest path on interdicted graph - check objectives match.
  double objective = DBL_MAX;
  int sp_result;
  ApplyInterdiction(x_bar);
  for (int q = 0; q < scenarios_; ++q) {
    sp_result = Dijkstra(q, G);
    if (sp_result < objective) {
      objective = sp_result;
    }
  }
  ApplyInterdiction(x_bar, true);
  return objective;
}

void SetPartitioningModel::ConfigureSolver(const Graph& G,
                                           AdaptiveInstance& instance) {
  try {
    // ------ Initialize solution to appropriate size -----
    current_solution_ =
        AdaptiveSolution(false, instance, policies_, scenarios_);
    // ------ Initialize model and environment ------
    sp_model_ = new GRBModel(*env_);
    sp_model_->set(GRB_IntParam_OutputFlag, 0);

    // ------ Decision variables ------
    std::vector<GRBVar> new_vector;
    // set partitioning variables
    for (int w = 0; w < policies_; ++w) {
      h_var_.push_back(new_vector);
      for (int q = 0; q < scenarios_; ++q) {
        std::string varname =
            "H_" + std::to_string(w) + "_" + std::to_string(q);
        h_var_[w].push_back(sp_model_->addVar(0, 1, 0, GRB_BINARY, varname));
      }
    }
    // interdiction policy on arcs 'x' - REMEMBER to check that this
    // initialization works
    for (int w = 0; w < policies_; ++w) {
      x_var_.push_back(new_vector);

      for (int a = 0; a < arcs_; a++) {
        std::string varname =
            "x_" + std::to_string(w) + "_" + std::to_string(a);
        x_var_[w].push_back(sp_model_->addVar(0, 1, 0, GRB_BINARY, varname));
      }
    }
    // objective func dummy 'z'
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
          for (size_t j = 0; j < G.adjacency_list()[i].size(); ++j) {
            jn = G.adjacency_list()[i][j];
            a = G.arc_index_hash()[i][j];
            // add constraint
            sp_model_->addConstr(
                (pi_var_[w][q][jn] - pi_var_[w][q][i] - lambda_var_[w][q][a]) <=
                instance.arc_costs()[q][a] +
                    (instance.interdiction_delta() * x_var_[w][a]));
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
              << e.getErrorCode() << "\n";
    std::cout << e.getMessage() << "\n";
  } catch (...) {
    std::cout << "Non-gurobi error during optimization "
                 "[SetPartitioningModel::ConfigureSolver]"
              << "\n";
  }
}

AdaptiveSolution SetPartitioningModel::Solve() {
  // Updates current solution to latest solution, or sets unbounded to true.
  // Returns current solution.
  // Compute objectives (of each partition) after solving with
  // AdaptiveSolution::ComputeAllObjectives.
  long begin = GetCurrentTime();
  sp_model_->optimize();
  long time = GetCurrentTime() - begin;
  current_solution_.set_solution_time(time);
  if (sp_model_->get(GRB_IntAttr_Status) == 2) {
    current_solution_.set_unbounded(false);
    current_solution_.set_worst_case_objective(
        -sp_model_->get(GRB_DoubleAttr_ObjVal));
    for (int w = 0; w < policies_; ++w) {
      std::vector<double> x_vector;
      for (int a = 0; a < arcs_; a++) {
        x_vector.push_back(x_var_[w][a].get(GRB_DoubleAttr_X));
      }
      Policy temp_policy = Policy(-1, x_vector);
      current_solution_.set_solution_policy(w, temp_policy);
      for (int q = 0; q < scenarios_; q++) {
        if (h_var_[w][q].get(GRB_DoubleAttr_X) > 0.5) {
          current_solution_.add_to_partition(w, q);
        }
      }
    }
  } else if (sp_model_->get(GRB_IntAttr_Status) == 5)
    current_solution_.set_unbounded(true);
  return current_solution_;
}

void BendersCallback::ConfigureIndividualSubModel(const Graph& G, int w,
                                                  int q) {
  submodels_[w][q]->set(GRB_IntParam_OutputFlag, 0);
  // Add decision variables.
  for (int a = 0; a < arcs_; ++a) {
    std::string varname = "y_" + std::to_string(w) + "_" + std::to_string(q) +
                          "_" + std::to_string(a);
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

void BendersCallback::ConfigureSubModels(const Graph& G,
                                         AdaptiveInstance& instance) {
  // Initialize current_solution_.
  current_solution_ = AdaptiveSolution(true, instance, policies_, scenarios_);
  // Initialize models (k x p).
  // submodels_ = std::vector<std::vector<GRBModel*>>(policies_,
  //         std::vector<GRBModel*>(scenarios_, new GRBModel(*env_)));
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
      ConfigureIndividualSubModel(G, w, q);
    }
  }
}

void BendersCallback::UpdateSubModels(bool rev) {
  for (int w = 0; w < policies_; ++w) {
    for (int a = 0; a < arcs_; ++a) {
      if (getSolution(x_var_[w][a]) > 0.5) {
        for (int q = 0; q < scenarios_; ++q) {
          int current_cost = y_var_[w][q][a].get(GRB_DoubleAttr_Obj);
          if (rev) {
            int new_cost = current_cost - interdiction_delta_;
            y_var_[w][q][a].set(GRB_DoubleAttr_Obj, new_cost);
          } else {
            int new_cost = current_cost + interdiction_delta_;
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
    // Updating upper_bound_ here, as this is equivalent to the "end" of the
    // loop iteration in the paper algorithm (i.e. gurobi just resolved the
    // master problem, which happens at the end of an iteration in the
    // algorithm, right after we do the subproblems and add cuts, and right
    // before we update the upper bound.
    upper_bound_ = getSolution(z_var_);
    if (upper_bound_ - lower_bound_ >= epsilon_) {
      UpdateSubModels();
      SolveSubModels();
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
        for (int w = 0; w < policies_; ++w) {
          std::vector<double> binary_policy(arcs_);
          for (int a = 0; a < arcs_; ++a) {
            if (getSolution(x_var_[w][a]) > 0.5)
              binary_policy[a] = 1;
            else
              binary_policy[a] = 0;
          }
          Policy policy(0, binary_policy);
          current_solution_.set_solution_policy(w, policy);
        }
      }
      AddLazyCuts();
      UpdateSubModels(true);
    }
  }
}

void SetPartitioningBenders::ConfigureSolver(const Graph& G,
                                             AdaptiveInstance& instance) {
  try {
    // Initialize Model/Environment.
    benders_model_ = new GRBModel(env_);
    benders_model_->set(GRB_IntParam_OutputFlag, 0);
    benders_model_->getEnv().set(GRB_IntParam_LazyConstraints, 1);
    // Add Decision variables - z, h(w)(q) and x(w).
    z_var_ = benders_model_->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS, "z");
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
    callback_ = BendersCallback(big_m_, instance, z_var_, h_var_, x_var_, env_);
    callback_.ConfigureSubModels(G, instance);
    // Add first Benders Cut "Manually" to immediately have an upper bound.
    callback_.SolveSubModels();
    for (int w = 0; w < policies_; ++w) {
      for (int q = 0; q < scenarios_; ++q) {
        GRBLinExpr lhs = z_var_ - big_m_ * (1 - h_var_[w][q]);
        GRBLinExpr rhs = 0;
        for (int a = 0; a < arcs_; ++a) {
          if (callback_.y_var_[w][q][a].get(GRB_DoubleAttr_X) > 0.5) {
            rhs +=
                instance.arc_costs()[q][a] + interdiction_delta_ * x_var_[w][a];
          }
        }
        benders_model_->addConstr(lhs <= rhs);
      }
    }
  } catch (GRBException e) {
    std::cout
        << "Gurobi error number [SetPartitioningBenders::ConfigureSolver]: "
        << e.getErrorCode() << "\n";
    std::cout << e.getMessage() << "\n";
  } catch (...) {
    std::cout << "Non-gurobi error during optimization "
                 "[SetPartitioningBenders::ConfigureSolver]"
              << "\n";
  }
}

AdaptiveSolution SetPartitioningBenders::Solve() {
  benders_model_->setCallback(&callback_);
  long begin = GetCurrentTime();
  benders_model_->optimize();
  long time = GetCurrentTime() - begin;
  if (benders_model_->get(GRB_IntAttr_Status) == 2) {
    callback_.current_solution_.set_unbounded(false);
    callback_.current_solution_.set_solution_time(time);
    callback_.current_solution_.set_worst_case_objective(
        -benders_model_->get(GRB_DoubleAttr_ObjVal));
    callback_.current_solution_.set_cuts(callback_.lazy_cuts_rounds());
    for (int w = 0; w < policies_; ++w) {
      std::vector<double> x_vector;
      for (int a = 0; a < arcs_; a++) {
        x_vector.push_back(x_var_[w][a].get(GRB_DoubleAttr_X));
      }
      Policy temp_policy = Policy(-1, x_vector);
      callback_.current_solution_.set_solution_policy(w, temp_policy);
      for (int q = 0; q < scenarios_; q++) {
        if (h_var_[w][q].get(GRB_DoubleAttr_X) > 0.5) {
          callback_.current_solution_.add_to_partition(w, q);
        }
      }
    }
  } else if (benders_model_->get(GRB_IntAttr_Status) == 5) {
    callback_.current_solution_.set_unbounded(true);
    callback_.current_solution_.set_solution_time(time);
  }
  return callback_.current_solution_;
}

std::vector<int> initKappa(int p, int k) {
  // initialize partition vector based on total number in set (p) and exact
  // number of partitions required
  std::vector<int> kappa(p, 0);

  for (int q = (p - k + 1); q < p; ++q) {
    kappa[q] = (q - (p - k));
  }

  return kappa;
}

void RobustAlgoModel::configureModel(const Graph& G, AdaptiveInstance& m3) {
  // Construct baseline model with no constraints
  algo_model = new GRBModel(env_);
  algo_model->set(GRB_IntParam_OutputFlag, 0);

  // Decision Variables
  std::string varname = "z";
  z = algo_model->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS,
                         varname);  // objective function dummy variable
  for (int a = 0; a < arcs; ++a) {
    varname = "x_" + std::to_string(a);
    x.push_back(algo_model->addVar(0, 1, 0, GRB_BINARY, varname));
  }  // interdiction policy on arcs

  std::vector<GRBVar> tempvector;
  for (int q = 0; q < scenarios; ++q) {
    pi.push_back(tempvector);  // post interdiction s-i best path (for every q)
    for (int i = 0; i < nodes; ++i) {
      varname = "pi_" + std::to_string(q) + "_" + std::to_string(i);
      pi[q].push_back(algo_model->addVar(-GRB_INFINITY, GRB_INFINITY, 0,
                                         GRB_CONTINUOUS, varname));
    }
  }

  for (int q = 0; q < scenarios; ++q) {
    lambda.push_back(tempvector);  // lambda variable on arcs (for every q)
    for (int a = 0; a < arcs; ++a) {
      varname = "lambda_" + std::to_string(q) + "_" + std::to_string(a);
      lambda[q].push_back(
          algo_model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
    }
  }

  // Add budget Constraint
  GRBLinExpr linexpr = 0;
  for (int a = 0; a < arcs; a++) {
    linexpr += x[a];
  }
  algo_model->addConstr(linexpr <= budget, "budget");
  algo_model->update();

  // constraint blocking interdiction of connecting arcs
  // for (int a=0; a<arcs; a++) {
  //     if (G.arcs[a].sub == -1) {
  //         algo_model->addConstr(x[a] == 0);
  //     }
  // }

  // Populate Global Constraints
  // Arc/Dual Constraints
  for (int q = 0; q < scenarios; ++q) {
    dual_constraints.push_back(std::vector<GRBTempConstr>());

    for (int i = 0; i < nodes; ++i) {
      for (size_t j = 0; j < G.adjacency_list()[i].size(); ++j) {
        int next = G.adjacency_list()[i][j];
        int a = G.arc_index_hash()[i][j];

        GRBTempConstr constraint =
            pi[q][next] - pi[q][i] - lambda[q][a] <=
            m3.arc_costs()[q][a] + (m3.interdiction_delta() * x[a]);
        dual_constraints[q].push_back(constraint);
      }
    }
  }

  // Objective Constraints
  for (int q = 0; q < scenarios; ++q) {
    GRBLinExpr linexpr = 0;
    linexpr += (pi[q][nodes - 1] - pi[q][0]);  // b^\top pi

    for (int a = 0; a < arcs; ++a) {
      linexpr += -lambda[q][a];  // u^\top \cdot lambda
    }

    z_constraints.push_back(z <= linexpr);
  }
}

void RobustAlgoModel::update(std::vector<int>& subset) {
  // add constraints to model for subset in partition
  for (int q : subset) {
    std::string zero_name = "zero_" + std::to_string(q);
    std::string z_name = "z_" + std::to_string(q);
    algo_model->addConstr(pi[q][0] == 0, zero_name);
    algo_model->addConstr(z_constraints[q], z_name);

    int a = 0;
    // std::cout << "about to loop through dual constraints" << std::endl;
    for (GRBTempConstr constraint : dual_constraints[q]) {
      std::string dual_name =
          "dual_" + std::to_string(q) + "_" + std::to_string(a);
      algo_model->addConstr(constraint, dual_name);
      ++a;
    }
  }

  algo_model->update();
}

void RobustAlgoModel::reverse_update(std::vector<int>& subset) {
  // remove all constraints from model - i.e. all constraints associated with
  // this subset
  for (int q : subset) {
    std::string zero_name = "zero_" + std::to_string(q);
    std::string z_name = "z_" + std::to_string(q);

    GRBConstr zero_constraint = algo_model->getConstrByName(zero_name);
    GRBConstr z_constraint = algo_model->getConstrByName(z_name);
    algo_model->remove(zero_constraint);
    algo_model->remove(z_constraint);

    int a = 0;
    for (GRBTempConstr constraint : dual_constraints[q]) {
      std::string dual_name =
          "dual_" + std::to_string(q) + "_" + std::to_string(a);
      GRBConstr dual_constraint = algo_model->getConstrByName(dual_name);
      algo_model->remove(dual_constraint);
      ++a;
    }
  }

  algo_model->update();
}

Policy RobustAlgoModel::Solve() {
  // solve static model, return policy and objective value
  algo_model->optimize();
  std::vector<double> binary_policy(arcs, 0);
  double objective = -algo_model->get(GRB_DoubleAttr_ObjVal);
  for (int a = 0; a < arcs; ++a) {
    binary_policy[a] = x[a].get(GRB_DoubleAttr_X);
  }
  algo_model->reset();
  Policy solution = Policy(objective, binary_policy);
  return solution;
}

int max_int(int a, int b) {
  if (a <= b) {
    return b;
  } else {
    return a;
  }
}

void AdaptiveSolution::ComputeAllObjectives(const Graph& G,
                                            AdaptiveInstance* instance,
                                            bool compute_adaptive_objective) {
  // Populate / recompute the all_objectives vector based on interdiction
  // policies and follower costs. Worst case objective for every policy will be
  // assigned to solution[w].objective (in the Policy struct). Note: this
  // requires knowing the optimal partition, which is part of the
  // AdaptiveSolution class. If the unbounded bool is true, it will set the
  // objective to DBL_MAX. THIS FUNCTION IS NOT STABLE.
  if (unbounded_) {
    std::cout << "unbounded" << std::endl;
    worst_case_objective_ = DBL_MAX;
  } else {
    std::vector<std::vector<int>> new_partition(policies_);
    std::vector<std::vector<double>> all_objectives(
        policies_, std::vector<double>(scenarios_));
    double min_max_objective = DBL_MAX;
    for (int q = 0; q < scenarios_; q++) {
      double max_objective_this_scenario = DBL_MIN;
      int subset_assignment = -1;
      for (int w = 0; w < policies_; w++) {
        double objective = instance->ComputeObjectiveOfPolicyForScenario(
            solution_[w].binary_policy(), G, q);
        all_objectives[w][q] = objective;
        if (objective > max_objective_this_scenario) {
          subset_assignment = w;
          max_objective_this_scenario = objective;
        }
      }
      if (max_objective_this_scenario < min_max_objective)
        min_max_objective = max_objective_this_scenario;
      new_partition[subset_assignment].push_back(q);
    }
    int w = 0;
    for (const std::vector<int>& subset : new_partition) {
      double min_objective_this_policy = DBL_MAX;
      for (int q : subset) {
        if (min_objective_this_policy > all_objectives[w][q])
          min_objective_this_policy = all_objectives[w][q];
      }
      solution_[w].set_objective(min_objective_this_policy);
      w++;
    }
    partition_ = new_partition;
    if (compute_adaptive_objective) worst_case_objective_ = min_max_objective;
  }
}

void AdaptiveSolution::LogSolution(const Graph& G, std::string title,
                                   bool policy) {
  std::cout << std::endl;
  if (title == "") {
    std::cout << "Solution, time (ms): " << solution_time_ << " ms";
  } else {
    std::cout << title << ", time: " << solution_time_ << " ms";
  }
  std::cout << ", k = " << policies_ << ", p = " << scenarios_;
  std::cout << ", nodes = " << nodes_ << ", arcs = " << arcs_ << std::endl;
  std::cout << "Worst Case Objective: " << worst_case_objective_ << std::endl;
  for (size_t w = 0; w < partition_.size(); ++w) {
    if (policy) {
      std::cout << "subset: { ";
      for (int q : partition_[w]) {
        std::cout << q << " ";
      }
      std::cout << "} - objective: " << solution_[w].objective() << std::endl;
      std::cout << "interdicted arcs: ";
      std::vector<std::vector<int>> arc_index_hash = G.arc_index_hash();
      std::vector<std::vector<int>> adjacency_list = G.adjacency_list();
      for (int i = 0; i < G.nodes(); ++i) {
        int index = 0;
        for (int a : arc_index_hash[i]) {
          if (solution_[w].binary_policy()[a] > 0.5) {
            G.PrintArc(a, i, index);
            std::cout << " ";
          }
          index++;
        }
      }
      std::cout << std::endl;
    }
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
        // std::cout << "second half pass loop" << std::endl;
        // std::cout << "u: " << u << std::endl;
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

std::pair<std::vector<std::vector<int>>, std::vector<Policy>> MergeEnumSols(
    std::pair<std::vector<std::vector<int>>, std::vector<Policy>>& sol1,
    std::pair<std::vector<std::vector<int>>, std::vector<Policy>>& sol2,
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

AdaptiveSolution enumSolve(AdaptiveInstance& m3, const Graph& G, GRBEnv* env) {
  // Pass a AdaptiveInstance and solve using enumeration algorithm
  // Enumeration maintained as in Orlov paper
  // OLD Return: pair:
  //      first: vector of vector of ints - optimal partition (optimal objective
  //      is with policies) second: vector of Policy - k interdiction policies
  //  NEW Return: Adaptive Solution Struct (basic info, nested vector of
  //  partitions, vector of policies)
  // ints and bool
  bool next = true;
  int p = m3.scenarios();
  int k = m3.policies();
  int m = m3.arcs();
  // initialize partitioning 'string' vector as in Orlov Paper
  // initialize corresponding max vector
  std::vector<int> kappa = initKappa(p, k);
  std::vector<int> max = kappa;

  // initialize solution vector
  std::vector<double> arc_vec(m, 0);
  std::vector<std::vector<double>> sol(k, arc_vec);

  // initialize static robust model
  try {
    long begin = GetCurrentTime();
    RobustAlgoModel static_robust = RobustAlgoModel(m3, env);
    static_robust.configureModel(G, m3);

    // objective value maintained here (and optimal partition)
    double best_worstcase_objective = 0;
    std::vector<double> final_objectives(k, 0);
    std::vector<std::vector<int>> best_worstcase_partition;

    // enumerate while not 'failing' to get next partition
    while (next) {
      std::vector<std::vector<double>> temp_sol(k, arc_vec);
      double temp_worst_objective = DBL_MAX;
      std::vector<double> temp_objectives(k, 0);
      std::vector<std::vector<int>> partition = kappa_to_partition(kappa, k, p);

      for (int w = 0; w < k; ++w) {
        // for every subset in partition, solve M2 for k=1
        static_robust.update(partition[w]);
        Policy temp_single_solution = static_robust.Solve();
        temp_objectives[w] = temp_single_solution.objective();
        static_robust.reverse_update(partition[w]);

        // update temp_worst_objective and temp_sol if this subset is worse
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

    std::vector<Policy> final_policies(k, Policy(m));

    for (int w = 0; w < k; ++w) {
      final_policies[w].set_policy(sol[w]);
      final_policies[w].set_objective(final_objectives[w]);
    }

    // std::vector<int> final_obj(1, best_worstcase_solution);
    // best_worstcase_partition.insert(best_worstcase_partition.begin(),
    // final_obj); OLD auto final_solution =
    // std::make_pair(best_worstcase_partition, final_policies); NEW
    AdaptiveSolution final_solution = AdaptiveSolution(
        false, k, p, m3, best_worstcase_partition, final_policies);
    final_solution.set_worst_case_objective(best_worstcase_objective);
    long time = GetCurrentTime() - begin;
    final_solution.set_solution_time(time);
    return final_solution;
  } catch (GRBException e) {
    std::cout << "Gurobi error number [EnumSolve]: " << e.getErrorCode()
              << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch (...) {
    std::cout
        << "Non Gurobi error during construction of static robust model object"
        << std::endl;
  }

  AdaptiveSolution dummy_solution;
  return dummy_solution;
}

void AdaptiveSolution::ExtendByOne(AdaptiveInstance& m3, const Graph& G,
                                   GRBEnv* env, bool mip_subroutine) {
  // Use an optimal solution found for k, to find a good solution for k+1
  // Take the worst subset in the optimal partition and "split it in 2"
  // "Split it in two": solve that subset for k = 2
  // Return format will always be exactly the same as enum solve, just for k
  // greater than one Naming - m3_prime is the copy of m3, and anything_prime is
  // associated with m3_prime int values
  std::cout << "heuristic - extend by one policy" << std::endl;
  long begin = GetCurrentTime();
  // find worst subset in optimal partition
  double min_subset_obj = GRB_INFINITY;
  int min_subset_windex;

  for (int w = 0; w < policies_; ++w) {
    if (solution_[w].objective() < min_subset_obj) {
      min_subset_obj = solution_[w].objective();
      min_subset_windex = w;
    }
  }

  int p_prime = partition_[min_subset_windex]
                    .size();  // the new p - number of scenarios in the subset
                              // that we will work on

  if (p_prime > 1) {
    std::cout << "Subset to split: { ";
    for (int q : partition_[min_subset_windex]) {
      std::cout << q << " ";
    }
    std::cout << "}" << std::endl;
    AdaptiveInstance m3_prime =
        AdaptiveInstance(&m3, partition_[min_subset_windex]);
    // careful - m3_prime always has k=2 because we are just extending by ONE
    m3_prime.set_policies(2);
    // std::cout << std::endl << std::endl;
    // std::cout << "m3 prime: " << std::endl;
    // m3_prime.printInstance(G);

    AdaptiveSolution k_prime_solution;

    if (mip_subroutine) {
      SetPartitioningModel m3prime_model =
          SetPartitioningModel(500, m3_prime, env);
      m3prime_model.ConfigureSolver(G, m3_prime);
      m3prime_model.Solve();
      k_prime_solution = m3prime_model.current_solution();
      k_prime_solution.ComputeAllObjectives(G, &m3_prime);
    } else {
      k_prime_solution = enumSolve(m3_prime, G, env);
    }

    // auto final_solution = MergeEnumSols(k_solution, k_prime_solution,
    // min_subset_windex);
    long time = GetCurrentTime() - begin;

    MergeEnumSols(k_prime_solution, &m3_prime, min_subset_windex);
    solution_time_ = time;
  }
}

std::vector<Policy> InitializeKPolicies(AdaptiveInstance* m3, const Graph& G,
                                        GRBEnv* env) {
  std::unordered_set<int> centers;
  RobustAlgoModel spi_model = RobustAlgoModel(*m3, env);
  std::vector<Policy> initial_policies;
  spi_model.configureModel(G, *m3);
  // Choose the first scenario to solve SPI for arbitrarily (we'll just use
  // index 0).
  centers.insert(0);
  std::vector<int> update_vector(1);
  spi_model.update(update_vector);
  Policy single_policy = spi_model.Solve();
  spi_model.reverse_update(update_vector);
  initial_policies.push_back(single_policy);
  size_t policies = m3->policies();
  while (centers.size() < policies) {
    // Evaluate all non center scenarios against the current policies,
    // maintaining the minimum one.
    std::pair<int, double> min_subset = {-1, DBL_MAX};
    for (int q = 0; q < m3->scenarios(); q++) {
      if (centers.count(q) == 1) continue;
      for (const Policy& policy : initial_policies) {
        double objective = m3->ComputeObjectiveOfPolicyForScenario(
            policy.binary_policy(), G, q);
        if (objective < min_subset.second) min_subset = {q, objective};
      }
    }
    // Solve SPI for scenario that is currently performing worst (for the
    // interdictor), and add the scenario to the set of centers.
    update_vector[0] = min_subset.first;
    spi_model.update(update_vector);
    Policy single_policy = spi_model.Solve();
    spi_model.reverse_update(update_vector);
    initial_policies.push_back(single_policy);
    centers.insert(min_subset.first);
  }
  return initial_policies;
}

double UpdateCurrentObjectiveGivenSolution(AdaptiveSolution* current_solution,
                                           AdaptiveInstance* m3,
                                           const Graph& G) {
  int k = m3->policies();
  int p = m3->scenarios();
  std::vector<std::vector<int>> partition(k);
  double min_max_objective = DBL_MAX;
  for (int q = 0; q < p; q++) {
    double max_objective_this_scenario = DBL_MIN;
    int subset_assignment = -1;
    for (int w = 0; w < k; w++) {
      double objective = m3->ComputeObjectiveOfPolicyForScenario(
          current_solution->solution()[w].binary_policy(), G, q);
      // std::cout << "q: " << q << ", w: " << w << ", obj: "<< objective <<
      // std::endl;
      if (objective > max_objective_this_scenario) {
        subset_assignment = w;
        max_objective_this_scenario = objective;
      }
    }
    if (max_objective_this_scenario < min_max_objective)
      min_max_objective = max_objective_this_scenario;
    partition[subset_assignment].push_back(q);
  }
  current_solution->set_partition(partition);

  return min_max_objective;
}

AdaptiveSolution KMeansHeuristic(AdaptiveInstance* m3, const Graph& G,
                                 GRBEnv* env) {
  // Use the KMeans - Style heuristic to solve an Adaptive Instance.
  long begin = GetCurrentTime();
  double z_previous = DBL_MIN;
  // Initialize the first solution - solving the non-robust model for k
  // followers chosen at random out of p.
  std::vector<Policy> first_solution = InitializeKPolicies(m3, G, env);
  std::vector<std::vector<int>> partitions(m3->policies());
  AdaptiveSolution current_solution = AdaptiveSolution(
      false, m3->policies(), m3->scenarios(), *m3, partitions, first_solution);
  // Update the objective (and initialize the current objective).
  current_solution.ComputeAllObjectives(G, m3);
  double z_current = current_solution.worst_case_objective();
  RobustAlgoModel static_model = RobustAlgoModel(*m3, env);
  static_model.configureModel(G, *m3);
  int counter = 0;
  while (z_previous < z_current) {
    std::cout << "objective, iteration " << counter << ": " << z_current
              << std::endl;
    // ONLY UPDATE MIN OBJECTIVE POLICY
    double min_objective = DBL_MAX;
    int subset = -1;
    for (int w = 0; w < m3->policies(); w++) {
      if (current_solution.solution()[w].objective() < min_objective) {
        min_objective = current_solution.solution()[w].objective();
        subset = w;
      }
    }
    static_model.update(current_solution.partition()[subset]);
    Policy new_policy = static_model.Solve();
    current_solution.set_solution_policy(subset, new_policy);
    static_model.reverse_update(current_solution.partition()[subset]);
    // UPDATE All POLICIES
    // for (int subset=0; subset<m3->policies(); subset++) {
    //     static_model.update(current_solution.partition()[subset]);
    //     Policy new_policy = static_model.Solve();
    //     current_solution.set_solution_policy(subset, new_policy);
    //     static_model.reverse_update(current_solution.partition()[subset]);
    // }
    // REASSIGNMENT
    current_solution.ComputeAllObjectives(G, m3);
    z_previous = z_current;
    z_current = current_solution.worst_case_objective();
    counter++;
  }
  std::cout << "objective, iteration " << counter << ": " << z_current
            << std::endl;
  long runtime = GetCurrentTime() - begin;
  current_solution.set_solution_time(runtime);
  return current_solution;
}
