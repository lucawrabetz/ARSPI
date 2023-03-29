#ifndef SOLVERS_H
#define SOLVERS_H

#include <float.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <fstream>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
// #include "/home/luw28/gurobi950/linux64/include/gurobi_c++.h"
#include "/home/luchino/gurobi1001/linux64/include/gurobi_c++.h"
// #include "/Library/gurobi902/mac64/include/gurobi_c++.h"

enum ASPI_Solver { MIP, BENDERS, ENUMERATION, GREEDY };
const double EPSILON = 0.000000001;

class Graph {
  // Class to read an arc list from a file and represent it as linked lists.
 public:
  Graph(const std::string& filename, int nodes);
  void PrintArc(int a, int i, int index, bool rev = false) const;
  void PrintGraph(bool rev = false) const;
  void PrintGraphWithCosts(const std::vector<std::vector<int>>& costs,
                           int interdiction_delta) const;
  // Getters.
  int nodes() const { return nodes_; }
  int arcs() const { return arcs_; }
  std::vector<std::vector<int>> arc_index_hash() const {
    return arc_index_hash_;
  }
  std::vector<std::vector<int>> adjacency_list() const {
    return adjacency_list_;
  }
  std::vector<std::vector<int>> rev_arc_index_hash() const {
    return rev_arc_index_hash_;
  }
  std::vector<std::vector<int>> rev_adjacency_list() const {
    return rev_adjacency_list_;
  }

 private:
  int nodes_, arcs_;
  const std::string filename_;
  // Arc vectors - adjacency_list_ is a linked list representation of the arc
  // list (outgoing arcs). The vector arc_index_hash_ directly maps every arc j
  // = adjacency_list_[i] to its index a in 0,...,m-1.
  std::vector<std::vector<int>> arc_index_hash_;
  std::vector<std::vector<int>> adjacency_list_;
  // Similarly, rev_adjacency_list and rev_arc_index_hash represent a linked
  // list representation of incoming arcs.
  std::vector<std::vector<int>> rev_arc_index_hash_;
  std::vector<std::vector<int>> rev_adjacency_list_;
};

class AdaptiveInstance {
  // Full adaptive problem instance with cost data, (Graph G is always passed by
  // reference separately).
 public:
  AdaptiveInstance(int scenarios, int policies, int budget, const Graph& G,
                   const std::string& directory, const std::string& name)
      : nodes_(G.nodes()),
        arcs_(G.arcs()),
        scenarios_(scenarios),
        policies_(policies),
        budget_(budget),
        directory_(directory),
        name_(name){};
  // Copy constructor with a different U (only a subset of scenarios to keep).
  AdaptiveInstance(AdaptiveInstance* adaptive_instance,
                   std::vector<int>& keep_scenarios);
  void ReadCosts(int interdiction_delta);
  void PrintInstance(const Graph& G) const;
  int Dijkstra(int q, const Graph& G);
  void ApplyInterdiction(const std::vector<double>& x_bar, bool rev = false);
  double ValidatePolicy(std::vector<double>& x_bar, const Graph& G);
  double ComputeObjectiveOfPolicyForScenario(
      const std::vector<double>& binary_policy, const Graph& G, int q);
  int nodes() const { return nodes_; }
  int arcs() const { return arcs_; }
  int scenarios() const { return scenarios_; }
  int policies() const { return policies_; }
  int budget() const { return budget_; }
  int interdiction_delta() const { return interdiction_delta_; }
  std::string name() const { return name_; }
  std::vector<std::vector<int>> arc_costs() const { return arc_costs_; }
  std::vector<int> scenario_index_map() const { return scenario_index_map_; }
  // No setter for scenarios_ - functionality reserved for change copy
  // constructor.
  void set_policies(int policies) { policies_ = policies; }

 private:
  int nodes_, arcs_, scenarios_, policies_, budget_;
  int interdiction_delta_;
  const std::string directory_;
  const std::string name_;
  std::vector<std::vector<int>> arc_costs_;
  // Map of scenario indices (indexed 0-(p-1)) to their original indices, for
  // when an AdaptiveInstance is copied but only keeps a subset of U (to be able
  // to recover the original). An empty map indicates an original instance.
  std::vector<int> scenario_index_map_;
};

class Policy {
  // Policy class - represent a single interdiction policy with its objective.
  // The objective - is the worst (max) objective for the subsets that the
  // interdiction policy is assigned to.
 public:
  Policy() : objective_(0), binary_policy_(std::vector<double>()){};
  Policy(int size) : objective_(0), binary_policy_(std::vector<double>(size)){};
  Policy(double objective, const std::vector<double>& binary_policy)
      : objective_(objective), binary_policy_(binary_policy){};
  double objective() const { return objective_; }
  std::vector<double> binary_policy() const { return binary_policy_; }
  void set_policy(const std::vector<double>& binary_policy) {
    binary_policy_ = binary_policy;
  }
  void set_objective(double objective) { objective_ = objective; }

 private:
  double objective_;
  std::vector<double> binary_policy_;
};

class AdaptiveSolution {
  // Full solution - i.e. a vector of k policy objects with extra info,
  // including the partition.
  // TODO: are all the constructors needed?
  // TODO: default initializations of stuff for different solvers (approximation
  // and enumeration algorithm).
 public:
  AdaptiveSolution() : unbounded_(false), benders_(false){};
  AdaptiveSolution(bool unbounded, bool benders)
      : unbounded_(unbounded), benders_(benders){};
  AdaptiveSolution(bool benders, AdaptiveInstance& instance, int policies,
                   int scenarios)
      : unbounded_(false),
        benders_(benders),
        policies_(policies),
        scenarios_(scenarios),
        nodes_(instance.nodes()),
        arcs_(instance.arcs()),
        partition_(std::vector<std::vector<int>>(policies)),
        solution_(std::vector<Policy>(policies)){};
  AdaptiveSolution(bool benders, int policies, int scenarios,
                   AdaptiveInstance& instance,
                   const std::vector<std::vector<int>>& partition,
                   const std::vector<Policy>& solution)
      : unbounded_(false),
        benders_(benders),
        policies_(policies),
        scenarios_(scenarios),
        nodes_(instance.nodes()),
        arcs_(instance.arcs()),
        partition_(partition),
        solution_(solution){};
  void LogSolution(const Graph& G, std::string title, bool policy = false);
  void MergeEnumSols(AdaptiveSolution sol2, AdaptiveInstance* instance2,
                     int split_index);
  void ExtendByOne(AdaptiveInstance& instance, const Graph& G, GRBEnv* env,
                   bool mip_subroutine = true);
  void ComputeAllObjectives(const Graph& G, AdaptiveInstance& instance,
                            bool compute_adaptive_objective = false);
  int policies() const { return policies_; }
  double worst_case_objective() const { return worst_case_objective_; }
  std::vector<std::vector<int>> partition() const { return partition_; }
  std::vector<Policy> solution() const { return solution_; }
  long solution_time() const { return solution_time_; }
  void set_policies(int policies) { policies_ = policies; }
  void set_scenarios(int scenarios) { scenarios_ = scenarios; }
  void set_worst_case_objective(double worst_case_objective) {
    worst_case_objective_ = worst_case_objective;
  }
  void set_solution_time(long solution_time) { solution_time_ = solution_time; }
  void set_solution_policy(int index, Policy policy) {
    solution_[index] = policy;
  }
  void set_partition(const std::vector<std::vector<int>>& partition) {
    partition_ = partition;
  }
  void set_unbounded(bool unbounded) { unbounded_ = unbounded; }
  void set_cuts(int cuts) { lazy_cuts_rounds_ = cuts; }
  void add_to_partition(int index, int scenario) {
    partition_[index].push_back(scenario);
  }

 private:
  bool unbounded_, benders_;
  int policies_, scenarios_, nodes_, arcs_;
  std::vector<std::vector<int>> partition_;
  std::vector<Policy> solution_;
  long solution_time_;
  double worst_case_objective_;
  int lazy_cuts_rounds_;
};

class RobustAlgoModel {
  // Special Purpose Model for the Enumerative Algorithm - Static Robust Dual
  // Reformulation - can also be used to solve non-robust shortest path
  // interdiction problems - just call update with only the desired scenario.
 public:
  RobustAlgoModel() : scenarios_(0), budget_(0), nodes_(0), arcs_(0){};
  RobustAlgoModel(AdaptiveInstance& instance, GRBEnv* env)
      : scenarios_(instance.scenarios()),
        budget_(instance.budget()),
        nodes_(instance.nodes()),
        arcs_(instance.arcs()),
        env_(env){};
  void ConfigureModel(const Graph& G, AdaptiveInstance& instance);
  void Update(std::vector<int>& subset);
  void ReverseUpdate(std::vector<int>& subset);
  Policy Solve();

 private:
  // All decision variables and non-negativity constraints will stay throughout
  // the algorithm. Constraints are added in ::Update based on the subset of [p]
  // that we want to include.
  int scenarios_, budget_, nodes_, arcs_;
  GRBEnv* env_;
  GRBModel* algo_model_;
  GRBVar z_var_;
  std::vector<std::vector<GRBVar>>
      pi_var_;  // Decision variable (for every q, i).
  std::vector<std::vector<GRBVar>>
      lambda_var_;             // Decision variable (for every q, a).
  std::vector<GRBVar> x_var_;  // Interdiction decision variable (for every a).
  std::vector<GRBTempConstr> z_constraints_;
  std::vector<std::vector<GRBTempConstr>> dual_constraints_;
};

class SetPartitioningModel {
  // Set Partitioning MIP for solving an adaptive instance.
 public:
  SetPartitioningModel()
      : big_m_(0),
        scenarios_(0),
        policies_(0),
        budget_(0),
        nodes_(0),
        arcs_(0){};
  SetPartitioningModel(int big_m, AdaptiveInstance& instance, GRBEnv* env)
      : big_m_(big_m),
        scenarios_(instance.scenarios()),
        policies_(instance.policies()),
        budget_(instance.budget()),
        nodes_(instance.nodes()),
        arcs_(instance.arcs()),
        env_(env){};
  void ConfigureSolver(const Graph& G, AdaptiveInstance& instance);
  AdaptiveSolution Solve();
  AdaptiveSolution current_solution() const { return current_solution_; }

 private:
  const int big_m_;
  int scenarios_, policies_, budget_, nodes_, arcs_;
  GRBEnv* env_;
  GRBModel* sp_model_;
  GRBVar z_var_;  // Piecewise linearization dummy variable z.
  std::vector<std::vector<GRBVar>>
      h_var_;  // Set partitioning variables h[w][q]==1 if q in subset w.
  std::vector<std::vector<GRBVar>>
      x_var_;  // Interdiction variable (for every w, a).
  std::vector<std::vector<std::vector<GRBVar>>>
      pi_var_;  // Dual decision variable (for every w, q, i) pi.
  std::vector<std::vector<std::vector<GRBVar>>>
      lambda_var_;  // Dual decision variable (for every w, q, a) lambda.
  AdaptiveSolution
      current_solution_;  // Solution updated whenever ::Solve() is called.
};

class BendersCallback : public GRBCallback {
 public:
  BendersCallback() : upper_bound_(DBL_MAX){};
  BendersCallback(int big_m, AdaptiveInstance& instance, GRBVar& z_var,
                  std::vector<std::vector<GRBVar>>& h_var,
                  std::vector<std::vector<GRBVar>>& x_var, GRBEnv* env)
      : upper_bound_(DBL_MAX),
        lower_bound_(DBL_MIN),
        big_m_(big_m),
        nodes_(instance.nodes()),
        arcs_(instance.arcs()),
        budget_(instance.budget()),
        scenarios_(instance.scenarios()),
        policies_(instance.policies()),
        arc_costs_(instance.arc_costs()),
        lazy_cuts_rounds_(0),
        interdiction_delta_(instance.interdiction_delta()),
        z_var_(z_var),
        h_var_(h_var),
        x_var_(x_var),
        env_(env){};
  void ConfigureSubModels(const Graph& G, AdaptiveInstance& instance);
  void SolveSubModels();
  int lazy_cuts_rounds() const { return lazy_cuts_rounds_; }

 protected:
  void callback();

 private:
  double epsilon_ = EPSILON;
  double upper_bound_, lower_bound_;
  int big_m_, nodes_, arcs_, budget_, scenarios_, policies_;
  std::vector<std::vector<int>> arc_costs_;
  int lazy_cuts_rounds_, interdiction_delta_;
  void ConfigureIndividualSubModel(const Graph& G, int w, int q);
  void UpdateSubModels(bool rev = false);
  void AddLazyCuts();

 public:
  GRBVar z_var_;  // Decision variable - objective value.
  std::vector<std::vector<GRBVar>>
      h_var_;  // Decision variable for every (w, q), set partitioning variable.
  std::vector<std::vector<GRBVar>>
      x_var_;  // Decision variable for every (w, a), interdiction policies.
  GRBEnv* env_;
  std::vector<std::vector<GRBModel*>>
      submodels_;  // Submodel for every (w, q) policy and scenario.
  std::vector<std::vector<std::vector<GRBVar>>>
      y_var_;  // Shortest path decision variable for every (w, q, a).
  AdaptiveSolution current_solution_;
};

class SetPartitioningBenders {
 public:
  SetPartitioningBenders() : big_m_(0){};
  SetPartitioningBenders(int big_m, AdaptiveInstance& instance, GRBEnv* env)
      : big_m_(big_m),
        nodes_(instance.nodes()),
        arcs_(instance.arcs()),
        budget_(instance.budget()),
        scenarios_(instance.scenarios()),
        policies_(instance.policies()),
        interdiction_delta_(instance.interdiction_delta()),
        env_(env){};
  void ConfigureSolver(const Graph& G, AdaptiveInstance& instance);
  AdaptiveSolution Solve();

 private:
  const int big_m_;
  int nodes_, arcs_, budget_, scenarios_, policies_, interdiction_delta_;
  GRBEnv* env_;
  GRBModel* benders_model_;
  GRBVar z_var_;  // Decision variable - objective function.
  std::vector<std::vector<GRBVar>>
      h_var_;  // Decision variable for every (w, q), set partitioning variable.
  std::vector<std::vector<GRBVar>>
      x_var_;  // Decision variable for every (w, a), interdiction policies.
  BendersCallback callback_;
};

AdaptiveSolution EnumSolve(AdaptiveInstance& instance, const Graph& G, GRBEnv* env);

std::vector<Policy> InitializeKPolicies(AdaptiveInstance& instance, const Graph& G,
                                        GRBEnv* env);

AdaptiveSolution KMeansHeuristic(AdaptiveInstance& instance, const Graph& G,
                                 GRBEnv* env);

long GetCurrentTime();

#endif
