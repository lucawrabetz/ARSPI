#ifndef SOLVERS_H
#define SOLVERS_H

#include <dirent.h>
#include <float.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "/home/luw28/gurobi950/linux64/include/gurobi_c++.h"
// #include "/home/luchino/gurobi1001/linux64/include/gurobi_c++.h"
// #include "/Library/gurobi902/mac64/include/gurobi_c++.h"
// #include "/Library/gurobi1001/macos_universal2/include/gurobi_c++.h"

typedef std::numeric_limits<double> dbl;

enum ASPI_Solver { MIP, BENDERS, ENUMERATION, GREEDY };

// Manual MIP symmetry constraint parameters.
constexpr int MANUAL_SYMMETRY_NONE = 0;  // Default.
constexpr int MANUAL_SYMMETRY_ASSIGNMENT = 1;
constexpr int MANUAL_SYMMETRY_NONDECREASING = 2;

// Gurobi MIP symmetry parameters.
constexpr int GUROBI_SYMMETRY_AUTO = -1;  // Default.
constexpr int GUROBI_SYMMETRY_NONE = 0;
constexpr int GUROBI_SYMMETRY_CONSERVATIVE = 1;
constexpr int GUROBI_SYMMETRY_AGGRESSIVE = 2;

constexpr double EPSILON = 0.000001;
constexpr int DEBUG = 0;

constexpr long TIME_LIMIT_S = 3600;
constexpr long TIME_LIMIT_MS = TIME_LIMIT_S * 1000;

const std::string DATA_DIRECTORY = "dat/";
const std::string OPTIMAL = "OPTIMAL";
const std::string NOT_OPTIMAL = "NOT_OPTIMAL";
const std::string UNBOUNDED = "UNBOUNDED";
const std::string NOT_UNBOUNDED = "NOT_UNBOUNDED";

const std::string INSTANCE_INFO_COLUMN_HEADERS =
    "set_name,instance_name,nodes,arcs,k_zero,density,scenarios,budget,"
    "policies";
// Partitions will be encoded as follows in the csv column:
// For every follower (indexed 0-p-1), we write the index (0-k-1) of the policy
// that they are assigned to, in a string delimeted by '-'. For example, if we
// have p=3 and k=2, with followers 0,1 in partition 0, and follower 2 in
// partition 1, we have: "0-0-1".
const std::string OUTPUT_COLUMN_HEADERS =
    "solver,unbounded,optimal,objective,gap,time,cuts_"
    "rounds,cuts_added,avg_cbtime,avg_sptime,"
    "partition,m_sym,g_sym";

class AdaptiveSolution;

struct InstanceInput {
  InstanceInput(const std::string& set_name, int nodes, int k_zero,
                int scenarios, int id)
      : nodes_(nodes),
        k_zero_(k_zero),
        scenarios_(scenarios),
        id_(id),
        setname_(set_name),
        graph_name_(setname_ + "-" + std::to_string(nodes) + "_" +
                    std::to_string(k_zero)),
        costs_name_(graph_name_ + "-costs_" + std::to_string(scenarios) + "_" +
                    std::to_string(id)) {}
  std::string GraphFileName() const {
    return DATA_DIRECTORY + setname_ + "/" + graph_name_ + ".txt";
  }
  std::string CostFileName() const {
    return DATA_DIRECTORY + setname_ + "/" + costs_name_ + ".csv";
  }
  const int nodes_, k_zero_, scenarios_, id_;
  const std::string setname_;
  const std::string graph_name_;
  const std::string
      costs_name_;  // InstanceInput.costs_name_ and AdaptiveInstance.name_ are
                    // the same, except for that InstanceInput.costs_name_ has
                    // the word "costs" in it to be quickly used in
                    // InstanceInput::CostFileName() to get the right cost file
                    // path.
};

class Graph {
  // Class to read an arc list from a file and represent it as linked lists.
 public:
  Graph();
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
  // Arc vectors - adjacency_list_ is a classic linked list representation of
  // the arc list (outgoing arcs). The vector arc_index_hash_ directly maps
  // every arc (i, j), (so adjacency_list[i] includes j) to its index a in
  // 0,...,m-1. adjacency_list[i] and arc_index_hash[i] are aligned. For
  // example, if node i = 0 has 2 outgoing arcs, to nodes j = 1, 2, with arc
  // indices a = 0, 1, respectively, we will have:
  //    adjacency_list[0] = {1, 2}
  //    arc_index_hash[0] = {0, 1}
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
  AdaptiveInstance(const InstanceInput& instance_input)
      : nodes_(instance_input.nodes_),
        scenarios_(instance_input.scenarios_),
        name_(instance_input.graph_name_ + "-" +
              std::to_string(instance_input.scenarios_) + "_" +
              std::to_string(instance_input.id_)),
        costs_filename_(instance_input.CostFileName()){};
  // Copy constructor with a different U (only a subset of scenarios to keep).
  AdaptiveInstance(const AdaptiveInstance& adaptive_instance,
                   std::vector<int>& keep_scenarios);
  void ReadCosts();
  void PrintInstance(const Graph& G) const;
  double SPDijkstra(int q, const Graph& G) const;
  double SPModel(int q, const Graph& G, GRBEnv* env) const;
  double SolveASPIZeroPolicies(const Graph& G, GRBEnv* env) const;
  void ApplyInterdiction(const std::vector<double>& x_bar, bool rev = false);
  double ValidatePolicy(std::vector<double>& x_bar, const Graph& G,
                        GRBEnv* env);
  int nodes() const { return nodes_; }
  // int arcs() const { return arcs_; }
  int scenarios() const { return scenarios_; }
  // int policies() const { return policies_; }
  // int budget() const { return budget_; }
  int interdiction_delta() const { return interdiction_delta_; }
  int big_m() const { return big_m_; }
  std::string name() const { return name_; }
  std::vector<std::vector<int>> arc_costs() const { return arc_costs_; }
  std::vector<int> scenario_index_map() const { return scenario_index_map_; }
  // No setter for scenarios_ - functionality reserved for change copy
  // constructor.
  // void set_policies(int policies) { policies_ = policies; }

 private:
  int nodes_, /*arcs_,*/ scenarios_ /*, policies_, budget_*/;
  int interdiction_delta_, big_m_;
  const std::string name_;
  const std::string costs_filename_;
  std::vector<std::vector<int>> arc_costs_;
  // Map of scenario indices (indexed 0-(p-1)) to their original indices, for
  // when an AdaptiveInstance is copied but only keeps a subset of U (to be able
  // to recover the original). An empty map indicates an original instance.
  std::vector<int> scenario_index_map_;
};

class SingleRunInput {
 public:
  SingleRunInput(const InstanceInput& instance_input, int policies, int budget,
                 GRBEnv* env, int manual_symmetry_constraints,
                 int gurobi_symmetry_detection, double greedy_mip_gap_threshold)
      : G_(Graph(instance_input.GraphFileName(), instance_input.nodes_)),
        instance_(AdaptiveInstance(instance_input)),
        policies_(policies),
        budget_(budget),
        k_zero_(instance_input.k_zero_),
        manual_symmetry_constraints_(manual_symmetry_constraints),
        gurobi_symmetry_detection_(gurobi_symmetry_detection),
        greedy_mip_gap_threshold_(greedy_mip_gap_threshold),
        env_(env) {
    instance_.ReadCosts();
  }
  SingleRunInput(const SingleRunInput& problem_input,
                 std::vector<int>& keep_scenarios, int policies)
      : G_(problem_input.G_),
        instance_(AdaptiveInstance(problem_input.instance_, keep_scenarios)),
        policies_(policies),
        budget_(problem_input.budget_),
        k_zero_(problem_input.k_zero_),
        manual_symmetry_constraints_(
            problem_input.manual_symmetry_constraints_),
        gurobi_symmetry_detection_(problem_input.gurobi_symmetry_detection_),
        greedy_mip_gap_threshold_(problem_input.greedy_mip_gap_threshold_),
        env_(problem_input.env_) {}
  std::string InputCSVString(const std::string& set_name) const;
  AdaptiveSolution SolveASPIZeroPolicies() const;
  const Graph G_;
  AdaptiveInstance instance_;
  int policies_, budget_, k_zero_;
  int manual_symmetry_constraints_;
  int gurobi_symmetry_detection_;
  double greedy_mip_gap_threshold_;
  GRBEnv* env_;
};

struct BendersMetaStats {
  BendersMetaStats()
      : number_of_callbacks(0),
        number_of_spcalls(0),
        total_lazy_cuts_added(0),
        total_sp_separation_time(0),
        total_callback_time(0){};
  int number_of_callbacks;
  int number_of_spcalls;
  int total_lazy_cuts_added;
  long total_sp_separation_time;  // Time in ms.
  long total_callback_time;       // Time in ms.
};

class AdaptiveSolution {
  // Full solution - i.e. a vector of k policy objects with extra info,
  // including the partition.
  // TODO: are all the constructors needed?
  // TODO: default initializations of stuff for different solvers (approximation
  // and enumeration algorithm).
 public:
  AdaptiveSolution() = delete;
  AdaptiveSolution(ASPI_Solver solver, bool unbounded, bool optimal)
      : solver_(solver),
        unbounded_(unbounded),
        optimal_(optimal),
        number_of_callbacks_(0),
        total_lazy_cuts_added_(0),
        avg_callback_time_(0),
        avg_sp_separation_time_(0){};
  AdaptiveSolution(ASPI_Solver solver, bool optimal,
                   const SingleRunInput& problem)
      : solver_(solver),
        unbounded_(false),
        optimal_(optimal),
        policies_(problem.policies_),
        scenarios_(problem.instance_.scenarios()),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        budget_(problem.budget_),
        partition_(std::vector<std::vector<int>>(problem.policies_)),
        solution_(std::vector<std::vector<double>>(
            problem.policies_, std::vector<double>(problem.G_.arcs()))),
        objectives_(std::vector<std::vector<double>>(
            problem.policies_,
            std::vector<double>(problem.instance_.scenarios()))),
        worst_case_objective_(-1),
        mip_gap_(-1),
        number_of_callbacks_(0),
        total_lazy_cuts_added_(0),
        avg_callback_time_(0),
        avg_sp_separation_time_(0){};
  AdaptiveSolution(ASPI_Solver solver, bool optimal,
                   const SingleRunInput& problem,
                   const std::vector<std::vector<int>>& partition,
                   const std::vector<std::vector<double>>& solution)
      : solver_(solver),
        unbounded_(false),
        optimal_(optimal),
        policies_(problem.policies_),
        scenarios_(problem.instance_.scenarios()),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        budget_(problem.budget_),
        partition_(partition),
        solution_(solution),
        objectives_(std::vector<std::vector<double>>(
            problem.policies_,
            std::vector<double>(problem.instance_.scenarios()))),
        worst_case_objective_(-1),
        mip_gap_(-1),
        number_of_callbacks_(0),
        total_lazy_cuts_added_(0),
        avg_callback_time_(0),
        avg_sp_separation_time_(0){};
  void LogSolution(const SingleRunInput& problem, bool policy = false);
  void MergeEnumSols(AdaptiveSolution sol2, AdaptiveInstance* instance2,
                     int split_index);
  void ExtendByOne(AdaptiveInstance& instance, const Graph& G, GRBEnv* env,
                   bool mip_subroutine = true);
  void ComputeObjectiveMatrix(const SingleRunInput& problem);
  void SetObjectiveMatrix(
      const std::vector<std::vector<double>>& objective_matrix) {
    objectives_ = objective_matrix;
  }
  void ComputePartition();
  void ComputeAdaptiveObjective();
  ASPI_Solver solver() const { return solver_; }
  int policies() const { return policies_; }
  bool optimal() const { return optimal_; }
  double worst_case_objective() const { return worst_case_objective_; }
  double mip_gap() const { return mip_gap_; }
  std::vector<std::vector<int>> partition() const { return partition_; }
  std::vector<std::vector<double>> solution() const { return solution_; }
  long solution_time() const { return solution_time_; }
  int number_of_callbacks() const { return number_of_callbacks_; }
  int total_lazy_cuts_added() const { return total_lazy_cuts_added_; }
  long avg_callback_time() const { return avg_callback_time_; }
  long avg_sp_separation_time() const { return avg_sp_separation_time_; }
  void set_policies(int policies) { policies_ = policies; }
  void set_scenarios(int scenarios) { scenarios_ = scenarios; }
  void set_worst_case_objective(double worst_case_objective) {
    worst_case_objective_ = worst_case_objective;
  }
  void set_solution_time(long solution_time) { solution_time_ = solution_time; }
  void set_solution_policy(int index, std::vector<double>& policy) {
    solution_[index] = policy;
  }
  void set_optimal(bool optimal) {
    optimal_ = optimal;
    mip_gap_ = 0;
  }
  void set_mip_gap(double gap) { mip_gap_ = gap; }
  void set_partition(const std::vector<std::vector<int>>& partition) {
    partition_ = partition;
  }
  void set_unbounded(bool unbounded) { unbounded_ = unbounded; }
  void set_benders_stats(const BendersMetaStats& stats) {
    number_of_callbacks_ = stats.number_of_callbacks;
    total_lazy_cuts_added_ = stats.total_lazy_cuts_added;
    if (number_of_callbacks_ <= 0) {
      avg_callback_time_ = -1;
    } else {
      avg_callback_time_ = stats.total_callback_time / number_of_callbacks_;
    }
    if (stats.number_of_spcalls <= 0) {
      avg_sp_separation_time_ = -1;
    } else {
      avg_sp_separation_time_ =
          stats.total_sp_separation_time / stats.number_of_spcalls;
    }
  }
  std::string OutputCSVString(const SingleRunInput& problem) const;
  std::string SingleRunLogLine(const SingleRunInput& problem) const;
  std::string PartitionString(const SingleRunInput& problem) const;
  void add_to_partition(int index, int scenario) {
    partition_[index].push_back(scenario);
  }

 private:
  ASPI_Solver solver_;
  bool unbounded_, benders_, time_limit_, optimal_;
  int policies_, scenarios_, nodes_, arcs_, budget_;
  std::vector<std::vector<int>> partition_;
  std::vector<std::vector<double>> solution_;
  std::vector<std::vector<double>>
      objectives_;  // k by p matrix, objective value of each interdiction
                    // policy applied to each follower [w][q].
  double worst_case_objective_;
  double mip_gap_;
  long solution_time_;
  int number_of_callbacks_;
  int total_lazy_cuts_added_;
  long avg_callback_time_;
  long avg_sp_separation_time_;
};

class RobustAlgoModel {
  // Special Purpose Model for the Enumerative Algorithm - Static Robust Dual
  // Reformulation - can also be used to solve non-robust shortest path
  // interdiction problems - just call update with only the desired scenario.
 public:
  ~RobustAlgoModel() { delete algo_model_; }
  RobustAlgoModel() : scenarios_(0), budget_(0), nodes_(0), arcs_(0){};
  RobustAlgoModel(const SingleRunInput& problem)
      : scenarios_(problem.instance_.scenarios()),
        budget_(problem.budget_),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        env_(problem.env_){};
  void ConfigureModel(const SingleRunInput& problem);
  void Update(std::vector<int>& subset);
  void ReverseUpdate(std::vector<int>& subset);
  std::pair<double, std::vector<double>> Solve();

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
  ~SetPartitioningModel() { delete sp_model_; }
  SetPartitioningModel()
      : big_m_(0),
        scenarios_(0),
        policies_(0),
        budget_(0),
        nodes_(0),
        arcs_(0){};
  SetPartitioningModel(const SingleRunInput& problem)
      : big_m_(problem.instance_.big_m()),
        scenarios_(problem.instance_.scenarios()),
        policies_(problem.policies_),
        budget_(problem.budget_),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        env_(problem.env_),
        sp_model_(new GRBModel(env_)),
        z_var_(sp_model_->addVar(0, big_m_, -1, GRB_CONTINUOUS, "z")),
        h_var_(std::vector<std::vector<GRBVar>>(
            policies_, std::vector<GRBVar>(scenarios_))),
        x_var_(std::vector<std::vector<GRBVar>>(policies_,
                                                std::vector<GRBVar>(arcs_))),
        pi_var_(std::vector<std::vector<std::vector<GRBVar>>>(
            policies_, std::vector<std::vector<GRBVar>>(
                           scenarios_, std::vector<GRBVar>(nodes_)))),
        lambda_var_(std::vector<std::vector<std::vector<GRBVar>>>(
            policies_, std::vector<std::vector<GRBVar>>(
                           scenarios_, std::vector<GRBVar>(arcs_)))){};
  void ConfigureSolver(const SingleRunInput& problem);
  AdaptiveSolution Solve(const SingleRunInput& problem);

 private:
  void AddSetPartitioningVariables();
  void AddInterdictionPolicyVariables();
  void AddPiShortestPathVariable();
  void AddDualLambdaArcVariable();
  void AddBudgetConstraints();
  void AddObjectiveBoundingConstraint();
  void AddArcCostBoundingConstraint(const SingleRunInput& problem);
  void AddSetPartitioningConstraints();
  void AddRootNodeZeroConstraint();
  void AddAssignmentSymmetryConstraints();
  void AddNonDecreasingSymmetryConstraints();
  void ProcessInputSymmetryParameters(const SingleRunInput& problem);
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
};

struct ShortestPath {
  std::vector<int> arcs;
  double cost;
};

class BendersCallback : public GRBCallback {
 public:
  ~BendersCallback() = default;
  BendersCallback() : upper_bound_(DBL_MAX){};
  BendersCallback(const SingleRunInput& problem, GRBVar& z_var,
                  std::vector<std::vector<GRBVar>>& h_var,
                  std::vector<std::vector<GRBVar>>& x_var)
      : upper_bound_(DBL_MAX),
        lower_bound_(DBL_MIN),
        big_m_(problem.instance_.big_m()),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        scenarios_(problem.instance_.scenarios()),
        policies_(problem.policies_),
        budget_(problem.budget_),
        interdiction_delta_(problem.instance_.interdiction_delta()),
        arc_costs_(problem.instance_.arc_costs()),
        sparse_x_var_(policies_, std::vector<int>(budget_, -1)),
        clusters_h_var_(policies_),
        adjacency_list_(problem.G_.adjacency_list()),
        arc_index_hash_(problem.G_.arc_index_hash()),
        shortest_paths_(scenarios_),
        stats_(BendersMetaStats()),
        z_var_(z_var),
        h_var_(h_var),
        x_var_(x_var),
        env_(problem.env_) {
    for (std::vector<int>& cluster : clusters_h_var_) {
      cluster.reserve(scenarios_);
    }
    for (auto& p : shortest_paths_) {
      p.arcs.reserve(nodes_);
      p.cost = -1;
    }
  };
  void UpdateMasterVariables();
  void UpdateArcCostsForQWPair(int w, int q, bool rev = false);
  void SPDijkstra(int q);
  void ConfigureSubModels(const SingleRunInput& problem);
  void SolveSubModels();
  void ComputeInitialPaths();
  std::vector<ShortestPath> paths() { return shortest_paths_; }
  BendersMetaStats meta_stats() const { return stats_; }

 protected:
  void callback();

 private:
  double epsilon_ = EPSILON;
  double upper_bound_, lower_bound_;
  int big_m_, nodes_, arcs_, scenarios_, policies_, budget_,
      interdiction_delta_;
  int iteration_ = 0;
  std::vector<std::vector<int>> arc_costs_;
  std::vector<std::vector<int>> sparse_x_var_;
  std::vector<std::vector<int>> clusters_h_var_;
  std::vector<std::vector<int>> adjacency_list_;
  std::vector<std::vector<int>> arc_index_hash_;
  // Shortest paths for each follower (interdiction policy assigned to it by H
  // will be used). each path (inner vector) - lists the arc indices in the
  // path.
  // shortest_paths_[q].first - path represented by arc indices.
  // shortest_paths_[q].second - cost of path.
  // THESE PATHS ARE NOT IN ORDER - THEY ARE JUST THE LIST OF ARCS INCLUDED IN
  // THE PATH, AS WE JUST NEED THAT TO ADD THEM AS LAZY CUTS.
  std::vector<ShortestPath> shortest_paths_;
  BendersMetaStats stats_;
  void ConfigureIndividualSubModel(const Graph& G, int w, int q);
  void UpdateSubModels(bool rev = false);
  void AddLazyCuts();

 public:
  // Decision variables from master problem.
  GRBVar z_var_;  // Decision variable - objective value.
  std::vector<std::vector<GRBVar>>
      h_var_;  // Decision variable for every (w, q), set partitioning variable.
  std::vector<std::vector<GRBVar>>
      x_var_;  // Decision variable for every (w, a), interdiction policies.
  GRBEnv* env_;
};

class SetPartitioningBenders {
 public:
  ~SetPartitioningBenders() { delete benders_model_; }
  SetPartitioningBenders() : big_m_(0){};
  SetPartitioningBenders(const SingleRunInput& problem)
      : big_m_(problem.instance_.big_m()),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        budget_(problem.budget_),
        scenarios_(problem.instance_.scenarios()),
        policies_(problem.policies_),
        interdiction_delta_(problem.instance_.interdiction_delta()),
        env_(problem.env_),
        benders_model_(new GRBModel(env_)),
        z_var_(benders_model_->addVar(0, big_m_, -1, GRB_CONTINUOUS, "z")),
        h_var_(std::vector<std::vector<GRBVar>>(
            policies_, std::vector<GRBVar>(scenarios_))),
        x_var_(std::vector<std::vector<GRBVar>>(policies_,
                                                std::vector<GRBVar>(arcs_))){};
  void ConfigureSolver(const SingleRunInput& problem);
  AdaptiveSolution Solve(const SingleRunInput& problem);
  void SetMIPGap(double mip_gap) {
    benders_model_->set(GRB_DoubleParam_MIPGap, mip_gap);
  }

 private:
  void AddSetPartitioningVariables();
  void AddInterdictionPolicyVariables();
  void AddBudgetConstraints();
  void AddSetPartitioningConstraints();
  void AddAssignmentSymmetryConstraints();
  void AddNonDecreasingSymmetryConstraints();
  void ProcessInputSymmetryParameters(const SingleRunInput& problem);
  void AddInitialPathConstraints(const std::vector<ShortestPath>& paths);
  void InitialConstraints();
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

AdaptiveSolution EnumSolve(SingleRunInput& problem);

AdaptiveSolution GreedyAlgorithm(SingleRunInput& problem);

// Experiment main loops.

void SingleRunOnAllInstancesInSetDirectory(
    const int min_policies, const int max_policies, const int min_budget,
    const int max_budget, const std::string& set_name,
    std::ofstream& result_file, const ASPI_Solver& solver,
    int manual_symmetry_constraints, int gurobi_symmetry_detection,
    double greedy_mip_gap_threshold);

void TestOnAllInstancesInSetDirectory(const std::string& set_name,
                                      int manual_symmetry_constraints,
                                      int gurobi_symmetry_detection,
                                      double greedy_mip_gap_threshold);

void UninterdictedRunOnAllInstancesInsetDirectory(const std::string& set_name,
                                                  std::ofstream& result_file);

#endif
