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
#include <limits>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// #include "/home/luw28/gurobi950/linux64/include/gurobi_c++.h"
#include "/home/luchino/gurobi1001/linux64/include/gurobi_c++.h"
// #include "/Library/gurobi902/mac64/include/gurobi_c++.h"

typedef std::numeric_limits<double> dbl;

enum ASPI_Solver { MIP, BENDERS, ENUMERATION, GREEDY };

const double EPSILON = 0.000001;
const int DEBUG = 2;

const long TIME_LIMIT_S = 3600;
const long TIME_LIMIT_MS = TIME_LIMIT_S * 1000;

const std::string DATA_DIRECTORY = "dat/";

const std::string INSTANCE_INFO_COLUMN_HEADERS =
    "set_name,instance_name,nodes,arcs,k_zero,density,scenarios,budget,"
    "policies";
// Partitions will be encoded as follows in the csv column:
// For every follower (indexed 0-p-1), we write the index (0-k-1) of the policy
// that they are assigned to, in a string delimeted by '-'. For example, if we
// have p=3 and k=2, with followers 0,1 in partition 0, and follower 2 in
// partition 1, we have: "0-0-1".
const std::string MIP_COLUMN_HEADERS =
    "MIP_OPTIMAL,MIP_objective,MIP_gap,MIP_time,MIP_partition";
const std::string BENDERS_COLUMN_HEADERS =
    "BENDERS_OPTIMAL,BENDERS_objective,BENDERS_gap,BENDERS_time,BENDERS_cuts_"
    "rounds,BENDERS_partition";
const std::string ENUMERATION_COLUMN_HEADERS =
    "ENUMERATION_OPTIMAL,ENUMERATION_objective,ENUMERATION_time,ENUMERATION_"
    "partition";
const std::string GREEDY_COLUMN_HEADERS =
    "GREEDY_objective,GREEDY_time,GREEDY_partition";

struct GraphInput {
  GraphInput(const std::string& setname, int nodes, int k_zero)
      : setname_(setname),
        nodes_(nodes),
        k_zero_(k_zero),
        graph_name_(setname + "-" + std::to_string(nodes) + "_" +
                    std::to_string(k_zero)) {}
  std::string FileName() const {
    return DATA_DIRECTORY + setname_ + "/" + graph_name_ + ".txt";
  }
  const std::string setname_;
  const int nodes_, k_zero_;
  const std::string graph_name_;
};

struct InstanceInput {
  InstanceInput(const GraphInput& graph_input, int scenarios, int id)
      : graph_input_(graph_input),
        scenarios_(scenarios),
        id_(id),
        setname_(graph_input.setname_),
        costs_name_(graph_input_.graph_name_ + "-costs_" +
                    std::to_string(scenarios) + "_" + std::to_string(id)) {}
  std::string CostFileName() const {
    return DATA_DIRECTORY + setname_ + "/" + costs_name_ + ".csv";
  }
  const GraphInput graph_input_;
  const int scenarios_, id_;
  const std::string setname_;
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
  AdaptiveInstance(const InstanceInput& instance_input)
      : nodes_(instance_input.graph_input_.nodes_),
        scenarios_(instance_input.scenarios_),
        name_(instance_input.graph_input_.graph_name_ + "-" +
              std::to_string(instance_input.scenarios_) + "_" +
              std::to_string(instance_input.id_)),
        costs_filename_(instance_input.CostFileName()){};
  // Copy constructor with a different U (only a subset of scenarios to keep).
  AdaptiveInstance(AdaptiveInstance& adaptive_instance,
                   std::vector<int>& keep_scenarios);
  void ReadCosts();
  void PrintInstance(const Graph& G) const;
  double SPModel(int q, const Graph& G, GRBEnv* env);
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

class ProblemInput {
 public:
  ProblemInput(const InstanceInput& instance_input, int policies, int budget,
               GRBEnv* env)
      : G_(Graph(instance_input.graph_input_.FileName(),
                 instance_input.graph_input_.nodes_)),
        instance_(AdaptiveInstance(instance_input)),
        policies_(policies),
        budget_(budget),
        k_zero_(instance_input.graph_input_.k_zero_),
        env_(env) {
    instance_.ReadCosts();
  }
  ProblemInput(ProblemInput& problem_input, std::vector<int>& keep_scenarios,
               int policies)
      : G_(problem_input.G_),
        instance_(AdaptiveInstance(problem_input.instance_, keep_scenarios)),
        policies_(policies),
        budget_(problem_input.budget_),
        k_zero_(problem_input.k_zero_),
        env_(problem_input.env_) {}
  void WriteLineToLogFile();
  const Graph G_;
  AdaptiveInstance instance_;
  int policies_, budget_, k_zero_;
  GRBEnv* env_;
};

class AdaptiveSolution {
  // Full solution - i.e. a vector of k policy objects with extra info,
  // including the partition.
  // TODO: are all the constructors needed?
  // TODO: default initializations of stuff for different solvers (approximation
  // and enumeration algorithm).
 public:
  AdaptiveSolution() : unbounded_(false), benders_(false), optimal_(false){};
  AdaptiveSolution(bool unbounded, bool benders, bool optimal)
      : unbounded_(unbounded), benders_(benders), optimal_(optimal){};
  AdaptiveSolution(bool benders, bool optimal, const ProblemInput& problem)
      : unbounded_(false),
        benders_(benders),
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
        mip_gap_(-1){};
  AdaptiveSolution(bool benders, bool optimal, const ProblemInput& problem,
                   const std::vector<std::vector<int>>& partition,
                   const std::vector<std::vector<double>>& solution)
      : unbounded_(false),
        benders_(benders),
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
        mip_gap_(-1){};
  void LogSolution(const ProblemInput& problem, bool policy = false);
  void MergeEnumSols(AdaptiveSolution sol2, AdaptiveInstance* instance2,
                     int split_index);
  void ExtendByOne(AdaptiveInstance& instance, const Graph& G, GRBEnv* env,
                   bool mip_subroutine = true);
  void ComputeObjectiveMatrix(const ProblemInput& problem);
  void SetObjectiveMatrix(
      const std::vector<std::vector<double>>& objective_matrix) {
    objectives_ = objective_matrix;
  }
  void ComputePartition();
  void ComputeAdaptiveObjective();
  int policies() const { return policies_; }
  int lazy_cuts_rounds() const { return lazy_cuts_rounds_; }
  bool optimal() const { return optimal_; }
  double worst_case_objective() const { return worst_case_objective_; }
  double mip_gap() const { return mip_gap_; }
  std::vector<std::vector<int>> partition() const { return partition_; }
  std::vector<std::vector<double>> solution() const { return solution_; }
  long solution_time() const { return solution_time_; }
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
  void set_cuts(int cuts) { lazy_cuts_rounds_ = cuts; }
  void add_to_partition(int index, int scenario) {
    partition_[index].push_back(scenario);
  }

 private:
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
  int lazy_cuts_rounds_;
};

class RobustAlgoModel {
  // Special Purpose Model for the Enumerative Algorithm - Static Robust Dual
  // Reformulation - can also be used to solve non-robust shortest path
  // interdiction problems - just call update with only the desired scenario.
 public:
  RobustAlgoModel() : scenarios_(0), budget_(0), nodes_(0), arcs_(0){};
  RobustAlgoModel(const ProblemInput& problem)
      : scenarios_(problem.instance_.scenarios()),
        budget_(problem.budget_),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        env_(problem.env_){};
  void ConfigureModel(const ProblemInput& problem);
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
  SetPartitioningModel()
      : big_m_(0),
        scenarios_(0),
        policies_(0),
        budget_(0),
        nodes_(0),
        arcs_(0){};
  SetPartitioningModel(const ProblemInput& problem)
      : big_m_(problem.instance_.big_m()),
        scenarios_(problem.instance_.scenarios()),
        policies_(problem.policies_),
        budget_(problem.budget_),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        env_(problem.env_){};
  void ConfigureSolver(const ProblemInput& problem);
  AdaptiveSolution Solve(const ProblemInput& problem);

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
};

class BendersCallback : public GRBCallback {
 public:
  BendersCallback() : upper_bound_(DBL_MAX){};
  BendersCallback(const ProblemInput& problem, GRBVar& z_var,
                  std::vector<std::vector<GRBVar>>& h_var,
                  std::vector<std::vector<GRBVar>>& x_var)
      : upper_bound_(DBL_MAX),
        lower_bound_(DBL_MIN),
        big_m_(problem.instance_.big_m()),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        scenarios_(problem.instance_.scenarios()),
        policies_(problem.policies_),
        arc_costs_(problem.instance_.arc_costs()),
        lazy_cuts_rounds_(0),
        interdiction_delta_(problem.instance_.interdiction_delta()),
        z_var_(z_var),
        h_var_(h_var),
        x_var_(x_var),
        env_(problem.env_){};
  void ConfigureSubModels(const ProblemInput& problem);
  void SolveSubModels();
  int lazy_cuts_rounds() const { return lazy_cuts_rounds_; }

 protected:
  void callback();

 private:
  double epsilon_ = EPSILON;
  double upper_bound_, lower_bound_;
  int big_m_, nodes_, arcs_, scenarios_, policies_;
  int iteration_ = 0;
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
};

class SetPartitioningBenders {
 public:
  SetPartitioningBenders() : big_m_(0){};
  SetPartitioningBenders(const ProblemInput& problem)
      : big_m_(problem.instance_.big_m()),
        nodes_(problem.instance_.nodes()),
        arcs_(problem.G_.arcs()),
        budget_(problem.budget_),
        scenarios_(problem.instance_.scenarios()),
        policies_(problem.policies_),
        interdiction_delta_(problem.instance_.interdiction_delta()),
        env_(problem.env_){};
  void ConfigureSolver(const ProblemInput& problem);
  AdaptiveSolution Solve(const ProblemInput& problem);

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

AdaptiveSolution EnumSolve(ProblemInput& problem);

AdaptiveSolution GreedyAlgorithm(ProblemInput& problem);

long GetCurrentTime();

// Experiment helper functions.

bool IsCostFile(const std::string& name);

bool IsRunFile(const std::string& name);

int DigitsToInt(const std::vector<int> num_digits);

std::pair<int, int> GetNodesKZero(const std::string& name, int start);

std::pair<int, int> GetScenariosID(const std::string& name);

std::string SolutionPartitionToString(const AdaptiveSolution& solution,
                                      const ProblemInput& problem);

std::string SolveAndPrintTest(const std::string& set_name,
                              const ProblemInput& problem,
                              ProblemInput& problem_copyable,
                              const std::vector<ASPI_Solver>& solvers,
                              int debug = 0);

void RunAllInstancesInSetDirectory(const int min_policies,
                                   const int max_policies, const int min_budget,
                                   const int max_budget,
                                   const std::string& set_name,
                                   const std::vector<ASPI_Solver>& solvers);

#endif
