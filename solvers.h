#ifndef SOLVERS_H
#define SOLVERS_H

#include <climits>
#include <utility>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <cmath>
#include <float.h>
#include <unordered_set>
#include <queue>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <vector>
#include <algorithm>
#include <array>
#include <chrono>
// #include "/home/luw28/gurobi950/linux64/include/gurobi_c++.h"
#include "/home/luchino/gurobi1001/linux64/include/gurobi_c++.h"
// #include "/Library/gurobi902/mac64/include/gurobi_c++.h"

using std::array;
using std::shuffle;
using std::default_random_engine;
using std::abs;
using std::random_device;
using std::mt19937;
using std::uniform_int_distribution;
using std::normal_distribution;
using std::bernoulli_distribution;
using std::pair;
using std::make_pair;
using std::stringstream;
using std::endl;
using std::cout;
using std::string;
using std::vector;
using std::unordered_set;
using std::to_string;
using std::stoi;
using std::ifstream;
using std::ofstream;
using std::round;

// WHY:
// Can we replace some explicit setters with with "updaters" that set internally (e.g. Policy.objective)

class Graph
{
    // Class to read an arc list from a file and represent it as linked lists.
public:
    Graph(const string &filename, int nodes);
    void PrintArc(int a, int i, int index, bool rev=false) const;
    void PrintGraph(bool rev=false) const;
    void PrintGraphWithCosts(const vector<vector<int>>& costs, int interdiction_delta) const;
    // Getters.
    int nodes() const {return nodes_;}
    int arcs() const {return arcs_;}
    vector<vector<int>> arc_index_hash() const {return arc_index_hash_;}
    vector<vector<int>> adjacency_list() const {return adjacency_list_;}
    vector<vector<int>> rev_arc_index_hash() const {return rev_arc_index_hash_;}
    vector<vector<int>> rev_adjacency_list() const {return rev_adjacency_list_;}
private:
    int nodes_, arcs_;
    const string filename_;
    // Arc vectors - adjacency_list_ is a linked list representation of the arc list (outgoing arcs).
    // The vector arc_index_hash_ directly maps every arc j = adjacency_list_[i] to its index a in 0,...,m-1.
    vector<vector<int>> arc_index_hash_; 
    vector<vector<int>> adjacency_list_;
    // Similarly, rev_adjacency_list and rev_arc_index_hash represent a linked list representation of incoming arcs.
    vector<vector<int>> rev_arc_index_hash_;
    vector<vector<int>> rev_adjacency_list_;
};

class AdaptiveInstance
{
    // Full adaptive problem instance with cost data, (Graph G is always passed by reference separately).
public:
    AdaptiveInstance(int scenarios, int policies, int budget, const Graph &G, const string &directory, const string &name) :
        nodes_(G.nodes()), arcs_(G.arcs()), scenarios_(scenarios), policies_(policies), budget_(budget), directory_(directory), name_(name) {};
    // Copy constructor with a different U (only a subset of scenarios to keep).
    AdaptiveInstance(AdaptiveInstance* adaptive_instance, vector<int>& keep_scenarios);
    void ReadCosts(int interdiction_delta);
    void PrintInstance(const Graph &G) const;
    int Dijkstra(int q, const Graph &G);
    void ApplyInterdiction(const vector<double>& x_bar, bool rev=false);
    double ValidatePolicy(vector<double>& x_bar, const Graph& G);
    double ComputeObjectiveOfPolicyForScenario(const vector<double>& binary_policy, const Graph& G, int q); 
    // Getters.
    int nodes() const {return nodes_;}
    int arcs() const {return arcs_;}
    int scenarios() const {return scenarios_;}
    int policies() const {return policies_;}
    int budget() const {return budget_;}
    int interdiction_delta() const {return interdiction_delta_;}
    vector<vector<int>> arc_costs() const {return arc_costs_;} // ?
    vector<int> scenario_index_map() const {return scenario_index_map_;}
    // Setters.
    // No setter for scenarios_ - functionality reserved for change copy constructor.
    void set_policies(int policies){policies_ = policies;}
private:
    int nodes_, arcs_, scenarios_, policies_, budget_;
    int interdiction_delta_;
    const string directory_;
    const string name_;
    vector<vector<int>> arc_costs_;
    // Map of scenario indices (indexed 0-(p-1)) to their original indices, for when an AdaptiveInstance
    // is copied but only keeps a subset of U (to be able to recover the original). An empty
    // map indicates an original instance.
    vector<int> scenario_index_map_;
};

class Policy
{
    // Policy class - represent a single interdiction policy with its objective.
    // The objective - is the worst (max) objective for the subsets that the interdiction policy is assigned to.
public:
    Policy() : size_(0), objective_(0), binary_policy_(vector<double>()) {};
    Policy(int size) : size_(size), objective_(0), binary_policy_(vector<double>(size)) {};
    Policy(int size, double objective, const vector<double>& binary_policy) : size_(size), objective_(objective), binary_policy_(binary_policy) {};
    // Getters.
    double objective() const {return objective_;}
    vector<double> binary_policy() const {return binary_policy_;}
    // Setters.
    void set_policy(const vector<double>& binary_policy) {binary_policy_=binary_policy;}
    void set_objective(double objective) {objective_=objective;}
private:
    int size_;
    double objective_;
    vector<double> binary_policy_;
};

class AdaptiveSolution
{
    // Full solution - i.e. a vector of k policy objects with extra info, including the partition.
public:
    AdaptiveSolution() : policies_(0), scenarios_(0), partition_(vector<vector<int>>(0)), solution_(vector<Policy>(0)){};
    AdaptiveSolution(bool unbounded) : unbounded_(unbounded){};
    AdaptiveSolution(int policies, int scenarios) : policies_(policies), scenarios_(scenarios), partition_(vector<vector<int>>(policies)), solution_(vector<Policy>(policies)) {};
    AdaptiveSolution(int policies, int scenarios, const vector<vector<int>>& partition, const vector<Policy>& solution) : policies_(policies), scenarios_(scenarios), partition_(partition), solution_(solution), unbounded_(false) {};
    void LogSolution(const Graph& G, AdaptiveInstance& m3, string title, bool policy=true);
    void MergeEnumSols(AdaptiveSolution sol2, AdaptiveInstance* instance2, int split_index);
    void ExtendByOne(AdaptiveInstance& m3, const Graph& G, GRBEnv* env, bool mip_subroutine=true);
    void ComputeAllObjectives(const Graph& G, AdaptiveInstance* m3);
    // Getters.
    int policies() const {return policies_;}
    double worst_case_objective() const {return worst_case_objective_;}
    vector<vector<int>> partition() const {return partition_;}
    vector<Policy> solution() const {return solution_;}
    long most_recent_solution_time() const {return most_recent_solution_time_;}
    // Setters.
    void set_policies(int policies){policies_=policies;}
    void set_scenarios(int scenarios){scenarios_=scenarios;}
    void set_worst_case_objective(int worst_case_objective){worst_case_objective_=worst_case_objective;}
    void set_most_recent_solution_time(long most_recent_solution_time) {most_recent_solution_time_=most_recent_solution_time;}
    void set_solution_policy(int index, Policy policy) {solution_[index]=policy;}
    void set_partition(const vector<vector<int>>& partition) {partition_=partition;}
    void set_unbounded(bool unbounded) {unbounded_=unbounded;}
    void add_to_partition(int index, int scenario) {partition_[index].push_back(scenario);}
private:
    int policies_, scenarios_;
    double worst_case_objective_;
    double average_objective_;
    vector<Policy> solution_;
    vector<vector<int>> partition_;
    long most_recent_solution_time_;
    bool unbounded_;
};

class RobustAlgoModel
{
    // Special Purpose Model for the Enumerative Algorithm - Static Robust Dual Reformulation
    // Can also be used to solve non-robust shortest path interdiction problems - just call update with only
    // the desired scenario.
public:
    // This data will stay throughout the algorithm (all decision variables and non-negativity constraints):
    int scenarios, budget, nodes, arcs; 

    GRBEnv *env_;
    GRBModel *algo_model;

    GRBVar z;
    vector<vector<GRBVar> > pi; // decision variable (for every q, i)
    vector<vector<GRBVar> > lambda; // decision variable (for every q, a)
    vector<GRBVar> x; // interdiction decision variable (for every a)

    // constraints will be added in update based on the subset of [p] that we want to include
    // hold on to linexpr for all constraints
    vector<GRBTempConstr> z_constraints; 
    vector<vector<GRBTempConstr> > dual_constraints;

    RobustAlgoModel() : 
        scenarios(0), budget(0), nodes(0), arcs(0) {};
    RobustAlgoModel(AdaptiveInstance& m3, GRBEnv* env) : 
        scenarios(m3.scenarios()), budget(m3.budget()), nodes(m3.nodes()), arcs(m3.arcs()), env_(env) {};

    void configureModel(const Graph& G, AdaptiveInstance& m3);
    void update(vector<int>& subset);
    void reverse_update(vector<int>& subset);
    Policy Solve();
};


class SetPartitioningModel
{
    // Set Partitioning MIP for solving an adaptive instance.
public:
    SetPartitioningModel() : big_m_(0), scenarios_(0), policies_(0), budget_(0), nodes_(0), arcs_(0) {};
    SetPartitioningModel(int big_m, AdaptiveInstance& instance, GRBEnv* env) : 
        big_m_(big_m), scenarios_(instance.scenarios()), policies_(instance.policies()), budget_(instance.budget()), nodes_(instance.nodes()), arcs_(instance.arcs()), env_(env) {};
    void ConfigureSolver(const Graph& G, AdaptiveInstance& instance);
    void Solve();
    AdaptiveSolution current_solution() const {return current_solution_;}
private:
    const int big_m_;
    int scenarios_, policies_, budget_, nodes_, arcs_;
    GRBEnv *env_;
    GRBModel *sp_model_;
    GRBVar z_var_; // Piecewise linearization dummy variable z.
    vector<vector<GRBVar> > h_var_; // Set partitioning variables h[w][q]==1 if q in subset w.
    vector<vector<GRBVar> > x_var_; // Interdiction variable (for every w, a).
    vector<vector<vector<GRBVar> > > pi_var_; // Dual decision variable (for every w, q, i) pi.
    vector<vector<vector<GRBVar> > > lambda_var_; // Dual decision variable (for every w, q, a) lambda.
    AdaptiveSolution current_solution_; // Solution updated whenever ::Solve() is called.
};

class BendersCallback : public GRBCallback
{
public:
    double upper_bound_, lower_bound_, epsilon_, temp_bound_;
    vector<vector<GRBVar> > h_var_; // Decision variable for every (w, q), set partitioning variable.
    vector<vector<GRBVar> > x_var_; // Decision variable for every (w, a), interdiction policies.
    GRBEnv* env_;
    vector<vector<GRBModel*>> submodels_; // Submodel for every (w, q) policy and scenario.
    vector<vector<vector<GRBVar>>> y_var_; // Shortest path decision variable for every (w, q, a).
    BendersCallback() : scenarios_(0), policies_(0), budget_(0), nodes_(0), arcs_(0) {};
    BendersCallback(AdaptiveInstance& instance, vector<vector<GRBVar> >& h_var, vector<vector<GRBVar> >& x_var, GRBEnv* env) : scenarios_(instance.scenarios()), policies_(instance.policies()), budget_(instance.budget()), nodes_(instance.nodes()), arcs_(instance.arcs()), h_var_(h_var), x_var_(x_var), interdiction_delta_(), env_(env) {cout << "BendersCallback::BendersCallback" << endl;};
    void ConfigureSubModels(const Graph& G, AdaptiveInstance& instance);
    void SolveSubModels();
protected:
    void callback();
private:
    AdaptiveSolution current_solution_;
    int nodes_, arcs_, budget_, scenarios_, policies_;
    int interdiction_delta_;
    void ConfigureIndividualSubModel(const Graph& G, AdaptiveInstance& instance, int w, int q);
    void UpdateSubModels();
};

class SetPartitioningBenders {
public:
    SetPartitioningBenders() : big_m_(0), scenarios_(0), policies_(0), budget_(0), nodes_(0), arcs_(0) {};
    SetPartitioningBenders (int big_m, AdaptiveInstance& instance, GRBEnv* env) : 
        big_m_(big_m), scenarios_(instance.scenarios()), policies_(instance.policies()), budget_(instance.budget()), nodes_(instance.nodes()), arcs_(instance.arcs()), interdiction_delta_(instance.interdiction_delta()), env_(env) {cout << "SetPartitioningBenders::SetPartitioningBenders" << endl;};
    void ConfigureSolver(const Graph& G, AdaptiveInstance& instance);
    void Solve();
private:
    const int big_m_;
    int nodes_, arcs_, budget_, scenarios_, policies_, interdiction_delta_;
    GRBEnv* env_;
    GRBModel* benders_model_;
    GRBVar z_var_; // Decision variable - objective function.
    vector<vector<GRBVar> > h_var_; // Decision variable for every (w, q), set partitioning variable.
    vector<vector<GRBVar> > x_var_; // Decision variable for every (w, a), interdiction policies.
    BendersCallback callback_;
};

// -------- OLD BENDERS ----------
// class BendersSub
// {
// public:
//     int n;
//     int m;
//     int p;
//     vector<GRBEnv *> Subenvs;
//     vector<GRBModel *> Submodels;
// 
//     vector<vector<int> > c_bar; // this is the current objective function cost vector
//     // i.e. - the objective function is c_bar \cdot y
//     // computed during solution tree based on graph costs (c^q and d) and the current
//     // x_bar from the master problem
// 
//     vector<vector<int> > c; // base costs
//     vector<int> d;              // interdiction costs
// 
//     GRBConstr *obj_constr;              // array of constraints for the objective lower bounding constraints over the qs
//                                         // need this as an array to update it
//     vector<GRBVar> zeta_subs;      // dummy objective function variable because we have to argmin over q
//     vector<vector<GRBVar> > y; // main decision variable - arc path selection/flow (one y vector for every q)
//     vector<GRBVar> y_dummy;        // just to construct and push_back y
//     vector<float> y_dummy2;        // just to construct and push_back yhat
//     GRBLinExpr linexpr;                 // when adding the flow constraints, we use this one for outgoing arcs
//     GRBLinExpr linexpr2;                // when adding the flow constraints, we use this one for incoming arcs
//     int rhs;                            // use this also for generating flow constraints
//     string varname;
// 
//     BendersSub() : n(0), m(0), p(0) {};
//     BendersSub(AdaptiveInstance *adaptive_instance);
//     vector<vector<float> > solve(int counter); // now returns a vector of vectors of size p+1, where the first is a singleton with the obj value
//     void update(vector<int> &xhat);
// };
// 
// class BendersSeparation : public GRBCallback
// {
// public:
//     int n;
//     int m;
//     int p;
//     int counter = 0;
//     int cut_count = 0;
// 
//     vector<vector<int> > c; // base costs
//     vector<int> d;              // interdiction costs
// 
//     // a hat vector is numbers, a bar vector is GRBVars
//     GRBLinExpr new_cut;                   // linexpr object for new cut to add to master formulation
//     GRBVar zetabar;                       // 'connecting' GRBVar for zeta
//     vector<GRBVar> xbar;             // 'connecting' GRBVars for x
//     vector<int> xhat;                // current xhat to solve with, i.e. interdiction policy we are subject to
//     vector<float> xprime;            // current best interdiction policy (includes extra x[0] for obj)
//     vector<vector<float> > yhat; // yhat from subproblem, i.e. shortest path given xhat policy (includes extra y[q][0] for objective for each q), it is of size p (vectors), first element of each flow is the objective
// 
//     BendersSub subproblem;
// 
//     // upper and lower bounds and epsilon
//     float zeta_u = GRB_INFINITY;
//     float zeta_l = -GRB_INFINITY;
//     float zeta_temp;
//     float epsilon = 0.000001;
// 
//     BendersSeparation();
//     BendersSeparation(GRBVar &the_zetabar, vector<GRBVar> &the_xbar, AdaptiveInstance *adaptive_instance);
// 
// protected:
//     void callback();
//     void printSep();
// };
// 
// class SetPartitioningBenders
// {
// public:
//     int s = 0;
//     int n;
//     int m;
//     int p;
//     int r_0;
//     string instance_name;
//     string setname;
//     string modelname;
// 
//     float running_time;
//     float optimality_gap;
//     int cut_count;
// 
//     AdaptiveInstance* adaptive_instance_;
//     BendersSeparation sep;
// 
//     GRBEnv *M2Bendersenv;
//     GRBModel *M2Bendersmodel;
// 
//     vector<GRBVar> x; // decision variable - shortest path solution
//     GRBVar zeta;           // objective function
//     GRBLinExpr linexpr;
// 
//     SetPartitioningBenders();
//     SetPartitioningBenders(AdaptiveInstance *adaptive_instance);
//     vector<float> Solve();
// };

// -------- OLD BENDERS ----------

pair<vector<vector<int> >, vector<Policy> > mergeEnumSols(pair<vector<vector<int> >, vector<Policy> >& sol1, pair<vector<vector<int> >, vector<Policy> >& sol2, int w_index);

AdaptiveSolution enumSolve(AdaptiveInstance& m3, const Graph& G, GRBEnv* env);

double UpdateCurrentObjectiveGivenSolution(AdaptiveSolution* current_solution, AdaptiveInstance* m3, const Graph& G);

vector<Policy> InitializeKPolicies(AdaptiveInstance* m3, const Graph& G, GRBEnv* env);

AdaptiveSolution KMeansHeuristic(AdaptiveInstance* m3, const Graph& G, GRBEnv* env);

// pair<vector<vector<int> >, vector<Policy> > extendByOne(pair<vector<vector<int> >, vector<Policy> >& k_solution, AdaptiveInstance& m3, const Graph& G);

long getCurrentTime();

void printSolution(pair<vector<vector<int> >, vector<Policy> >& sol, string solname="");


#endif
