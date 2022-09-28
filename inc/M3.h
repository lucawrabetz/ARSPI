#ifndef M3_H
#define M3_H

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
#include "/home/luw28/gurobi950/linux64/include/gurobi_c++.h"
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
using std::to_string;
using std::stoi;
using std::ifstream;
using std::ofstream;
using std::round;

struct Arc
{
    // Arc struct for the layer graph (directed Arc).
private:
    const int i;
    const int j;
public:
    Arc() : i(0), j(0) {};
    Arc(int the_i, int the_j) : i(the_i), j(the_j) {};
    int get_i() const {return i;}
    int get_j() const {return j;}
};

class Graph
{
    // Graph class (to be read from Arc list)
private:
    int n, kbar, m;
    const string filename;
    // Arc vectors are contiguous 0-(m-1) vectors.
    // This will also apply for cost vectors and weights on arcs that are held in instance classes.
    vector<Arc> arcs; // WHY: to get i and j given an a, only used in logging of solution
    vector<int> subgraph; // WHY: for some cost stuff - lets go look at that soon
    // Linking vectors: when you need the index of an arc from adjacency information:
    //      arc_index_hash maintains a linked list representation 
    //      of indexes in the corresponding 0-(m-1) contiguous arc vectors.
    vector<vector<int> > arc_index_hash; // WHY:
    vector<vector<int> > adjacency_list; // WHY:
    vector<vector<int> > n_n_adjacency_list; // WHY:
public:
    Graph() : n(0), kbar(0), m(0) {};
    Graph(const string &filename, int the_n, int the_k_0);
    vector<Arc> get_arcs() const {return arcs;}
    int get_n() const {return n;}
    int get_kbar() const {return kbar;}
    int get_m() const {return m;}
    vector<int> get_subgraph() const {return subgraph;}
    vector<vector<int>> get_arc_index_hash() const {return arc_index_hash;}
    vector<vector<int>> get_adjacency_list() const {return adjacency_list;}
    void updateGraph(vector<int>& x_bar, bool rev=false);
    void printGraph(vector<vector<int> > costs, vector<int> interdiction_costs, bool is_costs=false) const;
};

class AdaptiveInstance
{
    // Full Adaptive Instance, cost structure data, graph is maintained and passed separately
public:
    int scenarios, policies, budget, nodes, arcs, kbar;
    vector<int> interdiction_costs;
    vector<vector<int> > arc_costs;
    const string directory;
    const string name;
    // When we copy an AdaptiveInstance but only keep a subset of U, we reset the index
    // of the scenarios we keep. Later, we will want the original indices back. 
    // For this reason, we keep a map of the scenarios indices (0-p-1) to their original
    // indices in the AdaptiveInstance it was copied from.
    // An empty map means that the instance is infact an "original" instance.
    vector<int> scenario_index_map;

    // default
    AdaptiveInstance() : scenarios(0), policies(0), budget(0) {}; 

    // main constructor
    AdaptiveInstance(int p, int k, int r_zero, const Graph &G, const string &directory, const string &name) :
        scenarios(p), policies(k), budget(r_zero), nodes(G.get_n()), arcs(G.get_m()), kbar(G.get_kbar()), directory(directory), name(name) {};

    // change U constructor
    AdaptiveInstance(AdaptiveInstance* m3, vector<int>& keep_scenarios);
    //     scenarios(keep_scenarios.size()), policies(m3->policies), budget(m3->budget), nodes(m3->nodes), arcs(m3->arcs), interdiction_costs(m3->interdiction_costs) 
    // {arc_costs=vector<vector<int> >(scenarios, vector<int>(arcs)); 
    //     for (int q=0; q<scenarios; ++q) {arc_costs[q] = m3->arc_costs[keep_scenarios[q]];}}

    // no mutator for scenarios - functionality reserved for change copy constructor
    void set_policies(int k){policies=k;}
    void set_budget(int r_zero){budget=r_zero;}
    void set_costs(vector<int>& interdiction_costs, vector<vector<int> >& arc_costs) 
    {interdiction_costs=interdiction_costs; arc_costs=arc_costs;}

    void printInstance(const Graph &G) const;
    vector<int> dijkstra(int q, const Graph &G);
    void writeCosts();
    void generateCosts(float interdiction, int a, int b, int dist, vector<int> subgraph);
    void readCosts();
    void initCosts(float interdiction, int a, int b, int dist, const Graph &G, bool gen);
    void applyInterdiction(vector<double>& x_bar, bool rev=false);
    float validatePolicy(vector<double>& x_bar, const Graph& G);
};

struct Policy
{
    // Interdiction Policy with Objective Value
public:
    int size;
    double objective;
    vector<double> binary_policy;

    // // default
    Policy() : size(0), objective(0), binary_policy(vector<double>()) {};
    // just m but no policy 
    Policy(int m) : size(m), objective(0), binary_policy(vector<double>(m)) {};
    // full constructor
    Policy(int m, vector<double>& policy, double value) : size(m), binary_policy(policy), objective(value) {};

    // mutators - accept vector or binary new policy
    void set_size(int m) {size=m;}
    void set_policy(vector<double>& policy) {binary_policy=policy;}
    void set_objective(double value) {objective=value;}
};



 
struct AdaptiveSolution
{
    // Full Solution for Instance
public:
    int policies, scenarios;
    double worst_case_objective;
    double average_objective;
    vector<vector<int> > partition;
    vector<Policy> solutions;
    long most_recent_solution_time;

    // default
    AdaptiveSolution() : policies(0), scenarios(0), partition(vector<vector<int>>(0)), solutions(vector<Policy>(0)){};
    // just k and p
    AdaptiveSolution(int k, int p) : policies(k), scenarios(p), partition(vector<vector<int>>(k)), solutions(vector<Policy>(k)) {};
    // full
    AdaptiveSolution(int k, int p, vector<vector<int>> parts, vector<Policy> sols) : policies(k), scenarios(p), partition(parts), solutions(sols) {};

    void logSolution(const Graph& G, AdaptiveInstance& m3, string title, bool policy=true);
    void mergeEnumSols(AdaptiveSolution sol2, AdaptiveInstance* instance2, int split_index);
    void extendByOne(AdaptiveInstance& m3, const Graph& G, bool mip_subroutine=true);

    void set_policies(int k){policies=k;}
    void set_scenarios(int p){scenarios=p;}
    void set_worst_case_objective(int z){worst_case_objective=z;}
    
    void computeAllObjectives(const Graph& G, AdaptiveInstance& m3);
};


class RobustAlgoModel
{
    // Special Purpose Model for the Enumerative Algorithm - Static Robust Dual Reformulation
public:
    // This data will stay throughout the algorithm (all decision variables and non-negativity constraints):
    int scenarios, budget, nodes, arcs; 

    GRBEnv *algo_env;
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
    RobustAlgoModel(AdaptiveInstance& m3) : 
        scenarios(m3.scenarios), budget(m3.budget), nodes(m3.nodes), arcs(m3.arcs) {};

    void configureModel(const Graph& G, AdaptiveInstance& m3);
    void update(vector<int>& subset);
    void reverse_update(vector<int>& subset);
    vector<double> solve(int counter); // returns solution value at 0, arc interdiction policy at 1-m
};


class SetPartitioningModel
{
    // Linear MIP for solving an adaptive instance
public:
    const int M; // set to 500 before initializing in constructor
    int scenarios, policies, budget, nodes, arcs;

    GRBEnv *m3_env;
    GRBModel *m3_model;
    GRBVar z; // objective dummy
    vector<vector<GRBVar> > h_matrix; // set partitioning variables H[w][q]==1 if q in subset w 
    vector<vector<vector<GRBVar> > > pi;     // decision variable (for every w, q, i)
    vector<vector<vector<GRBVar> > > lambda; // decision variable (for every w, q, a)
    vector<vector<GRBVar> > x;                   // interdiction variable (for every w, a)
    vector<vector<float> > x_prime; // final solution, [0] is objective (for every w, a)
    AdaptiveSolution current_solution; // final solution updated whenever ::solve() is called

    SetPartitioningModel() : M(0), scenarios(0), policies(0), budget(0), nodes(0), arcs(0) {};
    SetPartitioningModel(int M, AdaptiveInstance& m3) : 
        M(M), scenarios(m3.scenarios), policies(m3.policies), budget(m3.budget), nodes(m3.nodes), arcs(m3.arcs) {};

    void configureModel(const Graph& G, AdaptiveInstance& m3);
    void solve();
};

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
//     BendersSub();
//     BendersSub(M2ProblemInstance *the_M2Instance);
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
//     BendersSeparation(GRBVar &the_zetabar, vector<GRBVar> &the_xbar, M2ProblemInstance *the_M2Instance);
// 
// protected:
//     void callback();
//     void printSep();
// };
// 
// class AdaptiveBenders
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
//     M2ProblemInstance *M2Instance;
//     BendersSeparation sep;
// 
//     GRBEnv *M2Bendersenv;
//     GRBModel *M2Bendersmodel;
// 
//     vector<GRBVar> x; // decision variable - shortest path solution
//     GRBVar zeta;           // objective function
//     GRBLinExpr linexpr;
// 
//     M2Benders();
//     M2Benders(M2ProblemInstance *the_M2Instance);
//     vector<float> solve();
// };

pair<vector<vector<int> >, vector<Policy> > mergeEnumSols(pair<vector<vector<int> >, vector<Policy> >& sol1, pair<vector<vector<int> >, vector<Policy> >& sol2, int w_index);

AdaptiveSolution enumSolve(AdaptiveInstance& m3, const Graph& G);

// pair<vector<vector<int> >, vector<Policy> > extendByOne(pair<vector<vector<int> >, vector<Policy> >& k_solution, AdaptiveInstance& m3, const Graph& G);

long getCurrentTime();

void printSolution(pair<vector<vector<int> >, vector<Policy> >& sol, string solname="");


#endif
