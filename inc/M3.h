#pragma once
#include <utility>
#include <time.h>
#include <float.h>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include "/Library/gurobi902/mac64/include/gurobi_c++.h"

using std::random_device;
using std::mt19937;
using std::uniform_int_distribution;
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

struct Arc
{
    // Arc struct for the layer graph (directed Arc)
public:
    int i;
    int j;
    Arc();
    Arc(int the_i, int the_j);
};

class LayerGraph
{
    // Layer Graph class (to be read from Arc list)
public:
    int n, m;
    vector<Arc> arcs;

    // arc_index_hash maintains a linked list representation with indexes 
    // to the corresponding arc in the 0-(m-1) vectors, for example the cost
    // or arc object vectors
    vector<vector<int> > arc_index_hash;
    vector<vector<int> > adjacency_list;
    vector<vector<int> > n_n_adjacency_list;

    LayerGraph();
    LayerGraph(const string &filename, int the_n);
    void updateGraph(vector<int>& x_bar, bool rev=false);
    void printGraph(vector<vector<int> > costs, vector<int> interdiction_costs, bool is_costs=false) const;
};

class AdaptiveInstance
{
    // Full Adaptive Instance, cost structure data, graph is maintained and passed separately
public:
    int scenarios, policies, budget, nodes, arcs;
    vector<int> interdiction_costs;
    vector<vector<int> > arc_costs;

    // default
    AdaptiveInstance() : scenarios(0), policies(0), budget(0) {}; 

    // main constructor
    AdaptiveInstance(int p, int k, int r_zero, const LayerGraph &G) :
        scenarios(p), policies(k), budget(r_zero), nodes(G.n), arcs(G.m) {};

    // change U constructor
    AdaptiveInstance(const AdaptiveInstance &m3, vector<int>& keep_scenarios) :
        scenarios(keep_scenarios.size()), policies(m3.policies), budget(m3.budget), nodes(m3.nodes), arcs(m3.arcs), interdiction_costs(m3.interdiction_costs) 
    {arc_costs=vector<vector<int> >(scenarios, vector<int>(arcs)); 
        for (int q : keep_scenarios) {arc_costs[q] = m3.arc_costs[q];}}

    // no mutator for scenarios - functionality reserved for change copy constructor
    void set_policies(int k){policies=k;}
    void set_budget(int r_zero){budget=r_zero;}
    void set_costs(vector<int>& interdiction_costs, vector<vector<int> >& arc_costs) 
    {interdiction_costs=interdiction_costs; arc_costs=arc_costs;}

    void printInstance(const LayerGraph&G) const;
    vector<int> dijkstra(int q, const LayerGraph &G);
    void generateCosts(int interdiction, int min, int max);
    void applyInterdiction(vector<float>& x_bar, bool rev=false);
    float validatePolicy(vector<float>& x_bar, const LayerGraph& G);
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

    void configureModel(const LayerGraph& G, AdaptiveInstance& m3);
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
    vector<vector<GRBVar> > h_matrix; // set partitioning variables H[q][w]==1 if q in subset k 
    vector<vector<vector<GRBVar> > > pi;     // decision variable (for every w, q, i)
    vector<vector<vector<GRBVar> > > lambda; // decision variable (for every w, q, a)
    vector<vector<GRBVar> > x;                   // interdiction variable (for every w, a)
    vector<vector<float> > x_prime; // final solution, [0] is objective (for every w, a)

    SetPartitioningModel() : M(0), scenarios(0), policies(0), budget(0), nodes(0), arcs(0) {};
    SetPartitioningModel(int M, AdaptiveInstance& m3) : 
        M(M), scenarios(m3.scenarios), policies(m3.policies), budget(m3.budget), nodes(m3.nodes), arcs(m3.arcs) {};

    void configureModel(const LayerGraph& G, AdaptiveInstance& m3);
    vector<vector<float> > solve();
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

pair<vector<vector<int> >, vector<vector<double> > > enumSolve(AdaptiveInstance& m3, const LayerGraph& G);

vector<vector<double> > extendByOne(pair<vector<vector<int> >, vector<vector<double> > >& k_solution, AdaptiveInstance& m3);
