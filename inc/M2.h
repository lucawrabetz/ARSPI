#pragma once
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include "/Library/gurobi902/mac64/include/gurobi_c++.h"

using std::stringstream;
using std::endl;
using std::cout;
using std::string;
using std::vector;
using std::to_string;
using std::stoi;
using std::ifstream;


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
    vector<vector<int>> adjacency_list;
    vector<vector<int>> reverse_list;

    LayerGraph();
    LayerGraph(const string &filename, int the_n);
    void printGraph();
};

class M2ProblemInstance
{
    // Full Instance of an M2 problem (a Layer Graph + arc costs and interdiction costs)
    // TODO: randomized cost generation in this class
public:
    int n;
    int m;
    int p;
    int r_0;

    vector<vector<int>> arc_costs;
    vector<int> interdiction_costs;

    LayerGraph G;

    M2ProblemInstance();
    M2ProblemInstance(const LayerGraph &the_G, int min, int max, int the_p, int the_r0);
};


class M2ModelLinear
{
    // Linear MIP for solving an M2ProblemInstance
public:
    int s = 0;
    int n;
    int m;
    int p;
    int r_0;
    float running_time;
    float optimality_gap;

    M2ProblemInstance *M2Instance;

    GRBEnv *M2env;
    GRBModel *M2model;

    GRBLinExpr linexpr;

    vector<vector<GRBVar>> pi;     // decision variable;
    vector<vector<GRBVar>> lambda; // decision variable;
    GRBVar z;                                // decision variable; objective func dummy
    vector<GRBVar> x;                   // decision variable; interdiction variable
    vector<float> x_prime; // float vector for x_final, [0] is objective

    M2ModelLinear();
    M2ModelLinear(M2ProblemInstance *the_M2Instance);

    vector<float> solve();
};

class BendersSub
{
public:
    int n;
    int m;
    int p;
    vector<GRBEnv *> Subenvs;
    vector<GRBModel *> Submodels;

    vector<vector<int>> c_bar; // this is the current objective function cost vector
    // i.e. - the objective function is c_bar \cdot y
    // computed during solution tree based on graph costs (c^q and d) and the current
    // x_bar from the master problem

    vector<vector<int>> c; // base costs
    vector<int> d;              // interdiction costs

    GRBConstr *obj_constr;              // array of constraints for the objective lower bounding constraints over the qs
                                        // need this as an array to update it
    vector<GRBVar> zeta_subs;      // dummy objective function variable because we have to argmin over q
    vector<vector<GRBVar>> y; // main decision variable - arc path selection/flow (one y vector for every q)
    vector<GRBVar> y_dummy;        // just to construct and push_back y
    vector<float> y_dummy2;        // just to construct and push_back yhat
    GRBLinExpr linexpr;                 // when adding the flow constraints, we use this one for outgoing arcs
    GRBLinExpr linexpr2;                // when adding the flow constraints, we use this one for incoming arcs
    int rhs;                            // use this also for generating flow constraints
    string varname;

    BendersSub();
    BendersSub(M2ProblemInstance *the_M2Instance);
    vector<vector<float>> solve(int counter); // now returns a vector of vectors of size p+1, where the first is a singleton with the obj value
    void update(vector<int> &xhat);
};

class BendersSeparation : public GRBCallback
{
public:
    int n;
    int m;
    int p;
    int counter = 0;
    int cut_count = 0;

    vector<vector<int>> c; // base costs
    vector<int> d;              // interdiction costs

    // a hat vector is numbers, a bar vector is GRBVars
    GRBLinExpr new_cut;                   // linexpr object for new cut to add to master formulation
    GRBVar zetabar;                       // 'connecting' GRBVar for zeta
    vector<GRBVar> xbar;             // 'connecting' GRBVars for x
    vector<int> xhat;                // current xhat to solve with, i.e. interdiction policy we are subject to
    vector<float> xprime;            // current best interdiction policy (includes extra x[0] for obj)
    vector<vector<float>> yhat; // yhat from subproblem, i.e. shortest path given xhat policy (includes extra y[q][0] for objective for each q), it is of size p (vectors), first element of each flow is the objective

    BendersSub subproblem;

    // upper and lower bounds and epsilon
    float zeta_u = GRB_INFINITY;
    float zeta_l = -GRB_INFINITY;
    float zeta_temp;
    float epsilon = 0.000001;

    BendersSeparation();
    BendersSeparation(GRBVar &the_zetabar, vector<GRBVar> &the_xbar, M2ProblemInstance *the_M2Instance);

protected:
    void callback();
    void printSep();
};

class M2Benders
{
public:
    int s = 0;
    int n;
    int m;
    int p;
    int r_0;
    float running_time;
    float optimality_gap;
    int cut_count;

    M2ProblemInstance *M2Instance;
    BendersSeparation sep;

    GRBEnv *M2Bendersenv;
    GRBModel *M2Bendersmodel;

    vector<GRBVar> x; // decision variable - shortest path solution
    GRBVar zeta;           // objective function
    GRBLinExpr linexpr;

    M2Benders();
    M2Benders(M2ProblemInstance *the_M2Instance);
    vector<float> solve();
};
