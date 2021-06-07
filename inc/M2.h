#pragma once
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include "/Library/gurobi902/mac64/include/gurobi_c++.h"

using std::cout;
using std::string;
using std::to_string;

struct Arc
{
    // little Arc struct for the layer graph (directed Arc)
public:
    int i;
    int j;
    Arc();
    Arc(int the_i, int the_j);
};

class LayerGraph
{
    //layerGraph class (to be read from Arc list)
public:
    int n, m;
    std::vector<Arc> arcs;
    std::vector<std::vector<int>> adjacency_list;
    std::vector<std::vector<int>> reverse_list;

    LayerGraph();
    LayerGraph(const std::string &filename, int the_n);
    void printGraph();
};

class M2ProblemInstance
{
public:
    int n;
    int m;
    int l;
    int r_0;

    std::vector<std::vector<int>> arc_costs;
    std::vector<int> interdiction_costs;

    LayerGraph G;

    M2ProblemInstance();
    M2ProblemInstance(const LayerGraph &the_G, int min, int max, int the_l, int the_r0); // 'normal' constructor
};

class M2ModelBilinear
{
public:
    int s = 0;
    int n;
    int m;
    int l;
    int r_0;
    float running_time;

    M2ProblemInstance *M2Instance;

    GRBEnv *M2env;
    GRBModel *M2model;

    GRBLinExpr linexpr;
    GRBQuadExpr quadexpr;

    // std::vector<std::vector<int>> arc_costs;
    // std::vector<int> interdiction_costs;

    // LayerGraph G;

    std::vector<GRBVar> pi;     // decision variable; post interdiction s-i path
    std::vector<GRBVar> lambda; // decision variable; convex combination of scenario costs
    std::vector<GRBVar> x;      // decision variable; interdiction variable

    M2ModelBilinear(M2ProblemInstance *the_M2Instance); // 'normal' constructor
    float solve();
};

class M2ModelLinear
{
public:
    int s = 0;
    int n;
    int m;
    int l;
    int r_0;
    float running_time;

    M2ProblemInstance *M2Instance;

    GRBEnv *M2env;
    GRBModel *M2model;

    GRBLinExpr linexpr;

    // std::vector<std::vector<int>> arc_costs;
    // std::vector<int> interdiction_costs;

    // LayerGraph G;

    std::vector<std::vector<GRBVar>> pi; // decision variable; post interdiction s-i path for each q
    GRBVar z;                            // decision variable; objective func dummy
    std::vector<GRBVar> x;               // decision variable; interdiction variable

    M2ModelLinear(M2ProblemInstance *the_M2Instance);

    float solve();
};

class BendersSub
{
public:
    int n;
    int m;
    int l;
    GRBEnv *Subenv;
    GRBModel *Submodel;

    std::vector<std::vector<int>> c_bar; // this is the current objective function cost vector
    // i.e. - the objective function is c_bar \cdot y
    // computed during solution tree based on graph costs (c^q and d) and the current
    // x_bar from the master problem

    std::vector<std::vector<int>> c; // base costs
    std::vector<int> d;              // interdiction costs

    GRBConstr *obj_constr;              // array of constraints for the objective lower bounding constraints over the qs
                                        // need this as an array to update it
    GRBVar zeta_sub;                    // dummy objective function variable because we have to argmin over q
    std::vector<std::vector<GRBVar>> y; // main decision variable - arc path selection/flow (one y vector for every q)
    std::vector<GRBVar> y_dummy;        // just to construct and push_back y
    std::vector<float> y_dummy2;        // just to construct and push_back yhat
    GRBLinExpr linexpr;                 // when adding the flow constraints, we use this one for outgoing arcs
    GRBLinExpr linexpr2;                // when adding the flow constraints, we use this one for incoming arcs
    int rhs;                            // use this also for generating flow constraints
    string varname;

    BendersSub();
    BendersSub(M2ProblemInstance *the_M2Instance);
    std::vector<std::vector<float>> solve(int counter); // now returns a vector of vectors of size l+1, where the first is a singleton with the obj value
    void update(std::vector<int> &xhat);
};

class BendersSeparation : public GRBCallback
{
public:
    int n;
    int m;
    int l;
    int counter = 0;

    std::vector<std::vector<int>> c; // base costs
    std::vector<int> d;              // interdiction costs

    // a hat vector is numbers, a bar vector is GRBVars
    GRBLinExpr new_cut;                   // linexpr object for new cut to add to master formulation
    GRBVar zetabar;                       // 'connecting' GRBVar for zeta
    std::vector<GRBVar> xbar;             // 'connecting' GRBVars for x
    std::vector<int> xhat;                // current xhat to solve with, i.e. interdiction policy we are subject to
    std::vector<float> xprime;            // current best interdiction policy (includes extra x[0] for obj)
    std::vector<std::vector<float>> yhat; // yhat from subproblem, i.e. shortest path given xhat policy (includes extra y[0][0] for objective), it is of size l+1, first is singleton obj, next l are the y vectors for every q

    BendersSub subproblem;

    // upper and lower bounds and epsilon
    float zeta_u = GRB_INFINITY;
    float zeta_l = -GRB_INFINITY;
    float zeta_temp;
    float epsilon = 0.0001;

    BendersSeparation();
    BendersSeparation(GRBVar &the_zetabar, std::vector<GRBVar> &the_xbar, M2ProblemInstance *the_M2Instance);

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
    int l;
    int r_0;
    float running_time;

    M2ProblemInstance *M2Instance;
    BendersSeparation sep;

    GRBEnv *M2Bendersenv;
    GRBModel *M2Bendersmodel;

    std::vector<GRBVar> x; // decision variable - shortest path solution
    GRBVar zeta;           // objective function
    GRBLinExpr linexpr;

    M2Benders(M2ProblemInstance *the_M2Instance);
    std::vector<float> solve();
};
