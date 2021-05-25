#pragma once
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include "/Library/gurobi902/mac64/include/gurobi_c++.h"

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

    GRBLinExpr linexpr;

    BendersSub();
    BendersSub(M2ProblemInstance *the_M2Instance);
};

class BendersSeparation : public GRBCallback
{
public:
    int n;
    int m;
    int l;
    GRBLinExpr new_cut;        // linexpr object for new cut to add to master formulation
    std::vector<int> xhat;     // current xhat to solve with, i.e. interdiction policy we are subject to
    std::vector<float> xprime; // current best interdiction policy (includes extra x[0] for obj)
    std::vector<float> yhat;   // yhat from subproblem, i.e. shortest path given xhat policy (includes extra y[0] for objective)

    BendersSub subproblem;

    // upper and lower bounds and epsilon
    float zeta_u = 1000000;
    float zeta_l = -1000000;
    float zeta_temp;
    float epsilon = 0.001;

    BendersSeparation();
    BendersSeparation(M2ProblemInstance *the_M2Instance);

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
