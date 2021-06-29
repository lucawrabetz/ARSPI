#pragma once
#include "/Library/gurobi902/mac64/include/gurobi_c++.h"

/* 
    rspiModel - class to store all ARSPI instance - solve() only solves 1 policy MIP
    - x and v are stored once and updated as the subsets are enumerated in main
*/

class M2ModelBilinear
{
public:
    int s = 0;
    int n;
    int m;
    int p;
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

class rspiModel
{
private:
    GRBEnv *env;
    GRBModel *model;

    const int s = 0;
    int n, m, k, t; // nodes, arcs

    std::vector<GRBVar> pi;     // decision variable; post interdiction s-i path
    std::vector<GRBVar> lambda; // decision variable; convex combination of scenario costs
    std::vector<GRBVar> x;      // decision variable; interdiction variable

    double z;                         // objective value to be compared to v
    double v = -1;                    // algorithm current best objective
    std::vector<int> x_bar;           // need to store a current copy of x because the GRB version is unstable
    std::vector<double> c_costs = {}; // original arc costs (all scenarios available)
    std::vector<double> d_costs = {}; // interdiction increment costs
    std::vector<int> q_indices = {};  // starting index for each scenario

    std::vector<int> current_scenarios; // subset of scenarios currently being solved
    int *a;                             // pointer to keep position in data array
    int next;                           // next value to add to current scenarios at next loop

public:
    rspiModel();                            // default constructor
    rspiModel(const std::string &filename); // constructor from file name
    void update();                          // update Gurobi Model with current subsets
    void solve(std::vector<int> &data);     // solve for the subsets in data (void, just update x and v)
};
