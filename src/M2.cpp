#include "../inc/M2.h"

Arc::Arc()
{
    i = 0;
    j = 0;
}

Arc::Arc(int the_i, int the_j)
{
    i = the_i;
    j = the_j;
}

LayerGraph::LayerGraph()
{
    n = 0;
    m = 0;
}

LayerGraph::LayerGraph(const std::string &filename, int the_n)
{
    std::string line;
    std::ifstream myfile(filename);

    if (myfile.is_open())
    {

        n = the_n;
        const char *cline;
        int i;
        int j;
        std::vector<int> new_vector;

        for (int i = 0; i < n; i++)
        {
            new_vector = {};
            adjacency_list.push_back(new_vector);
            reverse_list.push_back(new_vector);
        }

        while (getline(myfile, line))
        {
            // read line
            // assign i and j, create arc
            // add to arc set
            // add to adjacency list
            cline = line.c_str();
            sscanf(cline, "%d %d", &i, &j);
            arcs.push_back(Arc(i, j));
            adjacency_list[i].push_back(j);
            reverse_list[j].push_back(i);
        }
        m = arcs.size();
    }
}

void LayerGraph::printGraph()
{
    std::cout << "n: " << n << ", m: " << m << "\n";
    for (int a = 0; a < m; a++)
    {
        std::cout << "(" << arcs[a].i << "," << arcs[a].j << ")\n";
    }
}

M2ProblemInstance::M2ProblemInstance()
{
    r_0 = 0;
}

M2ProblemInstance::M2ProblemInstance(const LayerGraph &the_G, int min, int max, int the_l, int the_r0)
{
    // ------ Assign graph and random costs ------
    // ------ Variables and int parameters ------
    G = the_G;
    n = G.n;
    m = G.m;
    l = the_l;
    r_0 = the_r0;

    for (int a = 0; a < m; a++)
    {
        // interdiction_costs.push_back((max - min));
        // for simplegraph.txt
        interdiction_costs.push_back(100);
    }

    // std::random_device rd;                           // obtain a random number from hardware
    // std::mt19937 gen(rd());                          // seed the generator
    // std::uniform_int_distribution<> distr(min, max); // define the range

    // for (int q = 0; q < l; q++)
    // {
    //     std::vector<int> new_vector = {};
    //     arc_costs.push_back(new_vector);
    //     for (int a = 0; a < m; a++)
    //     {
    //         arc_costs[q].push_back(distr(gen)); // assign arc cost between min and max
    //     }
    // }

    // hardcoded example "simplegraph.txt"
    std::vector<int> costs1 = {9, 12, 3};
    std::vector<int> costs2 = {9, 1, 10};
    std::vector<int> costs3 = {9, 12, 10};
    arc_costs.push_back(costs1);
    arc_costs.push_back(costs2);
    arc_costs.push_back(costs3);
}

// ------ 'Straightup' MIP Formulations for M2 ------

M2ModelBilinear::M2ModelBilinear(M2ProblemInstance *the_M2Instance)
{
    try
    {
        // ------ Assign Instance ------
        M2Instance = the_M2Instance;

        // // ------ Assign graph and random costs ------
        // // ------ Variables and int parameters ------
        n = M2Instance->G.n;
        m = M2Instance->G.m;
        l = M2Instance->l;
        r_0 = M2Instance->r_0;

        // for (int a = 0; a < m; a++)
        // {
        //     // interdiction_costs.push_back((max - min));
        //     // for simplegraph.txt
        //     interdiction_costs.push_back(100);
        // }

        // // std::random_device rd;                           // obtain a random number from hardware
        // // std::mt19937 gen(rd());                          // seed the generator
        // // std::uniform_int_distribution<> distr(min, max); // define the range

        // // for (int q = 0; q < l; q++)
        // // {
        // //     std::vector<int> new_vector = {};
        // //     arc_costs.push_back(new_vector);
        // //     for (int a = 0; a < m; a++)
        // //     {
        // //         arc_costs[q].push_back(distr(gen)); // assign arc cost between min and max
        // //     }
        // // }

        // // hardcoded example "simplegraph.txt"
        // std::vector<int> costs1 = {3, 4, 3, 4, 3, 4, 3};
        // std::vector<int> costs2 = {3, 1, 3, 1, 3, 1, 10};
        // std::vector<int> costs3 = {3, 4, 3, 4, 3, 4, 10};
        // arc_costs.push_back(costs1);
        // arc_costs.push_back(costs2);
        // arc_costs.push_back(costs3);

        // ------ Initialize model and environment ------
        M2env = new GRBEnv();
        M2model = new GRBModel(*M2env);

        // ------ Decision variables ------
        std::string varname;

        // interdiction policy on arcs 'x'
        for (int a = 0; a < m; a++)
        {
            varname = "x_" + std::to_string(a);
            x.push_back(M2model->addVar(0, 1, 0, GRB_BINARY, varname));
        }

        // convex combination lambda on scenarios
        for (int q = 0; q < l; q++)
        {
            varname = "lambda_" + std::to_string(q);
            lambda.push_back(M2model->addVar(0, 1, 0, GRB_CONTINUOUS, varname));
        }

        // post interdiction shortest path (s-i) 'pi'
        for (int i = 0; i < n; i++)
        {
            varname = "pi_" + std::to_string(i);
            pi.push_back(M2model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
        }

        // ------ Constraints ------
        // convex combination (sum of lambdas = 1)
        linexpr = 0;
        for (int q = 0; q < l; q++)
        {
            linexpr += lambda[q];
        }
        M2model->addConstr(linexpr == 1);

        // budget constraint
        linexpr = 0;
        for (int a = 0; a < m; a++)
        {
            linexpr += x[a];
        }
        M2model->addConstr(linexpr <= r_0);

        // main constraint for each arc
        for (int a = 0; a < m; a++)
        {
            quadexpr = 0;
            for (int q = 0; q < l; q++)
            {
                quadexpr += (M2Instance->arc_costs[q][a] + M2Instance->interdiction_costs[a] * x[a]) * lambda[q];
            }
            M2model->addQConstr((pi[M2Instance->G.arcs[a].j] - pi[M2Instance->G.arcs[a].i]) <= quadexpr);
        }

        // pi[0] = 0
        M2model->addConstr(pi[0] == 0);
        linexpr = 0;
        linexpr += pi[n - 1];

        M2model->setObjective(linexpr, GRB_MAXIMIZE);
        M2model->update();
    }
    catch (GRBException e)
    {
        std::cout << "Gurobi error number [M2Model, constructor]: " << e.getErrorCode() << "\n";
        std::cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        std::cout << "Non-gurobi error during optimization [M2Model]"
                  << "\n";
    }
}

float M2ModelBilinear::solve()
{
    // try
    // {
    //     // ------ Initialize model and environment ------
    //     M2env = new GRBEnv();
    //     M2model = new GRBModel(*M2env);

    //     // ------ Decision variables ------
    //     std::string varname;

    //     // interdiction policy on arcs
    //     for (int a = 0; a < m; a++)
    //     {
    //         varname = "x_" + std::to_string(a);
    //         x.push_back(M2model->addVar(0, 1, 0, GRB_BINARY, varname));
    //     }

    //     // convex combination lambda on scenarios
    //     for (int q = 0; q < l; q++)
    //     {
    //         varname = "lambda_" + std::to_string(q);
    //         lambda.push_back(M2model->addVar(0, 1, 0, GRB_CONTINUOUS, varname));
    //     }

    //     // post interdiction shortest path (s-i)
    //     for (int i = 0; i < n; i++)
    //     {
    //         varname = "pi_" + std::to_string(i);
    //         pi.push_back(M2model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
    //     }

    //     // ------ Constraints ------
    //     // convex combination (sum of lambdas = 1)
    //     linexpr = 0;
    //     for (int q = 0; q < l; q++)
    //     {
    //         linexpr += lambda[q];
    //     }
    //     M2model->addConstr(linexpr == 1);

    //     // budget constraint
    //     linexpr = 0;
    //     for (int a = 0; a < m; a++)
    //     {
    //         linexpr += x[a];
    //     }
    //     M2model->addConstr(linexpr <= r_0);

    //     // convex combination (sum of lambdas = 1)
    //     linexpr = 0;
    //     for (int q = 0; q < l; q++)
    //     {
    //         linexpr += lambda[q];
    //     }
    //     M2model->addConstr(linexpr <= r_0);

    //     // main constraint for each arc
    //     for (int a = 0; a < m; a++)
    //     {
    //         quadexpr = 0;
    //         for (int q = 0; q < l; q++)
    //         {
    //             quadexpr += (arc_costs[q][a] + interdiction_costs[a] * x[a]) * lambda[q];
    //         }
    //         M2model->addQConstr((pi[G.arcs[a].j] - pi[G.arcs[a].i]) <= quadexpr);
    //     }

    //     // pi[0] = 0
    //     M2model->addConstr(pi[0] == 0);
    //     linexpr = 0;
    //     linexpr += pi[n - 1];

    //     M2model->setObjective(linexpr, GRB_MAXIMIZE);
    //     M2model->update();
    // }
    // catch (GRBException e)
    // {
    //     std::cout << "Gurobi error number [M2Model, constructor]: " << e.getErrorCode() << "\n";
    //     std::cout << e.getMessage() << "\n";
    // }
    // catch (...)
    // {
    //     std::cout << "Non-gurobi error during optimization [M2Model]"
    //               << "\n";
    // }
    try
    {
        clock_t model_begin = clock();
        M2model->optimize();
        running_time = float(clock() - model_begin) / CLOCKS_PER_SEC;
        std::cout << "Objective: " << M2model->get(GRB_DoubleAttr_ObjVal) << "\n";
        std::cout << "Running time: " << running_time << "\n";
        for (int a = 0; a < m; a++)
        {
            std::cout << "x_" << a << "(" << M2Instance->G.arcs[a].i << "," << M2Instance->G.arcs[a].j << ")"
                      << ": " << x[a].get(GRB_DoubleAttr_X) << "\n";
        }

        return running_time;
    }
    catch (GRBException e)
    {
        std::cout << "Gurobi error number [M2Model, solveMIP]: " << e.getErrorCode() << "\n";
        std::cout << e.getMessage() << "\n";
        return 0;
    }
    catch (...)
    {
        std::cout << "Non-gurobi error during optimization [M2Model]"
                  << "\n";
        return 0;
    }
}

M2ModelLinear::M2ModelLinear(M2ProblemInstance *the_M2Instance)
{
    try
    {
        // ------ Assign Instance ------
        M2Instance = the_M2Instance;

        // // ------ Assign graph and random costs ------
        // ------ Variables and int parameters ------
        n = M2Instance->G.n;
        m = M2Instance->G.m;
        l = M2Instance->l;
        r_0 = M2Instance->r_0;

        // for (int a = 0; a < m; a++)
        // {
        //     // interdiction_costs.push_back((max - min));
        //     // for simplegraph.txt
        //     interdiction_costs.push_back(100);
        // }

        // // std::random_device rd;                           // obtain a random number from hardware
        // // std::mt19937 gen(rd());                          // seed the generator
        // // std::uniform_int_distribution<> distr(min, max); // define the range

        // // for (int q = 0; q < l; q++)
        // // {
        // //     std::vector<int> new_vector = {};
        // //     arc_costs.push_back(new_vector);
        // //     for (int a = 0; a < m; a++)
        // //     {
        // //         arc_costs[q].push_back(distr(gen)); // assign arc cost between min and max
        // //     }
        // // }

        // // hardcoded example "simplegraph.txt"
        // std::vector<int> costs1 = {3, 4, 3, 4, 3, 4, 3};
        // std::vector<int> costs2 = {3, 1, 3, 1, 3, 1, 10};
        // std::vector<int> costs3 = {3, 4, 3, 4, 3, 4, 10};
        // arc_costs.push_back(costs1);
        // arc_costs.push_back(costs2);
        // arc_costs.push_back(costs3);

        // ------ Initialize model and environment ------
        M2env = new GRBEnv();
        M2model = new GRBModel(*M2env);

        // ------ Decision variables ------
        std::string varname;

        // interdiction policy on arcs 'x'
        for (int a = 0; a < m; a++)
        {
            varname = "x_" + std::to_string(a);
            x.push_back(M2model->addVar(0, 1, 0, GRB_BINARY, varname));
        }

        // objective func dummy 'z'
        varname = "z";
        z = M2model->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS, varname);

        // post interdiction shortest path (s-i) 'pi'
        for (int q = 0; q < l; q++)
        {
            std::vector<GRBVar> new_vector = {};
            pi.push_back(new_vector);
            for (int i = 0; i < n; i++)
            {
                varname = "pi_" + std::to_string(q) + "_" + std::to_string(i);
                pi[q].push_back(M2model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
            }
        }

        // ------ Constraints ------
        // budget constraint
        linexpr = 0;
        for (int a = 0; a < m; a++)
        {
            linexpr += x[a];
        }
        M2model->addConstr(linexpr <= r_0);

        // z constraints
        linexpr = 0;
        for (int q = 0; q < l; q++)
        {
            M2model->addConstr(z <= pi[q][n - 1] - pi[q][s]);
        }

        // main constraint for each arc
        for (int a = 0; a < m; a++)
        {
            for (int q = 0; q < l; q++)
            {
                M2model->addConstr((pi[q][M2Instance->G.arcs[a].j] - pi[q][M2Instance->G.arcs[a].i]) <= M2Instance->arc_costs[q][a] + M2Instance->interdiction_costs[a] * x[a]);
            }
        }

        // pi[0] = 0
        for (int q = 0; q < l; q++)
        {
            M2model->addConstr(pi[q][0] == 0);
        }
        M2model->update();
    }
    catch (GRBException e)
    {
        std::cout << "Gurobi error number [M2Model, constructor]: " << e.getErrorCode() << "\n";
        std::cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        std::cout << "Non-gurobi error during optimization [M2Model]"
                  << "\n";
    }
}

float M2ModelLinear::solve()
{
    try
    {
        clock_t model_begin = clock();
        M2model->optimize();
        running_time = float(clock() - model_begin) / CLOCKS_PER_SEC;
        std::cout << "Objective: " << M2model->get(GRB_DoubleAttr_ObjVal) << "\n";
        std::cout << "Running time: " << running_time << "\n";
        for (int a = 0; a < m; a++)
        {
            std::cout << "x_" << a << "(" << M2Instance->G.arcs[a].i << "," << M2Instance->G.arcs[a].j << ")"
                      << ": " << x[a].get(GRB_DoubleAttr_X) << "\n";
        }

        return running_time;
    }
    catch (GRBException e)
    {
        std::cout << "Gurobi error number [M2Model, solveMIP]: " << e.getErrorCode() << "\n";
        std::cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        std::cout << "Non-gurobi error during optimization [M2Model]"
                  << "\n";
    }
}

// ------ Bender's Schemes for M2 ------
BendersSub::BendersSub()
{
    n = 0;
    m = 0;
    l = 0;
}

BendersSub::BendersSub(M2ProblemInstance *the_M2Instance)
{
    // ------ Initialize Basic Parameters ------
    n = the_M2Instance->n;
    m = the_M2Instance->m;
    l = the_M2Instance->l;
}

BendersSeparation::BendersSeparation()
{
    n = 0;
    m = 0;
    l = 0;
}

BendersSeparation::BendersSeparation(M2ProblemInstance *the_M2Instance)
{
    try
    {
        // ------ Initialize Basic Parameters ------
        n = the_M2Instance->n;
        m = the_M2Instance->m;
        l = the_M2Instance->l;

        // ------ Initialize Submodel ------
        subproblem = BendersSub(the_M2Instance);

        for (int i = 0; i < m; ++i)
        {
            xhat.push_back(0);
            yhat.push_back(0);
            xhat.push_back(0);
        }
    }
    catch (GRBException e)
    {
        std::cout << "Gurobi error number [BendersSeparation, constructor]: " << e.getErrorCode() << "\n";
        std::cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        std::cout << "Non-gurobi error during optimization [BendersSeparation]"
                  << "\n";
    }
}

void BendersSeparation::callback()
{
    if (where == GRB_CB_MIPSOL)
    {
        // xhat = current solution from master problem
        // zeta_u = current objective from master problem
        // BendersSub.update(xhat);
        // yhat = BendersSub.Solve;
        // zeta_temp = yhat[0]; // first element of the yhat vector is the objective
        // use yhat[1-m] to create new cut from LinExpr
        // add lazy cut to main model
        if (zeta_l < zeta_temp)
        {
            // xprime = xhat;
        }
    }
}

M2Benders::M2Benders(M2ProblemInstance *the_M2Instance)
{
    try
    {
        // ------ Assign Instance ------
        M2Instance = the_M2Instance;

        // ------ Initialize model and environment ------
        M2Bendersenv = new GRBEnv();
        M2Bendersmodel = new GRBModel(*M2Bendersenv);

        // ------ Initialize separation/callback object ------
        sep = BendersSeparation(M2Instance);

        // ------ Variables and int parameters ------
        n = M2Instance->G.n;
        m = M2Instance->G.m;
        l = M2Instance->l;
        r_0 = M2Instance->r_0;

        // ------ Decision variables ------
        std::string varname;

        // interdiction variable 'x'
        for (int a = 0; a < m; a++)
        {
            varname = "x_" + std::to_string(a);
            x.push_back(M2Bendersmodel->addVar(0, 1, 0, GRB_BINARY, varname));
        }

        // objective function variable 'zeta'
        varname = "zeta";
        zeta = M2Bendersmodel->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS, varname);

        // ------ No Constraints Initially! ------
    }
    catch (GRBException e)
    {
        std::cout << "Gurobi error number [M2Benders, constructor]: " << e.getErrorCode() << "\n";
        std::cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        std::cout << "Non-gurobi error during optimization [BendersSPSub]"
                  << "\n";
    }
}

std::vector<float> M2Benders::solve()
{
    // ------ Set Callback on Master Model
    M2Bendersmodel->setCallback(&sep);

    // ------ Optimize Inside Benders Scheme -------
    while (sep.zeta_u - sep.zeta_l >= sep.epsilon)
    {
        M2Bendersmodel->optimize();
    }

    delete sep.subproblem.Subenv;
    delete sep.subproblem.Submodel;
    return sep.xprime;
}

// floar BendersSPSub::solve()
// {
// }
