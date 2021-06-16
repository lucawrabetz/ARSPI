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

    // ------ Initialize d and c costs ------
    for (int a = 0; a < m; a++)
    {
        d.push_back(the_M2Instance->interdiction_costs[a]);
    }

    for (int q = 0; q < l; q++)
    {
        c.push_back(the_M2Instance->arc_costs[q]);
        c_bar.push_back(the_M2Instance->arc_costs[q]);
    }

    // ------ Initialize Environment and Model ------
    Subenv = new GRBEnv();
    Submodel = new GRBModel(*Subenv);

    // ------ Decision Variables ------
    varname = "zeta_sub";
    zeta_sub = Submodel->addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS, varname);
    for (int q = 0; q < l; ++q)
    {
        y_dummy = {};
        y.push_back(y_dummy);
        for (int a = 0; a < m; ++a)
        {
            varname = "y_" + to_string(q) + "_" + to_string(a);
            y[q].push_back(Submodel->addVar(0, 1, 0, GRB_CONTINUOUS, varname));
        }
    }

    // ------ Constraints ------
    // constraints to bound objective value over q
    // obj_constr = new GRBConstr[l]
    obj_constr = new GRBConstr[l];
    for (int q = 0; q < l; ++q)
    {
        linexpr = 0;
        for (int a = 0; a < m; ++a)
        {
            linexpr += c_bar[q][a] * y[q][a];
        }
        obj_constr[q] = Submodel->addConstr(zeta_sub >= linexpr);
    }

    // flow constraints
    for (int q = 0; q < l; ++q)
    {
        for (int i = 0; i < n; ++i)
        {
            linexpr = 0;
            if (i == 0)
            {
                rhs = 1;
            }
            else if (i == n - 1)
            {
                rhs = -1;
            }
            else
            {
                rhs = 0;
            }
            for (int a = 0; a < m; a++)
            {
                if (the_M2Instance->G.arcs[a].i == i)
                {
                    // this arc is an outgoing arc for i
                    linexpr += (1) * y[q][a];
                }
                else if (the_M2Instance->G.arcs[a].j == i)
                {
                    // this arc is an outgoing arc for i
                    linexpr += (-1) * y[q][a];
                }
                else
                {
                    // this arc does not include i as an endpoint
                    continue;
                }
            }
            Submodel->addConstr(linexpr == rhs);
        }
    }
    Submodel->update();
}

void BendersSub::update(std::vector<int> &xhat)
{
    cout << "\nsubmodel: updating based on xbar, new interdiction policy: \n";
    for (int a = 0; a < m; ++a)
    {
        cout << "xbar_" << a << ": " << xhat[a] << "\n";
    }

    // update array of constraints instead of only the parameter vector
    for (int q = 0; q < l; ++q)
    {
        for (int a = 0; a < m; ++a)
        {
            // model.chgCoeff(obj_constr[q], variable, new cost)
            if (xhat[a] > 0.5)
            {
                Submodel->chgCoeff(obj_constr[q], y[q][a], -(c[q][a] + d[a]));
                // c_bar[q][a] = c[q][a] + d[a];
            }
            else
            {
                Submodel->chgCoeff(obj_constr[q], y[q][a], -(c[q][a]));
                // c_bar[q][a] = c[q][a];
            }
            // constraint[q].set...
        }
    }
    Submodel->update();
}

std::vector<std::vector<float>> BendersSub::solve(int counter)
{
    try
    {
        std::vector<std::vector<float>> yhat;
        string modelname = "submodel_" + to_string(counter) + ".lp";
        Submodel->write(modelname);
        Submodel->optimize();

        y_dummy2 = {};
        yhat.push_back(y_dummy2);
        yhat[0].push_back(Submodel->get(GRB_DoubleAttr_ObjVal));
        cout << "submodel_obj: " << yhat[0][0];
        cout << "\narc values: \n";

        for (int q = 0; q < l; ++q)
        {
            y_dummy2 = {};
            yhat.push_back(y_dummy2);
            cout << "q = " << q << "\n";
            for (int a = 0; a < m; ++a)
            {
                yhat[q + 1].push_back(y[q][a].get(GRB_DoubleAttr_X));
                cout << "y_" << q << "_" << a << ": " << y[q][a].get(GRB_DoubleAttr_X) << "\n";
            }
        }
        return yhat;
    }
    catch (GRBException e)
    {
        std::cout << "Gurobi error number [BendersSub::solve]: " << e.getErrorCode() << "\n";
        std::cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        std::cout << "Non-gurobi error during optimization [BendersSub::Solve]"
                  << "\n";
    }
}

BendersSeparation::BendersSeparation()
{
    n = 0;
    m = 0;
    l = 0;
}

BendersSeparation::BendersSeparation(GRBVar &the_zetabar, std::vector<GRBVar> &the_xbar, M2ProblemInstance *the_M2Instance)
{
    try
    {
        // ------ Initialize Basic Parameters ------
        n = the_M2Instance->n;
        m = the_M2Instance->m;
        l = the_M2Instance->l;

        // ------ Initialize Submodel ------
        subproblem = BendersSub(the_M2Instance);

        // ------ Initialize d and c costs ------
        for (int a = 0; a < m; a++)
        {
            d.push_back(the_M2Instance->interdiction_costs[a]);
        }

        for (int q = 0; q < l; q++)
        {
            c.push_back(the_M2Instance->arc_costs[q]);
        }

        // ------ Initialize Variable containers ------
        xbar = the_xbar;
        zetabar = the_zetabar;
        xprime.push_back(0);

        std::vector<float> y_dummy = {0};
        yhat.push_back(y_dummy);

        for (int q = 0; q < l; ++q)
        {
            if (q == 0)
            {
                y_dummy = {};
                yhat.push_back(y_dummy);
                for (int a = 0; a < m; ++a)
                {
                    xhat.push_back(0);
                    xprime.push_back(0);
                    yhat[q + 1].push_back(0);
                }
            }
            else
            {
                y_dummy = {};
                yhat.push_back(y_dummy);
                for (int a = 0; a < m; ++a)
                {
                    yhat[q + 1].push_back(0);
                }
            }
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
        if (zeta_u - zeta_l >= epsilon)
        {
            counter++;
            // xhat = current solution from master problem, then update subproblem
            for (int a = 0; a < m; ++a)
            {
                xhat[a] = getSolution(xbar[a]);
            }
            subproblem.update(xhat);
            zeta_u = GRB_CB_MIPSOL_OBJBST; // best obj found so far (entire tree)

            cout << "\n\n\n\nsolving sub from callback: \n\n\n\n";
            yhat = subproblem.solve(counter);
            // cout << "\nsubobjective: " << yhat[0][0] << "\n";

            for (int q = 1; q < l + 1; ++q)
            {
                for (int a = 0; a < m; ++a)
                {
                    // cout << "\nyhat[" << q - 1 << "][" << a << "]: " << yhat[q][a] << "\n";
                }
            }
            zeta_temp = yhat[0][0]; // first element of the yhat vector is the objective
            // use yhat[(q-l)+1][1-m] to create new cut from LinExpr
            for (int q = 0; q < l; ++q)
            {
                new_cut = 0;
                for (int a = 0; a < m; ++a)
                {
                    new_cut += (c[q][a] + d[a] * xbar[a]) * yhat[q + 1][a];
                }
                // add lazy cut to main model
                try
                {
                    addLazy(zetabar <= new_cut);
                    cout << "\nallegedly added cut: "
                         << "zetabar"
                         << "<=" << new_cut << "\n";
                }
                catch (GRBException e)
                {
                    std::cout << "Gurobi error number [BendersSeparation, addLazy]: " << e.getErrorCode() << "\n";
                    std::cout << e.getMessage() << "\n";
                }
            }

            if (zeta_l < zeta_temp)
            {
                zeta_l = zeta_temp;
            }
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

        M2Bendersmodel->getEnv().set(GRB_IntParam_LazyConstraints, 1);

        // ------ Variables and int parameters ------
        n = M2Instance->n;
        m = M2Instance->m;
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

        // ------ Initialize separation/callback object ------
        sep = BendersSeparation(zeta, x, M2Instance);

        // ------ Only Budget Constraints Initially! ------
        linexpr = 0;
        for (int a = 0; a < m; a++)
        {
            linexpr += x[a];
        }
        M2Bendersmodel->addConstr(linexpr <= r_0);

        // ------ Trying to add the first lazy contraint ------
        cout << "\n\n\n\nsolving sub from constructor: \n\n\n\n";
        sep.yhat = sep.subproblem.solve(0);
        for (int q = 0; q < l; ++q)
        {
            linexpr = 0;
            for (int a = 0; a < m; ++a)
            {
                linexpr += (sep.c[q][a] + (sep.d[a] * x[a])) * sep.yhat[q + 1][a];
                // cout << "\nyhat[" << q << "][" << a << "]: " << sep.yhat[q + 1][a] << "\n";
            }
            M2Bendersmodel->addConstr(zeta <= linexpr);
        }
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
    M2Bendersmodel->write("simplegraph1_benders_before.lp");

    try
    {
        M2Bendersmodel->optimize();
    }
    catch (GRBException e)
    {
        std::cout << "Gurobi error number [M2Benders.optimize()]: " << e.getErrorCode() << "\n";
        std::cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        std::cout << "Non-gurobi error during optimization [M2Benders.optimize()]"
                  << "\n";
    }
    M2Bendersmodel->write("simplegraph1_benders_after.lp");
    // while (sep.zeta_u - sep.zeta_l >= sep.epsilon)
    // {
    //     cout << "\n Iteration: " << i << "\n";
    //     cout << "\n Upper: " << sep.zeta_u << "\n";
    //     cout << "\n Lower: " << sep.zeta_l << "\n";

    //     M2Bendersmodel->optimize();

    //     ++i;
    // }

    sep.xprime[0] = M2Bendersmodel->get(GRB_DoubleAttr_ObjVal);
    for (int a = 0; a < m; ++a)
    {
        try
        {
            sep.xprime[a + 1] = x[a].get(GRB_DoubleAttr_X);
        }
        catch (GRBException e)
        {
            std::cout << "x_" << a << "\n";
            std::cout << "Gurobi error number [M2Benders - retrieve values for xprime]: " << e.getErrorCode() << "\n";
            std::cout << e.getMessage() << "\n";
        }
        catch (...)
        {
            std::cout << "Non-gurobi error during optimization [M2Benders - retrieve values for xprime]"
                      << "\n";
        }
        // sep.xprime[a] = 0;
    }

    // for submodel testing and shit
    std::vector<int> test_xhat = {1, 1, 0};

    // sep.subproblem.Submodel->write("spmodel.lp");
    // sep.yhat = sep.subproblem.solve(0);
    // sep.subproblem.update(test_xhat);
    // sep.subproblem.Submodel->write("spmodelupdated.lp");

    delete sep.subproblem.obj_constr;
    delete sep.subproblem.Subenv;
    delete sep.subproblem.Submodel;

    return sep.xprime;
}

// float BendersSPSub::solve()
// {
// }
