#include "../inc/rspi.h"

rspiModel::rspiModel()
{
    /* 
        Default constructor
    */
    n, m, k, t = 0;
}

rspiModel::rspiModel(const std::string &filename)
{
    /* 
        Read the graph from a file - construct RSPI class
        Construct base gurobi model
    */
}

void rspiModel::update()
{
    /*
        Add main constraints to the gurobi model (depend on subset of cost vectors)
    */
}

void rspiModel::solve(std::vector<int> &data)
{
    /*
        Solve the robust model for each subset of costs in the k-partition
        - If the minimum out of the k solutions is better than our current v, update v, and update x_bar
        - If not, do nothing
        - In the bigger scheme of this package, this function runs at the bottom of the partitioning recursion
    */
    a = data.data();
    next = 0;
    z = -1;

    for (int j = 0; j < k; j++)
    // there are k sets of scenarios in this partition
    {
        // update current_scenarios to reflect new set, based on data, which defines the partition split
        current_scenarios.clear();
        current_scenarios.push_back(next);
        next++;
        while (*a == 0)
        {
            current_scenarios.push_back(next);
            ++a;
            next++;
        }

        update(); // we need update the gurobi model with this subset of scenarios
        model->optimize();

        if (z < 0)
        {
            z = model->get(GRB_DoubleAttr_ObjVal);
        }
        else if (model->get(GRB_DoubleAttr_ObjVal) < z)
        {
            z = model->get(GRB_DoubleAttr_ObjVal);
        }
    }

    if (z > v)
    {
        v = z;
    }
}

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
        p = M2Instance->p;
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

        // // for (int q = 0; q < p; q++)
        // // {
        // //     vector<int> new_vector = {};
        // //     arc_costs.push_back(new_vector);
        // //     for (int a = 0; a < m; a++)
        // //     {
        // //         arc_costs[q].push_back(distr(gen)); // assign arc cost between min and max
        // //     }
        // // }

        // // hardcoded example "simplegraph.txt"
        // vector<int> costs1 = {3, 4, 3, 4, 3, 4, 3};
        // vector<int> costs2 = {3, 1, 3, 1, 3, 1, 10};
        // vector<int> costs3 = {3, 4, 3, 4, 3, 4, 10};
        // arc_costs.push_back(costs1);
        // arc_costs.push_back(costs2);
        // arc_costs.push_back(costs3);

        // ------ Initialize model and environment ------
        M2env = new GRBEnv();
        M2model = new GRBModel(*M2env);

        // ------ Decision variables ------
        string varname;

        // interdiction policy on arcs 'x'
        for (int a = 0; a < m; a++)
        {
            varname = "x_" + to_string(a);
            x.push_back(M2model->addVar(0, 1, 0, GRB_BINARY, varname));
        }

        // convex combination lambda on scenarios
        for (int q = 0; q < p; q++)
        {
            varname = "lambda_" + to_string(q);
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
        for (int q = 0; q < p; q++)
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
            for (int q = 0; q < p; q++)
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
