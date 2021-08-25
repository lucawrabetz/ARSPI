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

LayerGraph::LayerGraph(const string &filename, int the_n)
{
    string line;
    ifstream myfile(filename);

    if (myfile.is_open())
    {

        n = the_n;
        const char *cline;
        int i;
        int j;
        vector<int> new_vector;

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
    cout << "n: " << n << ", m: " << m << "\n";
    for (int a = 0; a < m; a++)
    {
        cout << "(" << arcs[a].i << "," << arcs[a].j << ")\n";
    }
}

M2ProblemInstance::M2ProblemInstance()
{
    // Default constructor
    r_0 = 0;
}

M2ProblemInstance::M2ProblemInstance(const LayerGraph &the_G, int min, int max, int the_p, int the_r0)
{
    // ------ Assign graph and random costs ------
    // ------ Variables and int parameters ------
    G = the_G;
    n = G.n;
    m = G.m;
    p = the_p;
    r_0 = the_r0;

    for (int a = 0; a < m; a++)
    {
        // ASSUMING FINITE INTERDICTION COST, A REASONABLE NUMBER IS MAX - MIN
        // interdiction_costs.push_back((max - min));
        // ASSUMING INFINITE INTERDICTION COST (REMOVING ARC) ADJUST SO LARGE ENOUGH
        interdiction_costs.push_back(100000);
    }

    std::random_device rd;                           // obtain a random number from hardware
    std::mt19937 gen(rd());                          // seed the generator
    std::uniform_int_distribution<> distr(min, max); // define the range

    for (int q = 0; q < p; q++)
    {
        vector<int> new_vector = {};
        arc_costs.push_back(new_vector);
        for (int a = 0; a < m; a++)
        {
            arc_costs[q].push_back(distr(gen)); // assign arc cost between min and max
        }
    }

    for (int q = 0; q < p; q++)
    {
        // cout << "q: " << q << endl;
        for (int a = 0; a < m; a++)
        {
            // cout << "a: " << a << endl;
            // cout << arc_costs[q][a] << "\n";
        }
    }

    // hardcoded example "simplegraph.txt"
    // vector<int> costs1 = {9, 12, 3};
    // vector<int> costs2 = {9, 1, 10};
    // vector<int> costs3 = {9, 12, 10};
    // arc_costs.push_back(costs1);
    // arc_costs.push_back(costs2);
    // arc_costs.push_back(costs3);
}

// ------ MIP Formulations for M2 ------

M2ModelLinear::M2ModelLinear(){
    int n = 0;
    int m = 0;
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
        p = M2Instance->p;
        r_0 = M2Instance->r_0;

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

        // objective func dummy 'z'
        varname = "z";
        z = M2model->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS, varname);

        vector<GRBVar> new_vector;
        // post interdiction flow 'pi'
        for (int q = 0; q < p; q++)
        {
            new_vector = {};
            pi.push_back(new_vector);
            for (int i = 0; i < n; i++)
            {
                varname = "pi_" + to_string(q) + "_" + to_string(i);
                pi[q].push_back(M2model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
            }
        }

        // arc variable 'lambda'
        for (int q = 0; q < p; q++)
        {
            new_vector = {};
            lambda.push_back(new_vector);
            for (int a = 0; a < m; a++)
            {
                varname = "lambda_" + to_string(q) + "_" + to_string(a);
                lambda[q].push_back(M2model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
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
        for (int q = 0; q < p; q++)
        {
            linexpr = 0;
            linexpr += (pi[q][n - 1] - pi[q][s]); // b^\top pi (our b is simply a source-sink single unit of flow)
            for (int a = 0; a < m; ++a)
            {
                linexpr += -lambda[q][a]; // u^\top \cdot lambda (our u is 1)
            }
            M2model->addConstr(z <= linexpr);
        }

        // main constraint for each arc
        int i;
        int j;
        for (int q = 0; q < p; ++q)
        {
            linexpr = 0;
            for (int a = 0; a < m; ++a)
            {
                i = M2Instance->G.arcs[a].i;
                j = M2Instance->G.arcs[a].j;
                M2model->addConstr((pi[q][j] - pi[q][i] - lambda[q][a]) <= M2Instance->arc_costs[q][a] + (M2Instance->interdiction_costs[a] * x[a]));
            }
        }

        // pi[0] = 0
        for (int q = 0; q < p; q++)
        {
            M2model->addConstr(pi[q][0] == 0);
        }
        M2model->update();
    }
    catch (GRBException e)
    {
        cout << "Gurobi error number [M2Model, constructor]: " << e.getErrorCode() << "\n";
        cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        cout << "Non-gurobi error during optimization [M2Model]"
                  << "\n";
    }
}

vector<float> M2ModelLinear::solve()
{
    /*
     * Returns vector of objective value [0] and interdiction policy [1-m+1]
     * Note: the return values are stored in class variables so unassigned in some compexp runs 
     */

    try
    {
        clock_t model_begin = clock();
        M2model->optimize();
        running_time = float(clock() - model_begin) / CLOCKS_PER_SEC;
        optimality_gap = M2model->get(GRB_DoubleAttr_MIPGap);

        cout << "In MIP.solve() method" << endl;
        cout << "Objective: " << M2model->get(GRB_DoubleAttr_ObjVal) << "\n";
        x_prime.push_back(M2model->get(GRB_DoubleAttr_ObjVal));
        cout << "Running time: " << running_time << "\n";

        for (int a = 0; a < m; a++)
        {
            cout << "x_" << a << "(" << M2Instance->G.arcs[a].i << "," << M2Instance->G.arcs[a].j << ")"
                      << ": " << x[a].get(GRB_DoubleAttr_X) << "\n";
            x_prime.push_back(x[a].get(GRB_DoubleAttr_X));
        }

        return x_prime;

    }
    catch (GRBException e)
    {
        cout << "Gurobi error number [M2Model, solveMIP]: " << e.getErrorCode() << "\n";
        cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        cout << "Non-gurobi error during optimization [M2Model]"
                  << "\n";
    }
}

// ------ Bender's Schemes for M2 ------
BendersSub::BendersSub()
{
    n = 0;
    m = 0;
    p = 0;
}

BendersSub::BendersSub(M2ProblemInstance *the_M2Instance)
{
    // ------ Initialize Basic Parameters ------
    n = the_M2Instance->n;
    m = the_M2Instance->m;
    p = the_M2Instance->p;

    // ------ Initialize d and c costs ------
    for (int a = 0; a < m; a++)
    {
        d.push_back(the_M2Instance->interdiction_costs[a]);
    }

    for (int q = 0; q < p; q++)
    {
        c.push_back(the_M2Instance->arc_costs[q]);
        c_bar.push_back(the_M2Instance->arc_costs[q]);
    }

    // ------ Initialize Environment and Model ------
    for (int q = 0; q < p; ++q)
    {
        Subenvs.push_back(new GRBEnv());
        Submodels.push_back(new GRBModel(*Subenvs[q]));
    }

    // ------ Decision Variables ------
    for (int q = 0; q < p; ++q)
    {
        varname = "zeta_sub_" + to_string(q);
        zeta_subs.push_back(Submodels[q]->addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS, varname));
    }
    for (int q = 0; q < p; ++q)
    {
        y_dummy = {};
        y.push_back(y_dummy);
        for (int a = 0; a < m; ++a)
        {
            varname = "y_" + to_string(q) + "_" + to_string(a);
            y[q].push_back(Submodels[q]->addVar(0, 1, 0, GRB_CONTINUOUS, varname));
        }
    }

    // ------ Constraints ------
    // constraints to bound objective value over q
    obj_constr = new GRBConstr[p];
    for (int q = 0; q < p; ++q)
    {
        linexpr = 0;
        for (int a = 0; a < m; ++a)
        {
            linexpr += c_bar[q][a] * y[q][a];
        }
        obj_constr[q] = Submodels[q]->addConstr(zeta_subs[q] >= linexpr);
    }

    // flow constraints
    for (int q = 0; q < p; ++q)
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
            Submodels[q]->addConstr(linexpr == rhs);
        }
        Submodels[q]->update();
    }
}

void BendersSub::update(vector<int> &xhat)
{
    cout << "\nsubmodel: updating based on xbar, new interdiction policy: \n";
    for (int a = 0; a < m; ++a)
    {
        cout << "xbar_" << a << ": " << xhat[a] << "\n";
    }

    // update array of constraints instead of only the parameter vector
    for (int q = 0; q < p; ++q)
    {
        for (int a = 0; a < m; ++a)
        {
            // model.chgCoeff(obj_constr[q], variable, new cost)
            if (xhat[a] > 0.5)
            {
                Submodels[q]->chgCoeff(obj_constr[q], y[q][a], -(c[q][a] + d[a]));
            }
            else
            {
                Submodels[q]->chgCoeff(obj_constr[q], y[q][a], -(c[q][a]));
            }
        }
        Submodels[q]->update();
    }
}

vector<vector<float>> BendersSub::solve(int counter)
{
    try
    {
        // yhat has q elements - each one has m+1 elements (first is the objective, rest is the flow)
        vector<vector<float>> yhat;

        for (int q = 0; q < p; ++q)
        {
            y_dummy2 = {};
            yhat.push_back(y_dummy2);
            string modelname = "submodel_" + to_string(counter) + "q=" + to_string(q) + ".lp";
            Submodels[q]->write(modelname);
            Submodels[q]->optimize();

            yhat[q].push_back(Submodels[q]->get(GRB_DoubleAttr_ObjVal));
            cout << "\nq = " + to_string(q);
            cout << "\nsubmodel obj: " << yhat[q][0];
            cout << "\narc values: \n";
            for (int a = 0; a < m; ++a)
            {
                yhat[q].push_back(y[q][a].get(GRB_DoubleAttr_X));
                cout << "y_" << q << "_" << a << ": " << y[q][a].get(GRB_DoubleAttr_X) << "\n";
            }
        }
        return yhat;
    }
    catch (GRBException e)
    {
        cout << "Gurobi error number [BendersSub::solve]: " << e.getErrorCode() << "\n";
        cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        cout << "Non-gurobi error during optimization [BendersSub::Solve]"
                  << "\n";
    }
}

BendersSeparation::BendersSeparation()
{
    n = 0;
    m = 0;
    p = 0;
}

BendersSeparation::BendersSeparation(GRBVar &the_zetabar, vector<GRBVar> &the_xbar, M2ProblemInstance *the_M2Instance)
{
    try
    {
        // ------ Initialize Basic Parameters ------
        n = the_M2Instance->n;
        m = the_M2Instance->m;
        p = the_M2Instance->p;

        // ------ Initialize Submodel ------
        subproblem = BendersSub(the_M2Instance);

        // ------ Initialize d and c costs ------
        for (int a = 0; a < m; a++)
        {
            d.push_back(the_M2Instance->interdiction_costs[a]);
        }

        for (int q = 0; q < p; q++)
        {
            c.push_back(the_M2Instance->arc_costs[q]);
        }

        // ------ Initialize Variable containers ------
        xbar = the_xbar;
        zetabar = the_zetabar;
        xprime.push_back(0);

        for (int q = 0; q < p; ++q)
        {
            if (q == 0)
            {
                for (int a = 0; a < m; ++a)
                {
                    xhat.push_back(0);
                    xprime.push_back(0);
                }
            }
        }
    }
    catch (GRBException e)
    {
        cout << "Gurobi error number [BendersSeparation, constructor]: " << e.getErrorCode() << "\n";
        cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        cout << "Non-gurobi error during optimization [BendersSeparation]"
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
            zeta_temp = GRB_INFINITY;
            for (int q = 0; q < p; ++q)
            {
                if (zeta_temp > yhat[q][0])
                {
                    zeta_temp = yhat[q][0]; // first element of yhat[q] is the objective
                }
            }
            // use yhat[(q-l)+1][1-m] to create new cut from LinExpr
            for (int q = 0; q < p; ++q)
            {
                new_cut = 0;
                for (int a = 0; a < m; ++a)
                {
                    new_cut += (c[q][a] + d[a] * xbar[a]) * yhat[q][a + 1];
                }
                // add lazy cut to main model
                try
                {
                    addLazy(zetabar <= new_cut);
                    cout << "\nadded cut: "
                         << "zetabar"
                         << "<=" << new_cut << "\n";
                    ++cut_count;
                }
                catch (GRBException e)
                {
                    cout << "Gurobi error number [BendersSeparation, addLazy]: " << e.getErrorCode() << "\n";
                    cout << e.getMessage() << "\n";
                }
            }

            if (zeta_l < zeta_temp)
            {
                zeta_l = zeta_temp;
            }
        }
    }
}

M2Benders::M2Benders()
{
    int n = 0;
    int m = 0;
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
        p = M2Instance->p;
        r_0 = M2Instance->r_0;

        // ------ Decision variables ------
        string varname;

        // interdiction variable 'x'
        for (int a = 0; a < m; a++)
        {
            varname = "x_" + to_string(a);
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
        for (int q = 0; q < p; ++q)
        {
            linexpr = 0;
            for (int a = 0; a < m; ++a)
            {
                linexpr += (sep.c[q][a] + (sep.d[a] * x[a])) * sep.yhat[q][a + 1];
            }
            M2Bendersmodel->addConstr(zeta <= linexpr);
        }
    }
    catch (GRBException e)
    {
        cout << "Gurobi error number [M2Benders, constructor]: " << e.getErrorCode() << "\n";
        cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        cout << "Non-gurobi error during optimization [BendersSPSub]"
                  << "\n";
    }
}

vector<float> M2Benders::solve()
{
    /*
     * Currently the return vector returns the objective value [0] and interdiction policy [1-m+1]
     * Again for computational experiments this can be ignored and not assigned as run stats are class vars
     */ 

    // ------ Set Callback on Master Model
    M2Bendersmodel->setCallback(&sep);

    // ------ Optimize Inside Benders Scheme -------
    M2Bendersmodel->write("simplegraph1_benders_before.lp");

    try
    {
        clock_t model_begin = clock();
        M2Bendersmodel->optimize();
        running_time = float(clock() - model_begin) / CLOCKS_PER_SEC;
        optimality_gap = M2Bendersmodel->get(GRB_DoubleAttr_MIPGap);
        cut_count = sep.cut_count;
    }
    catch (GRBException e)
    {
        cout << "Gurobi error number [M2Benders.optimize()]: " << e.getErrorCode() << "\n";
        cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        cout << "Non-gurobi error during optimization [M2Benders.optimize()]"
                  << "\n";
    }
    M2Bendersmodel->write("simplegraph1_benders_after.lp");

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

    // for submodel testing and stuff
    // vector<int> test_xhat = {1, 1, 0};

    // sep.subproblem.Submodel->write("spmodel.lp");
    // sep.yhat = sep.subproblem.solve(0);
    // sep.subproblem.update(test_xhat);
    // sep.subproblem.Submodel->write("spmodelupdated.lp");

    delete sep.subproblem.obj_constr;
    for (int q = 0; q < p; ++q)
    {
        delete sep.subproblem.Subenvs[q];
        delete sep.subproblem.Submodels[q];
    }
    return sep.xprime;
}

