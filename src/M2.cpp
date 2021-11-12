#include "../inc/M2.h"

Arc::Arc(){i=0; j=0;}

Arc::Arc(int the_i, int the_j){i=the_i; j=the_j;}

LayerGraph::LayerGraph(){n=0; m=0;}

LayerGraph::LayerGraph(const string &filename, int the_n)
{
    // LayerGraph constructor from file
    string line;
    ifstream myfile(filename);

    if (myfile.is_open())
    {
        n = the_n;
        int counter = 0;
        int i;
        int j;
        vector<int> new_vector;
        vector<int> zeros_vector;
        const char *cline;

        for (int i = 0; i < n; i++) {zeros_vector.push_back(0);}

        for (int i = 0; i < n; i++)
        {
            new_vector = {};
            adjacency_list.push_back(new_vector);
            arc_index_hash.push_back(new_vector);
            n_n_adjacency_list.push_back(zeros_vector);
        }

        while (getline(myfile, line))
        {
            cline = line.c_str();
            sscanf(cline, "%d %d", &i, &j);

            arcs.push_back(Arc(i, j));
            adjacency_list[i].push_back(j);
            arc_index_hash[i].push_back(counter);
            n_n_adjacency_list[i][j] = 1;
            ++counter; // to track the index in the 0-(m-1) vectors
        }

        m = arcs.size();
    }
}

void LayerGraph::printGraph(vector<vector<int> > costs, vector<int> interdiction_costs, bool is_costs) const
{
    // Print arc summary of a graph, with costs if called from M2 Instance (is_costs)
    cout << "n: " << n << ", m: " << m << endl;
    int p = costs.size();

    if (is_costs){
        for (int a = 0; a < m; a++)
        {
            cout << "(" << arcs[a].i << "," << arcs[a].j << ")" << endl;
            cout << "   interdiction_cost: " << interdiction_costs[a] << endl;

            for (int q = 0; q<p; ++q){
                cout << "   arc_cost_" << q << ": " << costs[q][a] << endl;
            }
        }
    }

    else {
        for (int a = 0; a < m; a++)
        {
            cout << "(" << arcs[a].i << "," << arcs[a].j << ")\n";
        }
    }
}

M2ProblemInstance::M2ProblemInstance(){r_0=0;}

M2ProblemInstance::M2ProblemInstance(const LayerGraph &the_G, int min, int max, int the_p, int the_k, int the_r0, string& the_instance_name, string& the_setname)
{
    // ------ Assign graph and random costs ------
    // ------ Variables and int parameters ------
    G = the_G;
    n = G.n;
    m = G.m;
    p = the_p;
    k = the_k;
    r_0 = the_r0;
    instance_name = the_instance_name;
    setname = the_setname;

    for (int a = 0; a < m; a++)
    {
        // ASSUMING FINITE INTERDICTION COST, A REASONABLE NUMBER IS MAX - MIN
        // interdiction_costs.push_back((max - min));
        // ASSUMING INFINITE INTERDICTION COST (REMOVING ARC) ADJUST SO LARGE ENOUGH
        interdiction_costs.push_back(10);
    }

    std::random_device rd;                           // obtain a random number from hardware
    std::mt19937 gen(rd());                          // seed the generator
    std::uniform_int_distribution<> distr(min, max); // define the range

    for (int q = 0; q < p; q++)
    {
        instance_name = the_instance_name;
        vector<int> new_vector = {};
        arc_costs.push_back(new_vector);
        for (int a = 0; a < m; a++)
        {
            arc_costs[q].push_back(distr(gen)); // assign arc cost between min and max
        }
    }

    // for (int q = 0; q < p; q++)
    // {
    //     // cout << "q: " << q << endl;
    //     for (int a = 0; a < m; a++)
    //     {
    //         // cout << "a: " << a << endl;
    //         // cout << arc_costs[q][a] << "\n";
    //     }
    // }

    // hardcoded example "simplegraph.txt"
    // vector<int> costs1 = {6, 4, 2, 2, 1, 2, 7, 1, 3};
    // arc_costs.push_back(costs1);
}

void M2ProblemInstance::printInstance() const {
    // Print Summary of Problem Instance
    cout << "k: " << k << endl;
    cout << "p: " << p << endl;
    G.printGraph(arc_costs, interdiction_costs, true);
}

vector<int> M2ProblemInstance::Dijkstra(int q)
{
    vector<int> bar_S;
    vector<int> S;
    vector<int> dist;
    vector<int> pred;
    vector<int> result;
    
    int min;
    int node;
    int j_node;
    int arc;
    int final_cost;

    // d(s) = 0
    // pred(s) = 0
    bar_S.push_back(0);
    dist.push_back(0);
    pred.push_back(0);

    for (int a=0; a<m+1; ++a){
        // result[0] = path objective
        // result[a+1] = 1 if arc in path
        result.push_back(0);
    }

    for (int i=1; i<n; ++i) {
        // d(i) = inf
        // pred(s) = -1
        bar_S.push_back(i);
        dist.push_back(interdiction_costs[0]*10);
        pred.push_back(-1);
    }

    while (S.size()<n) {
        min = interdiction_costs[0]*10 + 1;
        int erase_index;

        for (int i=0; i<bar_S.size(); ++i) {
            if (dist[bar_S[i]] < min) {
                erase_index = i;
                node=bar_S[i];
                min=dist[bar_S[i]];
            }
        }

        // remove node from bar_S
        bar_S.erase(bar_S.begin()+erase_index); 
        
        // add node to S
        S.push_back(node);
        // cout << "i: " << node << endl;

        for (int i=0; i<G.arc_index_hash[node].size(); ++i){
            arc=G.arc_index_hash[node][i];
            // cout << "arc: " << arc << " ; " << G.arcs[arc].i << " " << G.arcs[arc].j << endl; 
            j_node=G.adjacency_list[node][i];
            // cout << "j: " << j_node << endl;

            if (dist[j_node] > (dist[node]+arc_costs[q][arc])){
                dist[j_node] = dist[node]+arc_costs[q][arc];
                // cout << "d[j]: " << dist[j_node] << endl;
                pred[j_node] = node;
                // cout << "pred[j]: " << pred[j_node] << endl;
            }
        }
        // cout << "\n\n\n";
    }

    final_cost=dist[n-1];
    // cout << "final " << final_cost << endl;
    result[0]=final_cost;
    j_node=n-1;

    while(j_node!=0){
        // find the arc index for the arc i, j (looping backwards through sp tree using pred to determine i)
        node=pred[j_node];

        for (int i=0; i<G.adjacency_list[node].size(); ++i){
            if (G.adjacency_list[node][i]==j_node){
                arc=G.arc_index_hash[node][i];
                break;
            }
        }

        result[arc+1]=1;
        j_node=node;
    }

    return result;
}

void M2ProblemInstance::updateCosts(vector<float>& x_bar, bool rev){
    // Receive interdiction policy and update costs for the M2 instance
    // If rev, we are "removing" the interdiction policy and returning the instance to its original state
    for (int a=0; a<m; ++a){
        if (x_bar[a]==1) {
            for (int q=0; q<p; ++q){
                if (rev){
                    arc_costs[q][a]-=interdiction_costs[a];
                }
                else {
                    arc_costs[q][a]+=interdiction_costs[a];
                }
            }
        }
    }
}

float M2ProblemInstance::validatePolicy(vector<float>& x_bar)
{
    float objective=100000000;
    vector<int> sp_result;

    // Update M2 based on x_bar
    updateCosts(x_bar);

    // run dijstra on graph to get objective 
    for (int q=0; q<p; ++q){
        sp_result = Dijkstra(q);
        
        if (sp_result[0]<objective){
            objective=sp_result[0];
        }
    }
    
    // Revert M2 
    updateCosts(x_bar, true);

    return objective;
}

// ------ MIP Formulations for M2 ------
M2ModelLinear::M2ModelLinear(){n=0; m=0;}

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
        k = M2Instance->k;
        r_0 = M2Instance->r_0;
        instance_name = M2Instance->instance_name;
        setname = M2Instance->setname;

        // ------ Initialize model and environment ------
        M2env = new GRBEnv();
        M2model = new GRBModel(*M2env);

        // M2model->set(GRB_IntParam_OutputFlag, 0);
        M2model->set(GRB_DoubleParam_TimeLimit, 3600);

        // ------ Decision variables ------
        string varname;

        vector<GRBVar> new_vector;

        // set partitioning variables
        for (int w = 0; w < k; ++w){
            cout << "in model constructor" << endl;
            H.push_back(new_vector);
            for (int q =0; q<p; ++q){
                varname = "H_" + to_string(w) + "_" + to_string(q);
                H[w].push_back(M2model->addVar(0, 1, 0, GRB_BINARY, varname));
            }
        }

        // interdiction policy on arcs 'x'
        x = vector<vector<GRBVar> >(k, vector<GRBVar>(m, M2model->addVar(0, 1, 0, GRB_BINARY)));
//        for (int w = 0; w < k; ++w) {
         //   new_vector = {};
         //   x.push_back(new_vector);

         //   for (int a = 0; a < m; a++)
         //   {
         //       varname = "x_" + to_string(w) + "_" + to_string(a);
         //       x[w].push_back());
         //   }
         //   }
        //

        // objective func dummy 'z'
        z = M2model->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS);

        vector<vector<GRBVar> > newnew_vector;
        
        // post interdiction flow 'pi'
        for (int w = 0; w<k; ++w){
            newnew_vector = {};
            pi.push_back(newnew_vector);

            for (int q = 0; q < p; q++)
            {
                new_vector = {};
                pi[w].push_back(new_vector);

                for (int i = 0; i < n; i++)
                {
                    varname = "pi_" + to_string(w) + "_" + to_string(q) + "_" + to_string(i);
                    pi[w][q].push_back(M2model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
                }
            }
        }
        // arc variable 'lambda'
        for (int w = 0; w<k; ++w){
            newnew_vector = {};
            lambda.push_back(newnew_vector);

            for (int q = 0; q < p; q++)
            {
                new_vector = {};
                lambda[w].push_back(new_vector);
                
                for (int a = 0; a < m; a++)
                {
                    varname = "lambda_" + to_string(w) + "_"  + to_string(q) + "_" + to_string(a);
                    lambda[w][q].push_back(M2model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
                }
            }
        }

        // ------ Constraints ------
        // budget constraint
        for (int w = 0; w<k; ++w){
            linexpr = 0;
            for (int a = 0; a < m; a++)
            {
                linexpr += x[w][a];
            }
            M2model->addConstr(linexpr <= r_0);
        }

        // z constraints
        for (int w = 0; w<k; ++w){
            for (int q = 0; q < p; q++)
            {
                linexpr = 0;
                linexpr += (pi[w][q][n - 1] - pi[w][q][s]); // b^\top pi (our b is simply a source-sink single unit of flow)
                for (int a = 0; a < m; ++a)
                {
                    linexpr += -lambda[w][q][a]; // u^\top \cdot lambda (our u is 1)
                }

                linexpr += M * (1 - H[w][q]);
                M2model->addConstr(z <= linexpr);
            }
        }

        // main constraint for each arc
        int i;
        int j; // looping index 
        int jn; // node j
        int a;

        for (int w = 0; w<k; ++w){
            for (int q=0; q<p; ++q){
                for (i=0; i<n; ++i){
                    for (j=0; j<M2Instance->G.adjacency_list[i].size(); ++j){
                        jn=M2Instance->G.adjacency_list[i][j];
                        a=M2Instance->G.arc_index_hash[i][j];

                        // add constraint
                        M2model->addConstr((pi[w][q][jn] - pi[w][q][i] - lambda[w][q][a]) <= M2Instance->arc_costs[q][a] + (M2Instance->interdiction_costs[a] * x[w][a]));
                    }
                }
            }
        }

        // for (int q = 0; q < p; ++q)
        // {
        //     linexpr = 0;
        //     for (int a = 0; a < m; ++a)
        //     {
        //         i = M2Instance->G.arcs[a].i;
        //         j = M2Instance->G.arcs[a].j;
        //         M2model->addConstr((pi[q][j] - pi[q][i] - lambda[q][a]) <= M2Instance->arc_costs[q][a] + (M2Instance->interdiction_costs[a] * x[a]));
        //     }
        // }

        // set-partitioning constraint
        for (int q=0; q<p; ++q){
            linexpr=0;

            for (int w=0; w<k; ++w){
                linexpr+=H[w][q];
            }

            M2model->addConstr(linexpr == 1);
        }

        // pi[0] = 0
        for (int w=0; w<k; ++w) {
            for (int q = 0; q < p; q++)
            {
                M2model->addConstr(pi[w][q][0] == 0);
            }
        }
        M2model->update();

        modelname = "modelfiles/" + setname + "/" + instance_name + "-M2Model.lp";
        M2model->write(modelname);
        // cout << modelname << endl;

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

vector<vector<float> > M2ModelLinear::solve()
{
    /*
     * Returns vector of interdiction policy [0-m] (for every (w) for M3, plus an extra singleton as first element for the objective value)
     * Note: the return values are stored in class variables so unassigned in some compexp runs 
     */

    try
    {
        // cout << "In MIP.solve() method" << endl;

        cout << "in solve" << endl;

        clock_t model_begin = clock();
        M2model->optimize();
        running_time = float(clock() - model_begin) / CLOCKS_PER_SEC;

        cout << "just optimized" << endl;

        try {

            vector<float> x_vector; 
            optimality_gap = M2model->get(GRB_DoubleAttr_MIPGap);
            // cout << "Objective: " << M2model->get(GRB_DoubleAttr_ObjVal) << "\n";
            
            vector<float> objective_vec;
            objective_vec.push_back(M2model->get(GRB_DoubleAttr_ObjVal));
            x_prime.push_back(objective_vec);

            for (int w=0; w<k; ++w){
                x_vector = {};

                for (int a = 0; a < m; a++)
                {
                    // cout << "x_" << a << "(" << M2Instance->G.arcs[a].i << "," << M2Instance->G.arcs[a].j << ")"
                              // << ": " << x[a].get(GRB_DoubleAttr_X) << "\n";
                    x_vector.push_back(x[w][a].get(GRB_DoubleAttr_X));
                }
                x_prime.push_back(x_vector);
            }
        }

        catch (GRBException e) {

            if (e.getErrorCode() == 10005 || e.getErrorCode() == 10013) {
                // in this case the error is because the model was unbounded
                // set the optimality gap to -1 and we'll list as unbounded 
                // put -infinity as objective value (the objectives are negative, min -z)
                // set arc interdiction values as -1 - gurobi won't give us a solution one unboundedness is proven
                optimality_gap = -1;
                // cout << "Objective: unbounded" << "\n";
                vector<float> vec; 
                vec.push_back(-GRB_INFINITY);
                x_prime.push_back(vec);
            }

            else {
                cout << "Gurobi error number [M2Linear.optimize()]: " << e.getErrorCode() << "\n";
                cout << e.getMessage() << "\n";
            }
        }

        // cout << "Running time: " << running_time << "\n";


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
BendersSub::BendersSub(){n=0; m=0; p=0;}

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
        Submodels[q]->set(GRB_IntParam_OutputFlag, 0);
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
    // cout << "\nsubmodel: updating based on xbar, new interdiction policy: \n";
    // for (int a = 0; a < m; ++a)
    // {
    //     cout << "xbar_" << a << ": " << xhat[a] << "\n";
    // }

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

vector<vector<float> > BendersSub::solve(int counter)
{
    try
    {
        // yhat has q elements - each one has m+1 elements (first is the objective, rest is the flow)
        vector<vector<float> > yhat;

        for (int q = 0; q < p; ++q)
        {
            y_dummy2 = {};
            yhat.push_back(y_dummy2);
            // string modelname = "submodel_" + to_string(counter) + "q=" + to_string(q) + ".lp";
            // Submodels[q]->write(modelname);
            Submodels[q]->optimize();

            yhat[q].push_back(Submodels[q]->get(GRB_DoubleAttr_ObjVal));
            // cout << "\nq = " + to_string(q);
            // cout << "\nsubmodel obj: " << yhat[q][0];
            // cout << "\narc values: \n";
            for (int a = 0; a < m; ++a)
            {
                yhat[q].push_back(y[q][a].get(GRB_DoubleAttr_X));
                // cout << "y_" << q << "_" << a << ": " << y[q][a].get(GRB_DoubleAttr_X) << "\n";
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
        // zeta_u = -getDoubleInfo(GRB_CB_MIPSOL_OBJBST); // ??????? 

        if (zeta_u - zeta_l >= epsilon)
        {
            cout << "in callback" << endl;
            counter++;
            // xhat = current solution from master problem, then update subproblem
            for (int a = 0; a < m; ++a)
            {
                xhat[a] = getSolution(xbar[a]);
            }
            subproblem.update(xhat);
            // zeta_u = -getDoubleInfo(GRB_CB_MIPSOL_OBJ); // ??????? 

            // cout << "\n\n\n\nsolving sub from callback: \n\n\n\n";
            yhat = subproblem.solve(counter);
            zeta_temp = GRB_INFINITY;

            for (int q = 0; q < p; ++q)
            {
                cout << "sub q=" << q << ": " << yhat[q][0] << endl;
                if (zeta_temp > yhat[q][0])
                {
                    zeta_temp = yhat[q][0]; // first element of yhat[q] is the objective
                }
            }
            cout << "zeta_temp: " << zeta_temp << endl;

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
                    // cout << "\nadded cut: "
                         // << "zetabar"
                         // << "<=" << new_cut << "\n";
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
            cout << "zeta_l: " << zeta_l << endl;
            zeta_u = -getDoubleInfo(GRB_CB_MIPSOL_OBJBND); // ??????? 
            cout << "zeta_u: " << zeta_u << endl;
        }
    }
}

M2Benders::M2Benders(){int n, m=0;}

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
        M2Bendersmodel->set(GRB_IntParam_OutputFlag, 0);
        M2Bendersmodel->set(GRB_DoubleParam_TimeLimit, 3600);

        // ------ Variables and int parameters ------
        n = M2Instance->n;
        m = M2Instance->m;
        p = M2Instance->p;
        r_0 = M2Instance->r_0;
        instance_name = M2Instance->instance_name;
        setname = M2Instance->setname;

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
        // cout << "\n\n\n\nsolving sub from constructor: \n\n\n\n";
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

        modelname = "modelfiles/" + setname + "/" + instance_name + "-M2BendersModel.lp";
        M2Bendersmodel->write(modelname);

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
     * Return vector returns the objective value [0] and interdiction policy [1-m+1]
     * Again for computational experiments this can be ignored and not assigned as run stats are class vars
     */ 

    // ------ Set Callback on Master Model
    M2Bendersmodel->setCallback(&sep);


    // ------ Optimize Inside Benders Scheme -------

    try
    {
        clock_t model_begin = clock();
        M2Bendersmodel->optimize();
        running_time = float(clock() - model_begin) / CLOCKS_PER_SEC;

        try {
            optimality_gap = M2Bendersmodel->get(GRB_DoubleAttr_MIPGap);
            sep.xprime[0] = M2Bendersmodel->get(GRB_DoubleAttr_ObjVal);

            for (int a = 0; a < m; ++a)
            {
                sep.xprime[a + 1] = x[a].get(GRB_DoubleAttr_X);
            }

        }

        catch (GRBException e) {

            if (e.getErrorCode() == 10005 || e.getErrorCode() == 10013) {
                // in this case the error is because the model was unbounded
                // set the optimality gap to -1 and we'll list as unbounded 
                // put objective value as -infinity (objectives are negated min -z)
                // set arc interdiction values as -1 - gurobi won't give us a solution one unboundedness is proven
                optimality_gap = -1;
                sep.xprime[0] = -GRB_INFINITY;

                for (int a = 0; a < m; ++a)
                {
                    sep.xprime[a + 1] = -1;
                }
            }

            else {
                cout << "Gurobi error number [M2Benders.optimize()]: " << e.getErrorCode() << "\n";
                cout << e.getMessage() << "\n";
            }

        }
        
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

vector<int> initKappa(int p, int k) {
    // initialize partition vector based on total number in set (p) and exact number of partitions required
    vector<int> kappa(p, 0);

    for (int q=(p-k+1); q<p; ++q){
        kappa[q]=(q-(p-k));
    }

    return kappa;
}

RobustAlgoModel::RobustAlgoModel(){int n, m=0;}

RobustAlgoModel::RobustAlgoModel(M2ProblemInstance& M2) {
    // Construct baseline model with no constraints
    n = M2.n;
    m = M2.m;
    p = M2.p;
    r_0 = M2.r_0;
    AlgoEnv = new GRBEnv();
    AlgoModel = new GRBModel(*AlgoEnv);

    // Decision Variables
    z = AlgoModel->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS); // objective function dummy variable
    for (int a=0; a<m; ++a) {x.push_back(AlgoModel->addVar(0, 1, 0, GRB_BINARY));} // interdiction policy on arcs

    vector<GRBVar> tempvector;
    for (int q=0; q<p; ++q) {
        pi.push_back(tempvector); // post interdiction s-i best path (for every q)
        for (int i=0; i<n; ++i) {pi[q].push_back(AlgoModel->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS));}
    }

    for (int q=0; q<p; ++q) {
        lambda.push_back(tempvector); // lambda variable on arcs (for every q)
        for (int a=0; a<m; ++a) {lambda[q].push_back(AlgoModel->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS));}
    }
    AlgoModel->update();

    // Populate Global Constraints
    // Arc/Dual Constraints
    for (int q=0; q<p; ++q){
        dual_constraints.push_back(vector<GRBTempConstr>());

        for (int i=0; i<n; ++i) {
            for (int j=0; j<M2.G.adjacency_list[i].size(); ++j){
                int next = M2.G.adjacency_list[i][j];
                int a = M2.G.arc_index_hash[i][j];

                GRBTempConstr constraint = pi[q][next] - pi[q][i] - lambda[q][a] <= M2.arc_costs[q][a] + (M2.interdiction_costs[a]*x[a]);
                dual_constraints[q].push_back(constraint);
            }
        }
    }

    // Objective Constraints
    for (int q=0; q<p; ++q) {
        GRBLinExpr linexpr = 0;
        linexpr += (pi[q][n-1] - pi[q][s]); // b^\top pi

        for (int a=0; a<m; ++a){
            linexpr += -lambda[q][a]; // u^\top \cdot lambda
        }

        z_constraints.push_back(z <= linexpr);
    }
}

void RobustAlgoModel::update(vector<int>& subset) {
    // add constraints to model for subset in partition
    for (int q : subset) {
        string zero_name = "zero_" + to_string(q);
        string z_name = "z_" + to_string(q);
        AlgoModel->addConstr(pi[q][0]==0, zero_name);
        AlgoModel->addConstr(z_constraints[q], z_name);

        int a = 0;
        for (GRBTempConstr constraint : dual_constraints[q]) {
            string dual_name = "dual_" + to_string(q) + "_" + to_string(a);
            AlgoModel->addConstr(constraint, dual_name);
            ++a;
        }
    }
    
    AlgoModel->update();
}

void RobustAlgoModel::reverse_update(vector<int>& subset) {
    // remove all constraints from model - i.e. all constraints associated with this subset
    for (int q : subset) {
        string zero_name = "zero_" + to_string(q);
        string z_name = "z_" + to_string(q);

        GRBConstr zero_constraint = AlgoModel->getConstrByName(zero_name);
        GRBConstr z_constraint = AlgoModel->getConstrByName(z_name);
        AlgoModel->remove(zero_constraint);
        AlgoModel->remove(z_constraint);

        int a = 0;
        for (GRBTempConstr constraint : dual_constraints[q]) {
            string dual_name = "dual_" + to_string(q) + "_" + to_string(a);
            GRBConstr dual_constraint = AlgoModel->getConstrByName(dual_name);
            AlgoModel->remove(dual_constraint);
            ++a;
        }
    }

    AlgoModel->update();
}

vector<double> RobustAlgoModel::solve() {
    // solve static model, return policy and objective value
    AlgoModel->optimize();
    vector<double> sol(m+1, 0);
    
    sol[0] = AlgoModel->get(GRB_DoubleAttr_ObjVal);
    sol[0] = -sol[0];
    for (int a=1; a<m+1; ++a) {
        sol[a] = x[a-1].get(GRB_DoubleAttr_X);
    }

    return sol;
}

int max_int(int a, int b){
    if (a <= b) {return b;}
    else {return a;}
}

bool nextKappa(vector<int>& kappa, vector<int>& max, int k, int p){
    // update kappa and max place to next partition
    // if this is the last one, return false
    
    for (int q=(p-1); q>0; --q){
        if ((kappa[q] < k-1) && (kappa[q] <= max[q-1])){
            
            //cout << "q = " << q << endl; 

            ++kappa[q]; max[q] = max_int(max[q], kappa[q]);

            for (int u=q+1; u<=(p-(k-max[q])); ++u){
                //cout << "first half pass loop" << endl;
                //cout << "u: " << u << endl;
                kappa[u]=0; max[u]=max[q];
            }

            for (int u=(p-(k-max[q]))+1; u<p; ++u){
                // cout << "second half pass loop" << endl;
                // cout << "u: " << u << endl;
                kappa[u] = (k-(p-u));
                max[u] = kappa[u];
            }

            return true;
        }
    }
    return false;
}

vector<vector<int> > kappa_to_partition(vector<int>& kappa, int k, int p){
    // convert a kappa vector to a partition of subsets
    vector<int> subset;
    vector<vector<int> > partition(k, subset);

    for (int q=0; q<p; ++q) {
        int subset_index = kappa[q];
        partition[subset_index].push_back(q);
    }

    return partition;
}

pair<vector<vector<int> >, vector<vector<double> > > enumSolve(M2ProblemInstance& M2){
    // Pass a M2ProblemInstance and solve using enumeration algorithm
    // Initialize and maintain an H matrix representing partitioning 
    // Enumerate all possible partitions by enumerating H matrix
    // Enumeration is done through the 'string' method, and then the H matrix/problem instance are updated (??)
    // Return: pair of nested vectors:
    //      first: ints - optimal partition
    //      second: doubles - k interdiction policies

    // ints and bool
    int k = M2.k;
    int p = M2.p;
    int m = M2.G.m;
    bool next = true;
    
    // initialize partitioning 'string' vector as in Orlov Paper 
    // initialize corresponding max vector 
    vector<int> kappa = initKappa(p, k);
    vector<int> max = kappa;

    // initialize solution vector
    vector<double> arc_vec(m, 0);
    vector<vector<double> > sol(k, arc_vec);

    // initialize static robust model
    try {
        RobustAlgoModel static_robust = RobustAlgoModel(M2);

        // objective value maintained here (and optimal partition)
        double best_worstcase_solution = 0;
        vector<vector<int> > best_worstcase_partition; 

        // enumerate while not 'failing' to get next partition
        while (next) {
            vector<vector<double> > temp_sol(k, arc_vec);
            double temp_worst_sol = DBL_MAX;
            vector<vector<int> > partition = kappa_to_partition(kappa, k, p);

            for (int w=0; w<k; ++w) {
                // for every subset in partition, solve M2 for k=1
                static_robust.update(partition[w]);
                vector<double> temp_single_solution = static_robust.solve();
                static_robust.reverse_update(partition[w]);

                // update temp_worst_sol if this subset is worse
                if (temp_single_solution[0] < temp_worst_sol) {temp_worst_sol = temp_single_solution[0];}
            }

            if (temp_worst_sol > best_worstcase_solution) {
                best_worstcase_solution = temp_worst_sol; 
                sol = temp_sol;
                best_worstcase_partition = partition;
            }

            cout << "next: " << next << endl;
            cout << "[ ";
            for (int q=0; q<p; ++q){cout << kappa[q] << " ";}
            cout << "]" << endl;
            cout << "\nworst case solution (this partition): " << temp_worst_sol << endl;
            cout << "best worst case solution so far: " << best_worstcase_solution << endl;
            next = nextKappa(kappa, max, k, p);
        }


        auto final_solution = make_pair(best_worstcase_partition, sol);
        cout << "optimal solution: " << best_worstcase_solution << endl;
        return final_solution;
    }
    catch (GRBException e) {
        cout << "Gurobi error number [EnumSolve]: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Non Gurobi error during construction of static robust model object" << endl;
    }

    vector<vector<int> > vec1;
    vector<vector<double> > vec2;
    auto dummy_solution = make_pair(vec1, vec2);
    return dummy_solution;
}
