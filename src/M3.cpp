#include "../inc/M3.h"

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
    // Print arc summary of a graph, with costs if called from M3 Instance (is_costs)
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

void AdaptiveInstance::writeCosts() {
    // Write a costs file for the instance
    // Costs file has (p+1) lines - p sets of arc costs plus interdiction costs
    string filename = directory + name + "-costs_" + to_string(scenarios) + ".csv";
    ofstream myfile(filename);

    for (int q=0; q<scenarios; ++q) {
        stringstream ss;
        for (int a=0; a<arcs; ++a) {
            if (a!=0) {
                ss << ",";
            }
            ss << arc_costs[q][a];
        }

        string line = ss.str() + "\n";
        myfile << line;
    }

    stringstream ss;
    for (int a=0; a<arcs; ++a) {
        if (a!=0) {
            ss << ",";
        }
        ss << interdiction_costs[a];
    }
    string line = ss.str() + "\n";
    myfile << line;
}

void AdaptiveInstance::generateCosts(int interdiction, int min, int max) {
    /* 
     * Randomly Generate and set cost structure:
     * interdiction will populate interdiction_costs for every arc
     * min and max will bound uniform distribution to generate costs
     */
    interdiction_costs = vector<int>(arcs, interdiction);

    random_device rd;                           // obtain a random number from hardware
    mt19937 gen(rd());                          // seed the generator
    uniform_int_distribution<> distr(min, max); // define the range

    arc_costs = vector<vector<int> >(scenarios, vector<int>(arcs));

    for (int q=0; q<scenarios; ++q) {
        for (int a=0; a<arcs; ++a) {
            arc_costs[q][a] = distr(gen);
        }
    }
}

void AdaptiveInstance::readCosts() {
    // Read arc and interdiction costs from a file
    string line, word;
    string filename = directory + name + "-costs_" + to_string(scenarios) + ".csv";
    ifstream myfile(filename);
    int q = 0;
    vector<int> v;
    int cost;

    if (myfile.is_open())
    {
        while(getline(myfile, line)) {

            v.clear();
            stringstream str(line);

            while(getline(str, word, ',')) {
                cost = stoi(word);
                v.push_back(cost);
            }

            if (q < scenarios) {
                arc_costs.push_back(v);
            }

            else {
                interdiction_costs = v;
            }
            
            ++q;
        }
    }
}

void AdaptiveInstance::initCosts(int interdiction, int min, int max) {
    // If interdiction is not passed (<0), then read from file
    if (interdiction < 0) {
        readCosts();
    }
    else {
        generateCosts(interdiction, min, max);
        writeCosts();
    }
}

void AdaptiveInstance::printInstance(const LayerGraph& G) const {
    // Print Summary of Problem Instance
    cout << "k: " << policies << endl;
    cout << "p: " << scenarios << endl;
    G.printGraph(arc_costs, interdiction_costs, true);
}

vector<int> AdaptiveInstance::dijkstra(int q, const LayerGraph& G)
{
    // Compute shortest path 0-n-1
    vector<int> bar_S(nodes), dist(nodes), pred(nodes), result(arcs+1);
    vector<int> S;
    int min, node, j_node, arc, final_cost;

    for (int i=1; i<nodes; ++i) {
        // d(i) = i
        // pred(s) = -1
        bar_S[i] = i;
        dist[i] = interdiction_costs[0]*10;
        pred[i] = -1;
    }

    while (S.size()<nodes) {
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

    final_cost=dist[nodes-1];
    // cout << "final " << final_cost << endl;
    result[0]=final_cost;
    j_node=nodes-1;

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

void AdaptiveInstance::applyInterdiction(vector<float>& x_bar, bool rev){
    // Receive interdiction policy and update costs for the M3 instance
    // If rev, we are "removing" the interdiction policy and returning the instance to its original state
    for (int a=0; a<arcs; ++a){
        if (x_bar[a]==1) {
            for (int q=0; q<scenarios; ++q){
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

float AdaptiveInstance::validatePolicy(vector<float>& x_bar, const LayerGraph& G)
{
    // Solve shortest path on interdicted graph - check objectives match
    float objective=100000000;
    vector<int> sp_result;

    // Update M3 based on x_bar
    applyInterdiction(x_bar);

    // run dijstra on graph to get objective 
    for (int q=0; q<scenarios; ++q){
        sp_result = dijkstra(q, G);
        
        if (sp_result[0]<objective){
            objective=sp_result[0];
        }
    }
    
    // Revert M3 
    applyInterdiction(x_bar, true);

    return objective;
}

// ------ MIP Formulations for M3 ------
void SetPartitioningModel::configureModel(const LayerGraph& G, AdaptiveInstance& m3) {
    try
    {
        // ------ Initialize model and environment ------
        m3_env = new GRBEnv();
        m3_model = new GRBModel(*m3_env);
        m3_model->getEnv().set(GRB_DoubleParam_TimeLimit, 3600); // set time limit to 1 hour
        m3_model->set(GRB_IntParam_OutputFlag, 0);

        // ------ Decision variables ------
        vector<GRBVar> new_vector;

        // set partitioning variables
        for (int w = 0; w < policies; ++w){
            // cout << "in model constructor" << endl;
            h_matrix.push_back(new_vector);
            for (int q =0; q<scenarios; ++q){
                string varname = "H_" + to_string(w) + "_" + to_string(q);
                h_matrix[w].push_back(m3_model->addVar(0, 1, 0, GRB_BINARY, varname));
            }
        }

        // interdiction policy on arcs 'x' - REMEMBER to check that this initialization works
        for (int w = 0; w < policies; ++w) {
            x.push_back(new_vector);

            for (int a = 0; a < arcs; a++)
            {
                string varname = "x_" + to_string(w) + "_" + to_string(a);
                x[w].push_back(m3_model->addVar(0, 1, 0, GRB_BINARY, varname));
            }
        } 
       

        // objective func dummy 'z'
        z = m3_model->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS);

        vector<vector<GRBVar> > newnew_vector;
        
        // post interdiction flow 'pi'
        for (int w = 0; w<policies; ++w){
            pi.push_back(newnew_vector);

            for (int q = 0; q < scenarios; q++)
            {
                pi[w].push_back(new_vector);

                for (int i = 0; i < nodes; i++)
                {
                    string varname = "pi_" + to_string(w) + "_" + to_string(q) + "_" + to_string(i);
                    pi[w][q].push_back(m3_model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
                }
            }
        }

        // arc variable 'lambda'
        for (int w = 0; w<policies; ++w){
            lambda.push_back(newnew_vector);

            for (int q = 0; q < scenarios; q++)
            {
                lambda[w].push_back(new_vector);
                
                for (int a = 0; a < arcs; a++)
                {
                    string varname = "lambda_" + to_string(w) + "_"  + to_string(q) + "_" + to_string(a);
                    lambda[w][q].push_back(m3_model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));
                }
            }
        }
        // ------ Constraints ------
        // budget constraint
        for (int w = 0; w<policies; ++w){
            GRBLinExpr linexpr = 0;
            for (int a = 0; a < arcs; a++)
            {
                linexpr += x[w][a];
            }
            m3_model->addConstr(linexpr <= budget);
        }

        // z constraints
        for (int w = 0; w<policies; ++w){
            for (int q = 0; q < scenarios; q++)
            {
                GRBLinExpr linexpr = 0;
                linexpr += (pi[w][q][nodes - 1] - pi[w][q][0]); // b^\top pi (our b is simply a source-sink single unit of flow)
                for (int a = 0; a < arcs; ++a)
                {
                    linexpr += -lambda[w][q][a]; // u^\top \cdot lambda (our u is 1)
                }

                linexpr += M * (1 - h_matrix[w][q]);
                m3_model->addConstr(z <= linexpr);
            }
        }

        // main constraint for each arc
        int i;
        int j; // looping index 
        int jn; // node j
        int a;

        for (int w = 0; w<policies; ++w){
            for (int q=0; q<scenarios; ++q){
                for (i=0; i<nodes; ++i){
                    for (j=0; j<G.adjacency_list[i].size(); ++j){
                        jn=G.adjacency_list[i][j];
                        a=G.arc_index_hash[i][j];

                        // add constraint
                        m3_model->addConstr((pi[w][q][jn] - pi[w][q][i] - lambda[w][q][a]) <= m3.arc_costs[q][a] + (m3.interdiction_costs[a] * x[w][a]));
                    }
                }
            }
        }

        // for (int q = 0; q < p; ++q)
        // {
        //     linexpr = 0;
        //     for (int a = 0; a < m; ++a)
        //     {
        //         i = G.arcs[a].i;
        //         j = G.arcs[a].j;
        //         m3_model->addConstr((pi[q][j] - pi[q][i] - lambda[q][a]) <= M2Instance->arc_costs[q][a] + (M2Instance->interdiction_costs[a] * x[a]));
        //     }
        // }

        // set-partitioning constraint
        for (int q=0; q<scenarios; ++q){
            GRBLinExpr linexpr=0;

            for (int w=0; w<policies; ++w){
                linexpr+=h_matrix[w][q];
            }

            m3_model->addConstr(linexpr == 1);
        }

        // pi[0] = 0
        for (int w=0; w<policies; ++w) {
            for (int q = 0; q < scenarios; q++)
            {
                m3_model->addConstr(pi[w][q][0] == 0);
            }
        }
        m3_model->update();
    }
    catch (GRBException e)
    {
        cout << "Gurobi error number [m3_model, constructor]: " << e.getErrorCode() << "\n";
        cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        cout << "Non-gurobi error during optimization [m3_model]"
                  << "\n";
    }
}

vector<vector<float> > SetPartitioningModel::solve()
{
    /*
     * Returns vector of interdiction policy [0-m] (for every (w) for M3, plus an extra vector of size 2 as first element - [objective value, MIPGap])
     * Note: the return values are stored in class variables so unassigned in some compexp runs 
     */

    try
    {
        m3_model->optimize();

        try {
            vector<float> objective_vec;
            objective_vec.push_back(m3_model->get(GRB_DoubleAttr_ObjVal));
            objective_vec.push_back(m3_model->get(GRB_DoubleAttr_MIPGap));
            x_prime.push_back(objective_vec);

            for (int w=0; w<policies; ++w){
                vector<float> x_vector(arcs);
                for (int a = 0; a < arcs; a++)
                {
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
                vector<float> vec; 
                vec.push_back(-GRB_INFINITY);
                x_prime.push_back(vec);
            }

            else {
                cout << "Gurobi error number [M2Linear.optimize()]: " << e.getErrorCode() << "\n";
                cout << e.getMessage() << "\n";
            }
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

// // ------ Bender's Schemes for M2 ------
// BendersSub::BendersSub(){n=0; m=0; p=0;}
// 
// BendersSub::BendersSub(AdaptiveInstance *the_M2Instance)
// {
//     // ------ Initialize Basic Parameters ------
//     n = the_M2Instance->n;
//     m = the_M2Instance->m;
//     p = the_M2Instance->p;
// 
//     // ------ Initialize d and c costs ------
//     for (int a = 0; a < m; a++)
//     {
//         d.push_back(the_M2Instance->interdiction_costs[a]);
//     }
// 
//     for (int q = 0; q < p; q++)
//     {
//         c.push_back(the_M2Instance->arc_costs[q]);
//         c_bar.push_back(the_M2Instance->arc_costs[q]);
//     }
// 
//     // ------ Initialize Environment and Model ------
//     for (int q = 0; q < p; ++q)
//     {
//         Subenvs.push_back(new GRBEnv());
//         Submodels.push_back(new GRBModel(*Subenvs[q]));
//         Submodels[q]->set(GRB_IntParam_OutputFlag, 0);
//     }
// 
//     // ------ Decision Variables ------
//     for (int q = 0; q < p; ++q)
//     {
//         varname = "zeta_sub_" + to_string(q);
//         zeta_subs.push_back(Submodels[q]->addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS, varname));
//     }
//     for (int q = 0; q < p; ++q)
//     {
//         y_dummy = {};
//         y.push_back(y_dummy);
//         for (int a = 0; a < m; ++a)
//         {
//             varname = "y_" + to_string(q) + "_" + to_string(a);
//             y[q].push_back(Submodels[q]->addVar(0, 1, 0, GRB_CONTINUOUS, varname));
//         }
//     }
// 
//     // ------ Constraints ------
//     // constraints to bound objective value over q
//     obj_constr = new GRBConstr[p];
//     for (int q = 0; q < p; ++q)
//     {
//         linexpr = 0;
//         for (int a = 0; a < m; ++a)
//         {
//             linexpr += c_bar[q][a] * y[q][a];
//         }
//         obj_constr[q] = Submodels[q]->addConstr(zeta_subs[q] >= linexpr);
//     }
// 
//     // flow constraints
//     for (int q = 0; q < p; ++q)
//     {
//         for (int i = 0; i < n; ++i)
//         {
//             linexpr = 0;
//             if (i == 0)
//             {
//                 rhs = 1;
//             }
//             else if (i == n - 1)
//             {
//                 rhs = -1;
//             }
//             else
//             {
//                 rhs = 0;
//             }
//             for (int a = 0; a < m; a++)
//             {
//                 if (the_M2Instance->G.arcs[a].i == i)
//                 {
//                     // this arc is an outgoing arc for i
//                     linexpr += (1) * y[q][a];
//                 }
//                 else if (the_M2Instance->G.arcs[a].j == i)
//                 {
//                     // this arc is an outgoing arc for i
//                     linexpr += (-1) * y[q][a];
//                 }
//                 else
//                 {
//                     // this arc does not include i as an endpoint
//                     continue;
//                 }
//             }
//             Submodels[q]->addConstr(linexpr == rhs);
//         }
//         Submodels[q]->update();
//     }
// }
// 
// void BendersSub::update(vector<int> &xhat)
// {
//     // cout << "\nsubmodel: updating based on xbar, new interdiction policy: \n";
//     // for (int a = 0; a < m; ++a)
//     // {
//     //     cout << "xbar_" << a << ": " << xhat[a] << "\n";
//     // }
// 
//     // update array of constraints instead of only the parameter vector
//     for (int q = 0; q < p; ++q)
//     {
//         for (int a = 0; a < m; ++a)
//         {
//             // model.chgCoeff(obj_constr[q], variable, new cost)
//             if (xhat[a] > 0.5)
//             {
//                 Submodels[q]->chgCoeff(obj_constr[q], y[q][a], -(c[q][a] + d[a]));
//             }
//             else
//             {
//                 Submodels[q]->chgCoeff(obj_constr[q], y[q][a], -(c[q][a]));
//             }
//         }
//         Submodels[q]->update();
//     }
// }
// 
// vector<vector<float> > BendersSub::solve(int counter)
// {
//     try
//     {
//         // yhat has q elements - each one has m+1 elements (first is the objective, rest is the flow)
//         vector<vector<float> > yhat;
// 
//         for (int q = 0; q < p; ++q)
//         {
//             y_dummy2 = {};
//             yhat.push_back(y_dummy2);
//             // string modelname = "submodel_" + to_string(counter) + "q=" + to_string(q) + ".lp";
//             // Submodels[q]->write(modelname);
//             Submodels[q]->optimize();
// 
//             yhat[q].push_back(Submodels[q]->get(GRB_DoubleAttr_ObjVal));
//             // cout << "\nq = " + to_string(q);
//             // cout << "\nsubmodel obj: " << yhat[q][0];
//             // cout << "\narc values: \n";
//             for (int a = 0; a < m; ++a)
//             {
//                 yhat[q].push_back(y[q][a].get(GRB_DoubleAttr_X));
//                 // cout << "y_" << q << "_" << a << ": " << y[q][a].get(GRB_DoubleAttr_X) << "\n";
//             }
//         }
//         return yhat;
//     }
//     catch (GRBException e)
//     {
//         cout << "Gurobi error number [BendersSub::solve]: " << e.getErrorCode() << "\n";
//         cout << e.getMessage() << "\n";
//     }
//     catch (...)
//     {
//         cout << "Non-gurobi error during optimization [BendersSub::Solve]"
//                   << "\n";
//     }
// }
// 
// BendersSeparation::BendersSeparation()
// {
//     n = 0;
//     m = 0;
//     p = 0;
// }
// 
// BendersSeparation::BendersSeparation(GRBVar &the_zetabar, vector<GRBVar> &the_xbar, AdaptiveInstance *the_M2Instance)
// {
//     try
//     {
//         // ------ Initialize Basic Parameters ------
//         n = the_M2Instance->n;
//         m = the_M2Instance->m;
//         p = the_M2Instance->p;
// 
//         // ------ Initialize Submodel ------
//         subproblem = BendersSub(the_M2Instance);
// 
//         // ------ Initialize d and c costs ------
//         for (int a = 0; a < m; a++)
//         {
//             d.push_back(the_M2Instance->interdiction_costs[a]);
//         }
// 
//         for (int q = 0; q < p; q++)
//         {
//             c.push_back(the_M2Instance->arc_costs[q]);
//         }
// 
//         // ------ Initialize Variable containers ------
//         xbar = the_xbar;
//         zetabar = the_zetabar;
//         xprime.push_back(0);
// 
//         for (int q = 0; q < p; ++q)
//         {
//             if (q == 0)
//             {
//                 for (int a = 0; a < m; ++a)
//                 {
//                     xhat.push_back(0);
//                     xprime.push_back(0);
//                 }
//             }
//         }
//     }
//     catch (GRBException e)
//     {
//         cout << "Gurobi error number [BendersSeparation, constructor]: " << e.getErrorCode() << "\n";
//         cout << e.getMessage() << "\n";
//     }
//     catch (...)
//     {
//         cout << "Non-gurobi error during optimization [BendersSeparation]"
//                   << "\n";
//     }
// }
// 
// void BendersSeparation::callback()
// {
//     if (where == GRB_CB_MIPSOL)
//     {
//         // zeta_u = -getDoubleInfo(GRB_CB_MIPSOL_OBJBST); // ??????? 
// 
//         if (zeta_u - zeta_l >= epsilon)
//         {
//             cout << "in callback" << endl;
//             counter++;
//             // xhat = current solution from master problem, then update subproblem
//             for (int a = 0; a < m; ++a)
//             {
//                 xhat[a] = getSolution(xbar[a]);
//             }
//             subproblem.update(xhat);
//             // zeta_u = -getDoubleInfo(GRB_CB_MIPSOL_OBJ); // ??????? 
// 
//             // cout << "\n\n\n\nsolving sub from callback: \n\n\n\n";
//             yhat = subproblem.solve(counter);
//             zeta_temp = GRB_INFINITY;
// 
//             for (int q = 0; q < p; ++q)
//             {
//                 cout << "sub q=" << q << ": " << yhat[q][0] << endl;
//                 if (zeta_temp > yhat[q][0])
//                 {
//                     zeta_temp = yhat[q][0]; // first element of yhat[q] is the objective
//                 }
//             }
//             cout << "zeta_temp: " << zeta_temp << endl;
// 
//             // use yhat[(q-l)+1][1-m] to create new cut from LinExpr
//             for (int q = 0; q < p; ++q)
//             {
//                 new_cut = 0;
//                 for (int a = 0; a < m; ++a)
//                 {
//                     new_cut += (c[q][a] + d[a] * xbar[a]) * yhat[q][a + 1];
//                 }
//                 // add lazy cut to main model
//                 try
//                 {
//                     addLazy(zetabar <= new_cut);
//                     // cout << "\nadded cut: "
//                          // << "zetabar"
//                          // << "<=" << new_cut << "\n";
//                     ++cut_count;
//                 }
//                 catch (GRBException e)
//                 {
//                     cout << "Gurobi error number [BendersSeparation, addLazy]: " << e.getErrorCode() << "\n";
//                     cout << e.getMessage() << "\n";
//                 }
//             }
// 
//             if (zeta_l < zeta_temp)
//             {
//                 zeta_l = zeta_temp;
//             }
//             cout << "zeta_l: " << zeta_l << endl;
//             zeta_u = -getDoubleInfo(GRB_CB_MIPSOL_OBJBND); // ??????? 
//             cout << "zeta_u: " << zeta_u << endl;
//         }
//     }
// }
// 
// M2Benders::M2Benders(){int n, m=0;}
// 
// M2Benders::M2Benders(AdaptiveInstance *the_M2Instance)
// {
//     try
//     {
//         // ------ Assign Instance ------
//         M2Instance = the_M2Instance;
// 
//         // ------ Initialize model and environment ------
//         M2Bendersenv = new GRBEnv();
//         M2Bendersmodel = new GRBModel(*M2Bendersenv);
// 
//         M2Bendersmodel->getEnv().set(GRB_IntParam_LazyConstraints, 1);
//         M2Bendersmodel->set(GRB_IntParam_OutputFlag, 0);
//         M2Bendersmodel->set(GRB_DoubleParam_TimeLimit, 3600);
// 
//         // ------ Variables and int parameters ------
//         n = M2Instance->n;
//         m = M2Instance->m;
//         p = M2Instance->p;
//         r_0 = M2Instance->r_0;
//         instance_name = M2Instance->instance_name;
//         setname = M2Instance->setname;
// 
//         // ------ Decision variables ------
//         string varname;
// 
//         // interdiction variable 'x'
//         for (int a = 0; a < m; a++)
//         {
//             varname = "x_" + to_string(a);
//             x.push_back(M2Bendersmodel->addVar(0, 1, 0, GRB_BINARY, varname));
//         }
// 
//         // objective function variable 'zeta'
//         varname = "zeta";
//         zeta = M2Bendersmodel->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS, varname);
// 
//         // ------ Initialize separation/callback object ------
//         sep = BendersSeparation(zeta, x, M2Instance);
// 
//         // ------ Only Budget Constraints Initially! ------
//         linexpr = 0;
//         for (int a = 0; a < m; a++)
//         {
//             linexpr += x[a];
//         }
//         M2Bendersmodel->addConstr(linexpr <= r_0);
// 
//         // ------ Trying to add the first lazy contraint ------
//         // cout << "\n\n\n\nsolving sub from constructor: \n\n\n\n";
//         sep.yhat = sep.subproblem.solve(0);
//         for (int q = 0; q < p; ++q)
//         {
//             linexpr = 0;
//             for (int a = 0; a < m; ++a)
//             {
//                 linexpr += (sep.c[q][a] + (sep.d[a] * x[a])) * sep.yhat[q][a + 1];
//             }
//             M2Bendersmodel->addConstr(zeta <= linexpr);
//         }
// 
//         modelname = "modelfiles/" + setname + "/" + instance_name + "-M2BendersModel.lp";
//         M2Bendersmodel->write(modelname);
// 
//     }
//     catch (GRBException e)
//     {
//         cout << "Gurobi error number [M2Benders, constructor]: " << e.getErrorCode() << "\n";
//         cout << e.getMessage() << "\n";
//     }
//     catch (...)
//     {
//         cout << "Non-gurobi error during optimization [BendersSPSub]"
//                   << "\n";
//     }
// }
// 
// vector<float> M2Benders::solve()
// {
//     /*
//      * Return vector returns the objective value [0] and interdiction policy [1-m+1]
//      * Again for computational experiments this can be ignored and not assigned as run stats are class vars
//      */ 
// 
//     // ------ Set Callback on Master Model
//     M2Bendersmodel->setCallback(&sep);
// 
// 
//     // ------ Optimize Inside Benders Scheme -------
// 
//     try
//     {
//         clock_t model_begin = clock();
//         M2Bendersmodel->optimize();
//         running_time = float(clock() - model_begin) / CLOCKS_PER_SEC;
// 
//         try {
//             optimality_gap = M2Bendersmodel->get(GRB_DoubleAttr_MIPGap);
//             sep.xprime[0] = M2Bendersmodel->get(GRB_DoubleAttr_ObjVal);
// 
//             for (int a = 0; a < m; ++a)
//             {
//                 sep.xprime[a + 1] = x[a].get(GRB_DoubleAttr_X);
//             }
// 
//         }
// 
//         catch (GRBException e) {
// 
//             if (e.getErrorCode() == 10005 || e.getErrorCode() == 10013) {
//                 // in this case the error is because the model was unbounded
//                 // set the optimality gap to -1 and we'll list as unbounded 
//                 // put objective value as -infinity (objectives are negated min -z)
//                 // set arc interdiction values as -1 - gurobi won't give us a solution one unboundedness is proven
//                 optimality_gap = -1;
//                 sep.xprime[0] = -GRB_INFINITY;
// 
//                 for (int a = 0; a < m; ++a)
//                 {
//                     sep.xprime[a + 1] = -1;
//                 }
//             }
// 
//             else {
//                 cout << "Gurobi error number [M2Benders.optimize()]: " << e.getErrorCode() << "\n";
//                 cout << e.getMessage() << "\n";
//             }
// 
//         }
//         
//         cut_count = sep.cut_count;
//     }
// 
//     catch (GRBException e)
//     {
//         cout << "Gurobi error number [M2Benders.optimize()]: " << e.getErrorCode() << "\n";
//         cout << e.getMessage() << "\n";
//     }
//     catch (...)
//     {
//         cout << "Non-gurobi error during optimization [M2Benders.optimize()]"
//                   << "\n";
//     }
// 
// 
//     // for submodel testing and stuff
//     // vector<int> test_xhat = {1, 1, 0};
// 
//     // sep.subproblem.Submodel->write("spmodel.lp");
//     // sep.yhat = sep.subproblem.solve(0);
//     // sep.subproblem.update(test_xhat);
//     // sep.subproblem.Submodel->write("spmodelupdated.lp");
// 
//     delete sep.subproblem.obj_constr;
//     for (int q = 0; q < p; ++q)
//     {
//         delete sep.subproblem.Subenvs[q];
//         delete sep.subproblem.Submodels[q];
//     }
//     return sep.xprime;
// }

vector<int> initKappa(int p, int k) {
    // initialize partition vector based on total number in set (p) and exact number of partitions required
    vector<int> kappa(p, 0);

    for (int q=(p-k+1); q<p; ++q){
        kappa[q]=(q-(p-k));
    }

    return kappa;
}


void RobustAlgoModel::configureModel(const LayerGraph& G, AdaptiveInstance& m3) {
    // Construct baseline model with no constraints
    algo_env = new GRBEnv();
    algo_model = new GRBModel(*algo_env);
    algo_model->set(GRB_IntParam_OutputFlag, 0);

    // Decision Variables
    string varname = "z";
    z = algo_model->addVar(0, GRB_INFINITY, -1, GRB_CONTINUOUS, varname); // objective function dummy variable
    for (int a=0; a<arcs; ++a) {varname = "x_" + to_string(a); x.push_back(algo_model->addVar(0, 1, 0, GRB_BINARY, varname));} // interdiction policy on arcs

    vector<GRBVar> tempvector;
    for (int q=0; q<scenarios; ++q) {
        pi.push_back(tempvector); // post interdiction s-i best path (for every q)
        for (int i=0; i<nodes; ++i) {varname = "pi_" + to_string(q) + "_" + to_string(i); pi[q].push_back(algo_model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));}
    }

    for (int q=0; q<scenarios; ++q) {
        lambda.push_back(tempvector); // lambda variable on arcs (for every q)
        for (int a=0; a<arcs; ++a) {varname = "lambda_" + to_string(q) + "_" + to_string(a); lambda[q].push_back(algo_model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, varname));}
    }

    // Add budget Constraint
    GRBLinExpr linexpr = 0;
    for (int a = 0; a < arcs; a++)
    {
        linexpr += x[a];
    }
    algo_model->addConstr(linexpr <= budget, "budget");
    algo_model->update();

    // Populate Global Constraints
    // Arc/Dual Constraints
    for (int q=0; q<scenarios; ++q){
        dual_constraints.push_back(vector<GRBTempConstr>());

        for (int i=0; i<nodes; ++i) {
            for (int j=0; j<G.adjacency_list[i].size(); ++j){
                int next = G.adjacency_list[i][j];
                int a = G.arc_index_hash[i][j];

                GRBTempConstr constraint = pi[q][next] - pi[q][i] - lambda[q][a] <= m3.arc_costs[q][a] + (m3.interdiction_costs[a]*x[a]);
                dual_constraints[q].push_back(constraint);
            }
        }
    }

    // Objective Constraints
    for (int q=0; q<scenarios; ++q) {
        GRBLinExpr linexpr = 0;
        linexpr += (pi[q][nodes-1] - pi[q][0]); // b^\top pi

        for (int a=0; a<arcs; ++a){
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
        algo_model->addConstr(pi[q][0]==0, zero_name);
        algo_model->addConstr(z_constraints[q], z_name);

        int a = 0;
        for (GRBTempConstr constraint : dual_constraints[q]) {
            string dual_name = "dual_" + to_string(q) + "_" + to_string(a);
            algo_model->addConstr(constraint, dual_name);
            ++a;
        }
    }
    
    algo_model->update();
}

void RobustAlgoModel::reverse_update(vector<int>& subset) {
    // remove all constraints from model - i.e. all constraints associated with this subset
    for (int q : subset) {
        string zero_name = "zero_" + to_string(q);
        string z_name = "z_" + to_string(q);

        GRBConstr zero_constraint = algo_model->getConstrByName(zero_name);
        GRBConstr z_constraint = algo_model->getConstrByName(z_name);
        algo_model->remove(zero_constraint);
        algo_model->remove(z_constraint);

        int a = 0;
        for (GRBTempConstr constraint : dual_constraints[q]) {
            string dual_name = "dual_" + to_string(q) + "_" + to_string(a);
            GRBConstr dual_constraint = algo_model->getConstrByName(dual_name);
            algo_model->remove(dual_constraint);
            ++a;
        }
    }

    algo_model->update();
}

vector<double> RobustAlgoModel::solve(int counter) {
    // solve static model, return policy and objective value
    string lp_filename = "static_model_" + to_string(counter) + ".lp";
    algo_model->write(lp_filename);
    algo_model->optimize();
    vector<double> sol(arcs+1, 0);
    
    sol[0] = algo_model->get(GRB_DoubleAttr_ObjVal);
    sol[0] = -sol[0];
    for (int a=1; a<arcs+1; ++a) {
        sol[a] = x[a-1].get(GRB_DoubleAttr_X);
    }

    algo_model->reset();

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
            
            ++kappa[q]; max[q] = max_int(max[q], kappa[q]);

            for (int u=q+1; u<=(p-(k-max[q])); ++u){
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

pair<vector<vector<int> >, vector<vector<double> > > enumSolve(AdaptiveInstance& m3, const LayerGraph& G){
    // Pass a AdaptiveInstance and solve using enumeration algorithm
    // Initialize and maintain an H matrix representing partitioning 
    // Enumerate all possible partitions by enumerating H matrix
    // Enumeration is done through the 'string' method, and then the H matrix/problem instance are updated (??)
    // Return: pair of nested vectors:
    //      first: ints - optimal partition - obj value singleton at pos 0
    //      second: doubles - k interdiction policies

    // ints and bool
    bool next = true;
    int p = m3.scenarios;
    int k = m3.policies;
    int m = m3.arcs;
    
    // initialize partitioning 'string' vector as in Orlov Paper 
    // initialize corresponding max vector 
    vector<int> kappa = initKappa(p, k);
    vector<int> max = kappa;

    // initialize solution vector
    vector<double> arc_vec(m, 0);
    vector<vector<double> > sol(k, arc_vec);

    // initialize static robust model
    try {
        RobustAlgoModel static_robust = RobustAlgoModel(m3);
        static_robust.configureModel(G, m3);


        // objective value maintained here (and optimal partition)
        double best_worstcase_solution = 0;
        vector<vector<int> > best_worstcase_partition; 

        int counter = 0;

        // enumerate while not 'failing' to get next partition
        while (next) {
            vector<vector<double> > temp_sol(k, arc_vec);
            double temp_worst_sol = DBL_MAX;
            vector<vector<int> > partition = kappa_to_partition(kappa, k, p);

            for (int w=0; w<k; ++w) {
                // for every subset in partition, solve M2 for k=1
                static_robust.update(partition[w]);
                vector<double> temp_single_solution = static_robust.solve(counter);
                static_robust.reverse_update(partition[w]);
                ++counter;

                // update temp_worst_sol if this subset is worse
                if (temp_single_solution[0] < temp_worst_sol) {temp_worst_sol = temp_single_solution[0];}
            }

            if (temp_worst_sol > best_worstcase_solution) {
                best_worstcase_solution = temp_worst_sol; 
                sol = temp_sol;
                best_worstcase_partition = partition;
            }

            next = nextKappa(kappa, max, k, p);
        }
        

        vector<int> final_obj(1, best_worstcase_solution);
        best_worstcase_partition.insert(best_worstcase_partition.begin(), final_obj);
        auto final_solution = make_pair(best_worstcase_partition, sol);
        return final_solution;
    }
    catch (GRBException e) {
        cout << "Gurobi error number [EnumSolve]: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Non Gurobi error during construction of static robust model object" << endl;
    }

    vector<vector<int>> vec1;
    vector<vector<double>> vec2;
    auto dummy_solution = make_pair(vec1, vec2);
    return dummy_solution;
}

vector<vector<double> > extendByOne(pair<vector<vector<int> >, vector<vector<double> > >& k_solution, AdaptiveInstance& m3) {
    // Use an optimal solution found by the enumerative algorithm for k, to find a good solution for k+1
    // Take the worst subset in the optimal partition and "split it in 2"
    // "Split it in two": solve that subset for k = 2
    
    // int values
    int k = m3.policies;
    
    // find worst subset in optimal partition
    double min_subset_obj = GRB_INFINITY;
    int min_subset_windex;

    for (int w=0; w<k; ++w){
        if (k_solution.second[k][0] < min_subset_obj) {
            min_subset_obj = k_solution.second[k][0];
            min_subset_windex = w;
        }
    }

    int p_prime = k_solution.first[min_subset_windex].size(); // the new p - number of scenarios in the subset that we will work on

    if (p_prime == 1) {
        return k_solution.second; 
    }
    else if (p_prime == 2) {
        // create an m3Instance copy with 1 of the scenarios
        cout << "placeholder" << endl;

        // create an m3Instance copy with the other
        
        // solve each of them for k=1, p=1
    }
    else {
        // create an m3Instance copy with all scenarios in the subset
        // solve it for k=2
        cout << "placeholder" << endl;
    }
}

long getCurrentTime() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}
