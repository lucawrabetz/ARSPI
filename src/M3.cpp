#include "../inc/M3.h"

long getCurrentTime() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

Arc::Arc(){i=0; j=0;}

Arc::Arc(int the_i, int the_j){i=the_i; j=the_j;}

Graph::Graph(){n=0; m=0;}

Graph::Graph(const string &filename, int the_n, int the_k_0)
{
    // Graph constructor from file
    string line, word;
    ifstream myfile(filename);
    int max_sub = -2;

    if (myfile.is_open())
    {
        n = the_n;
        int arc_index = 0;
        int i;
        int j;
        // string sub;

        // cout << "n: " << n << endl;
        // cout << "filename: " << filename << endl;

        arc_index_hash = vector<vector<int> >(n);
        adjacency_list = vector<vector<int> >(n);
        n_n_adjacency_list = vector<vector<int> >(n, vector<int>(n, 0));

        while (getline(myfile, line))
        {
            stringstream str(line);
            int counter = 0;

            while(getline(str, word, ' ')){
                if (counter==0) {i = stoi(word);}
                else if (counter==1) {j = stoi(word);}
                // else if (counter==3) {
                //     sub = word;
                //     sub.pop_back();
                // }
                ++counter;
            }

            // int subval = stoi(sub)-1;
            
            // aif (subval > max_sub) {max_sub = subval;}

            arcs.push_back(Arc(i, j));
            // subgraph.push_back(subval);

            adjacency_list[i].push_back(j);
            arc_index_hash[i].push_back(arc_index);
            n_n_adjacency_list[i][j] = 1;

            ++arc_index; 
        }

        m = arcs.size();
        kbar = the_k_0;
    }
}

void Graph::printGraph(vector<vector<int> > costs, vector<int> interdiction_costs, bool is_costs) const
{
    // Print arc summary of a graph, with costs if called from M3 Instance (is_costs)
    cout << "n: " << n << ", m: " << m << endl;
    int p = costs.size();

    if (is_costs){
        for (int a = 0; a < m; a++)
        {
            cout << "(" << arcs[a].i << "," << arcs[a].j << ")";
            // cout << "   interdiction_cost: " << interdiction_costs[a];
            cout << " costs:";

            for (int q = 0; q<p; ++q){
                cout << " " << q << ": " << costs[q][a];
            }
            cout << endl;
        }
    }

    else {
        for (int a = 0; a < m; a++)
        {
            cout << "(" << arcs[a].i << "," << arcs[a].j << endl;
        }
    }
}

AdaptiveInstance::AdaptiveInstance(AdaptiveInstance* m3, vector<int>& keep_scenarios) {
    // Copy constructor only keeping a subset of the scenarios
    scenarios = keep_scenarios.size();
    policies = m3->policies;
    budget = m3->budget;
    nodes = m3->nodes;
    arcs = m3->arcs;
    interdiction_costs = m3->interdiction_costs;
    arc_costs = vector<vector<int> >(scenarios);
    scenario_index_map = vector<int>(scenarios);
    for (int q=0; q<scenarios; q++) {
        arc_costs[q] = m3->arc_costs[keep_scenarios[q]];
        scenario_index_map[q] = keep_scenarios[q];
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

void AdaptiveInstance::generateCosts(float interdiction, int a, int b, int dist, vector<int> subgraph) {
    /* 
     * Randomly Generate and set cost structure:
     * interdiction is the fractional multiplier on average cost used for (anti)-proportional interdiction costs
     * a and b are distribution parameters:
     *  a is the lowest mean
     *  b is the increase in every next mean
     */

    random_device rd;                           
    mt19937 gen(rd());
    
    arc_costs = vector<vector<int> >(scenarios, vector<int>(arcs));
    vector<int> scenario_sub_map(scenarios);

    int sub = 0;
    for (int q=0; q<scenarios; ++q) {
        if (sub==kbar) {sub=0;}
        scenario_sub_map[q] = sub;
        ++sub;
    }

    if (dist==3) {
        normal_distribution<> cheap(a, 20);
        normal_distribution<> expensive(b, 20);

        for (int a=0; a<arcs; ++a) {
            for (int q=0; q<scenarios; ++q) {
                if (scenario_sub_map[q] == subgraph[a]) {
                    // draw cheap
                    arc_costs[q][a] = abs(cheap(gen));
                }
                else {
                    // draw expensive
                    arc_costs[q][a] = abs(expensive(gen));
                }
            }
        }
    }

    if (dist==2) {
        // bernoulli_distribution binary(0.5);

        // define the distribution for q=0
        // normal_distribution<> norm_high(b, 50);
        // normal_distribution<> norm_low(a, 50);

        // define vector of p normal distributions
        // define random indexing map
        vector<normal_distribution<>> normals;
        vector<int> indices;
        uniform_int_distribution<> unif(1, kbar);

        for (int q=0; q<kbar; ++q) {
            normal_distribution<> norm(a + b*q, 50);
            normals.push_back(norm);
            indices.push_back(q);
        }
        for (int q=kbar; q < scenarios; ++q) {
            int index = unif(gen);
            indices.push_back(index);
        }


        for (int a=0; a<arcs; ++a) {

            // obtain seed and shuffle indices
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            shuffle (indices.begin(), indices.end(), default_random_engine(seed));

            for (int q=0; q<scenarios; ++q) {

                arc_costs[q][a] = abs(normals[indices[q]](gen));

                // if (ind) {
                //     if (q==0 || q==1) {
                //         arc_costs[q][a] = abs (norm_low(gen)); 
                //     } 
                //     else {
                //         arc_costs[q][a] = abs (norm_high(gen));
                //     }
                // }
                // else {
                //     if (q==0 || q==1) {
                //         arc_costs[q][a] = abs (norm_high(gen));
                //     } 
                //     else {
                //         arc_costs[q][a] = abs (norm_low(gen));
                //     }
                // }
            }
        }
    }

    else {
        for (int q=0; q<scenarios; ++q) {
            // define the distribution
            uniform_int_distribution<> unif(a+10, b+10); 
            normal_distribution<> norm(a+10, b);
            for (int a=0; a<arcs; ++a) {
                if (dist==0) {arc_costs[q][a] = unif(gen);}
                else if (dist==1) {arc_costs[q][a] = round(norm(gen));}
            }
        }
    }

    interdiction_costs = vector<int>(arcs, 10);
    // int total = 0;

    // for (int q=0; q<scenarios; ++q) {
    //     total += (a + q*b);
    // }

    // int interdiction_baseline = total / scenarios;

    // for (int a=0; a<arcs; ++a) {
    //     cout << "arc: " << a << endl;
    //     int average = 0;
    //     for (int q=0; q<scenarios; ++q) {
    //         average += arc_costs[q][a];
    //         cout << arc_costs[q][a] << " ";
    //     }
    //     average = average / scenarios;
    //     cout << "average: " << average << endl;
    //     cout << "baseline: " << interdiction_baseline << endl;
    //     cout << "difference: " << interdiction_baseline - average << endl;


    //     interdiction_costs[a] = (interdiction_baseline * interdiction) + (interdiction_baseline - average) * interdiction;
    //     cout << "final cost: " << interdiction_costs[a] << endl;
    // }

    // cout << "doing interdiction costs" << endl;
    // cout << "fraction: " << interdiction << endl;
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

void AdaptiveInstance::initCosts(float interdiction, int a, int b, int dist, const Graph &G, bool gen) {
    // If gen is false read from file
    // Distributions:
    //  - dist = 0: uniform, a = min, b = max
    //  - dist = 1: normal, a = mean, b = stddev
    if (gen == false) {
        readCosts();
    }
    else {
        generateCosts(interdiction, a, b, dist, G.subgraph);
        writeCosts();
    }
    // for (int i = 0; i < arcs; i++) {
    //     if (G.arcs[i].sub == -1) {
    //         arc_costs
    //     }
    // }
}

void AdaptiveInstance::printInstance(const Graph& G) const {
    // Print Summary of Problem Instance
    cout << "k: " << policies << endl;
    cout << "p: " << scenarios << endl;
    G.printGraph(arc_costs, interdiction_costs, true);
}

vector<int> AdaptiveInstance::dijkstra(int q, const Graph& G)
{
    // Compute shortest path 0-n-1
    vector<int> pred(nodes), result(arcs+1);
    std::priority_queue<pair<int, int>, vector<pair<int,int> >, std::greater<pair<int, int> > > dist;
    std::unordered_set<int> S, bar_S;

    bar_S.insert(0);
    dist.push(make_pair(0, 0));
    pred[0] = 0;
    for (int i=1; i<nodes; ++i) {
        bar_S.insert(i);
        pred[i] = -1;
    }

    while (S.size()<nodes) {
        int node = dist.top().second;
        int distance = dist.top().first;
        cout << dist.size() << endl;
        cout << "node: " << node << ", distance: " << distance << endl;
        dist.pop();
        cout << "distsize: " << dist.size() << endl;
        bar_S.erase(node); 
        S.insert(node);
        cout << "S size: " << S.size() << endl;
        cout << "bar S size: " << bar_S.size() << endl;


        for (int i=0; i<G.arc_index_hash[node].size(); ++i){
            int arc=G.arc_index_hash[node][i];
            int j_node=G.adjacency_list[node][i];
            cout << "i, j, arc" << node << ", " << j_node << ", " << arc << endl;
            dist.push(make_pair(distance+arc_costs[q][arc], j_node));
            pred[j_node]=node;
        }
        break;
    }

    return pred;

    result[0]=0;
    int j_node=nodes-1;

    while(j_node!=0){
        // find the arc index for the arc i, j (looping backwards through sp tree using pred to determine i)
        int node=pred[j_node];
        int arc;
        cout << "j_node: " << j_node << ", node: " << node << endl;

        for (int i=0; i<G.adjacency_list[node].size(); ++i){
            if (G.adjacency_list[node][i]==j_node){
                arc=G.arc_index_hash[node][i];
                break;
            }
        }

        result[0]+=arc_costs[q][arc];
        result[arc+1]=1;
        j_node=node;
    }

    return result;
}

void AdaptiveInstance::applyInterdiction(vector<double>& x_bar, bool rev){
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

float AdaptiveInstance::validatePolicy(vector<double>& x_bar, const Graph& G)
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
void SetPartitioningModel::configureModel(const Graph& G, AdaptiveInstance& m3) {
    try
    {
        // ------ Initialize solution to appropriate size -----
        current_solution = AdaptiveSolution(m3.policies, m3.scenarios);

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

        // constraint blocking interdiction of connecting arcs
        // for (int a=0; a<arcs; a++) {
        //     if (G.arcs[a].sub == -1) {
        //         for (int w=0; w<policies; w++) {
        //             m3_model->addConstr(x[w][a] == 0);
        //         }
        //     }
        // }

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

void SetPartitioningModel::solve()
{
    /*
     * Returns vector of interdiction policy [0-m] (for every (w) for M3, plus an extra vector of size 2 as first element - [objective value, MIPGap])
     * Note: the return values are stored in class variables so unassigned in some compexp runs 
     */

    try
    {
        long begin = getCurrentTime();
        m3_model->optimize();
        long time = getCurrentTime() - begin;
        current_solution.most_recent_solution_time = time;

        try {
            // vector<float> objective_vec;
            // objective_vec.push_back(m3_model->get(GRB_DoubleAttr_ObjVal));
            // objective_vec.push_back(m3_model->get(GRB_DoubleAttr_MIPGap));
            // x_prime.push_back(objective_vec);
            current_solution.worst_case_objective = m3_model->get(GRB_DoubleAttr_ObjVal);

            // cout << "policies: " << policies << endl;
            // cout << "scenarios: " << scenarios << endl;
            // cout << h_matrix.size() << endl;
            // cout << h_matrix[0].size() << endl;
            for (int w=0; w<policies; ++w){
                vector<double> x_vector;
                for (int a = 0; a < arcs; a++)
                {
                    x_vector.push_back(x[w][a].get(GRB_DoubleAttr_X));
                }
                // Policy.objective taken care of after in computeAllObjectives
                Policy temp_policy = Policy(arcs, x_vector, -1);
                current_solution.solutions[w] = temp_policy;

                // cout << "w: " << w << endl;
                // current_solution.partition[w].clear();
                for (int q=0; q<scenarios; q++) {
                    // cout << "q: " << q << endl;
                    if (h_matrix[w][q].get(GRB_DoubleAttr_X) > 0.5) {
                        current_solution.partition[w].push_back(q);
                    }
                }
                // x_prime.push_back(x_vector);
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
                cout << "Gurobi error number [SetPartitioningModel::solve]: " << e.getErrorCode() << "\n";
                cout << e.getMessage() << "\n";
            }
        }
    }
    catch (GRBException e)
    {
        cout << "Gurobi error number [SetPartitioningModel::solve]: " << e.getErrorCode() << "\n";
        cout << e.getMessage() << "\n";
    }
    catch (...)
    {
        cout << "Non-gurobi error during optimization [SetPartitioningModel::solve]"
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

void RobustAlgoModel::configureModel(const Graph& G, AdaptiveInstance& m3) {
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

    // constraint blocking interdiction of connecting arcs
    // for (int a=0; a<arcs; a++) {
    //     if (G.arcs[a].sub == -1) {
    //         algo_model->addConstr(x[a] == 0);
    //     }
    // }

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
        // cout << "q = " << q << endl;
        string zero_name = "zero_" + to_string(q);
        string z_name = "z_" + to_string(q);
        algo_model->addConstr(pi[q][0]==0, zero_name);
        algo_model->addConstr(z_constraints[q], z_name);

        int a = 0;
        // cout << "about to loop through dual constraints" << endl;
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
    // algo_model->write(lp_filename);
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

void printSolution(pair<vector<vector<int> >, vector<Policy> >& sol, string solname){
    // print adaptive solution
    cout << endl;
    if (solname == ""){cout << "solution:" << endl;}
    else {cout << solname << ":" << endl;} 
    int k = sol.first.size();

    for (int w=0; w<k; ++w){
        cout << "subset: ";
        for (int q : sol.first[w]){
            cout << q << " ";
        }
        cout << "- objective: " << sol.second[w].objective << endl;
    }
    cout << endl;
}

void AdaptiveSolution::computeAllObjectives(const Graph& G, AdaptiveInstance& m3) {
    // Populate / recompute the all_objectives vector based on interdiction policies and follower costs
    // Worst case objective for every policy will be assigned to solution[w].objective (in the Policy struct)
    RobustAlgoModel static_model = RobustAlgoModel(m3);
    static_model.configureModel(G, m3);
    double worst_objective = DBL_MAX;

    for (int w=0; w < policies; w++) {
        static_model.update(partition[w]);
        vector<double> sol = static_model.solve(0);
        solutions[w].objective = sol[0];
        if (sol[0] < worst_objective) {worst_objective = sol[0];}
        static_model.reverse_update(partition[w]);
    }
    worst_case_objective = worst_objective;
}

void AdaptiveSolution::logSolution(const Graph& G, AdaptiveInstance& m3, string title, bool policy) {
    int arcs = solutions[0].binary_policy.size();
    // print subsets and objectives
    cout << endl;
    if (title == ""){cout << "solution:" << endl;}
    else {cout << title << ":" << endl;} 
    cout << "Worst Case Objective: " << worst_case_objective << endl;
    int k = partition.size();

    for (int w=0; w<k; ++w){
        cout << "subset: { ";
        for (int q : partition[w]){
            cout << q << " ";
        }
        cout << "} - objective: " << solutions[w].objective << endl;
        if (policy) {
            cout << "interdicted arcs: ";
            for (int a = 0; a < arcs; a++) {
                if (solutions[w].binary_policy[a] > 0.5) {
                    cout << a << " (" << 
                        G.arcs[a].i << ", " <<
                        G.arcs[a].j << ")";
                        // G.arcs[a].sub << ") ";
                }
            }
            cout << endl << endl;
        }
    }
    cout << endl;
}

void AdaptiveSolution::mergeEnumSols(AdaptiveSolution sol2, AdaptiveInstance* instance2, int split_index) {
    // Merge 2 adaptive solutions
    // The AdaptiveSolution being passed has k=2
    // The AdaptiveSolution being worked on needs its split_index replaced:
    //   - in .partition, by the two subsets in sol2
    //   - in .solutions, by the two policies in sol2
    //   - in both cases, we can just replace split_index by the first, and add the 
    //   second to the end
    // Additionally, .policies must be increased by 1
    // Remember - change-of-scenario copy constructor reindexes, but we maintain
    // a map of the original indices (sol2 has this map)

    // policies++
    set_policies(policies+1);

    // reindexed copy of sol2.partition
    vector<vector<int> > reindexed_partition(sol2.policies);
    for (int w=0; w<sol2.policies; w++) {
        for (int i : sol2.partition[w]) {
            int q = instance2->scenario_index_map[i]; 
            reindexed_partition[w].push_back(q);
        }
    }

    // replace and push_back in partition
    partition[split_index] = reindexed_partition[0];
    for (int w=1; w<sol2.policies; w++) {
        partition.push_back(reindexed_partition[w]);
    }
    
    // replace in solutions
    solutions[split_index] = sol2.solutions[0];
    for (int w=1; w<sol2.policies; w++) {
        solutions.push_back(sol2.solutions[w]);
    }
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

pair<vector<vector<int> >, vector<Policy> > mergeEnumSols(pair<vector<vector<int> >, vector<Policy> >& sol1, pair<vector<vector<int> >, vector<Policy> >& sol2, int w_index){
    // Merge 2 solutions e.g. when extending an enum solve, combine the new split worst subset
    // w_index is the position of the split subset in the original solution (sol1)

    // place the first subset/policy of sol2 at position w_index of sol1, replacing the subset that was split
    sol1.first[w_index] = sol2.first[0];
    sol1.second[w_index] = sol2.second[0];

    // place the second subset/policy of sol2 at the end of sol1, extending it to k+1 
    sol1.first.push_back(sol2.first[1]);
    sol1.second.push_back(sol2.second[1]);

    return sol1;
}

AdaptiveSolution enumSolve(AdaptiveInstance& m3, const Graph& G){
    // Pass a AdaptiveInstance and solve using enumeration algorithm
    // Enumeration maintained as in Orlov paper
    // OLD Return: pair:
    //      first: vector of vector of ints - optimal partition (optimal objective is with policies)
    //      second: vector of Policy - k interdiction policies
    //  NEW Return: Adaptive Solution Struct (basic info, nested vector of partitions, vector of policies)
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
        long begin = getCurrentTime();
        RobustAlgoModel static_robust = RobustAlgoModel(m3);
        static_robust.configureModel(G, m3);

        // objective value maintained here (and optimal partition)
        double best_worstcase_objective = 0;
        vector<double> final_objectives(k, 0);
        vector<vector<int> > best_worstcase_partition; 

        int counter = 0;

        // enumerate while not 'failing' to get next partition
        while (next) {
            vector<vector<double> > temp_sol(k, arc_vec);
            double temp_worst_objective = DBL_MAX;
            vector<double> temp_objectives(k, 0);
            vector<vector<int> > partition = kappa_to_partition(kappa, k, p);

            for (int w=0; w<k; ++w) {
                // for every subset in partition, solve M2 for k=1
                static_robust.update(partition[w]);
                vector<double> temp_single_solution = static_robust.solve(counter);
                temp_objectives[w] = temp_single_solution[0];
                static_robust.reverse_update(partition[w]);
                ++counter;

                // update temp_worst_objective and temp_sol if this subset is worse
                if (temp_single_solution[0] < temp_worst_objective) {temp_worst_objective = temp_single_solution[0];}
            }

            if (temp_worst_objective > best_worstcase_objective) {
                best_worstcase_objective = temp_worst_objective; 
                sol = temp_sol;
                best_worstcase_partition = partition;
                final_objectives = temp_objectives;
            }

            next = nextKappa(kappa, max, k, p);
        }
        
        vector<Policy> final_policies(k, Policy(m));

        for (int w=0; w<k; ++w) {
            final_policies[w].set_policy(sol[w]);
            final_policies[w].set_objective(final_objectives[w]);
        }

        // vector<int> final_obj(1, best_worstcase_solution);
        // best_worstcase_partition.insert(best_worstcase_partition.begin(), final_obj);
        // OLD 
        // auto final_solution = make_pair(best_worstcase_partition, final_policies); 
        // NEW
        AdaptiveSolution final_solution = AdaptiveSolution(k, p, best_worstcase_partition, final_policies);
        final_solution.set_worst_case_objective(best_worstcase_objective);
        long time = getCurrentTime() - begin;
        final_solution.most_recent_solution_time = time;
        return final_solution;
    }
    catch (GRBException e) {
        cout << "Gurobi error number [EnumSolve]: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Non Gurobi error during construction of static robust model object" << endl;
    }

    AdaptiveSolution dummy_solution;
    return dummy_solution;
}

void AdaptiveSolution::extendByOne(AdaptiveInstance& m3, const Graph& G, bool mip_subroutine) {
    // Use an optimal solution found for k, to find a good solution for k+1
    // Take the worst subset in the optimal partition and "split it in 2"
    // "Split it in two": solve that subset for k = 2
    // Return format will always be exactly the same as enum solve, just for k greater than one
    // Naming - m3_prime is the copy of m3, and anything_prime is associated with m3_prime
    // int values
    cout << "heuristic - extend by one policy" << endl;
    // cout << "m3: " << endl;
    long begin = getCurrentTime();
    int k = policies;
    int p = scenarios;
    int m = m3.arcs;
    int k_prime = k+1;

    
    // find worst subset in optimal partition
    double min_subset_obj = GRB_INFINITY;
    int min_subset_windex;

    for (int w=0; w<k; ++w){
        if (solutions[w].objective < min_subset_obj) {
            min_subset_obj = solutions[w].objective;
            min_subset_windex = w;
        }
    }

    int p_prime = partition[min_subset_windex].size(); // the new p - number of scenarios in the subset that we will work on

    if (p_prime > 1) {
        cout << "Subset to split: { ";
        for (int q : partition[min_subset_windex]){cout << q << " ";}
        cout << "}" << endl;
        AdaptiveInstance m3_prime = AdaptiveInstance(&m3, partition[min_subset_windex]);
        // careful - m3_prime always has k=2 because we are just extending by ONE
        m3_prime.set_policies(2);
        // cout << endl << endl;
        // cout << "m3 prime: " << endl;
        // m3_prime.printInstance(G);
        
        AdaptiveSolution k_prime_solution;

        if (mip_subroutine) {
            SetPartitioningModel m3prime_model = SetPartitioningModel(500, m3_prime);
            m3prime_model.configureModel(G, m3_prime);
            m3prime_model.solve();
            k_prime_solution = m3prime_model.current_solution;
            k_prime_solution.computeAllObjectives(G, m3_prime);
        }
        else {
            k_prime_solution = enumSolve(m3_prime, G);
        }

        // auto final_solution = mergeEnumSols(k_solution, k_prime_solution, min_subset_windex);
        long time = getCurrentTime() - begin;

        mergeEnumSols(k_prime_solution, &m3_prime, min_subset_windex);
        most_recent_solution_time = time;
    }
}

