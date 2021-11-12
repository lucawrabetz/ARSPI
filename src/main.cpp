#include "../inc/M2.h"

void comp_exp_M2(vector<int>& sizes, vector<string>& graph_names, vector<int>& r_0s, vector<int>& followers_set, string& setname, string& outfile, string& experiment_logfile) {
    /*
     * Function to run a computational experiment for M2
     * Sample latex line for SFFP report: 
     * $n$ & Density & Followers & $r_0$ & MIP (s) & MIP Gap (\%) & Benders (s) & Benders Gap (\%) & Benders Cuts
     *
     */ 

    // ----- Temporary Variables -----
    vector<vector<float> > MIP_result;
    vector<float> benders_result;
    string instance_name;
    string exp_logline;

    // ----- File-Related Global Variables -----
    string graph_name;

    // ----- Latex Variables ----- 
    int n;
    int max_n = 1000000;
    float density;
    int followers; 
    float MIP_time;
    string MIP_gap;
    float benders_time;
    string benders_gap;
    int benders_cuts;
    int r_0;
    string latex_line;
    string objective_correct;
    bool policy_agree = true;
    string policy_agree_str;
    string MIP_policy = "";
    string benders_policy = "";

    // ----- Global Variables Related to Algorithm -----
    M2ProblemInstance M2;
    M2ModelLinear M2_L;
    M2Benders M2_B;
    int min = 30;
    int max = 80;

    cout << "starting experiment" << endl;
    ofstream out(outfile);
    ofstream out2(experiment_logfile);

    string csv_header = "graphname,followers,objectiveagree,policyagree,policies\n";
    out2 << csv_header;

    // for each graph 
    for (int i = 0; i<graph_names.size(); ++i){
        graph_name = graph_names[i];
        n = sizes[i];
        r_0 = r_0s[i];

        if (n>max_n){continue;}
        // read graph (set density)
        const LayerGraph G = LayerGraph(graph_name, n);
        
        // for SFFP report we just use G.m instead of density 
        // because computed it wrong in first experiment
        density = G.m;
        
        // later on: actual density is 
        // density = (float(G.m) / (n*(n-1) / 2) );
        // density = (int)(density * 100 + .5);
        // density = (float)density / 100; 

        // for each number of followers (set followers, n, r_0)
        for (int j = 0; j<followers_set.size(); ++j){
            policy_agree = true;
            MIP_policy = "";
            benders_policy = "";

            followers = followers_set[j];
            instance_name = to_string(n) + "_" + "p-" + to_string(followers);

            // set M2 instance with followers and generate costs 
            M2 = M2ProblemInstance(G, min, max, followers, 1, r_0, instance_name, setname);

            // solve with MIP, (set MIP stats: MIP, MIP Gap)
            M2_L = M2ModelLinear(&M2);
            MIP_result = M2_L.solve();
            MIP_time = M2_L.running_time;
            
            if (M2_L.optimality_gap == 0) {
                MIP_gap = "-";
            }
            else if (M2_L.optimality_gap == -1) {
                MIP_gap = "unb";
            }
            else {MIP_gap = to_string(M2_L.optimality_gap);}
            
            // solve with Benders, (set Benders stats: Benders, Benders Gap, Benders Cuts)
            M2_B = M2Benders(&M2);
            benders_result = M2_B.solve();
            benders_time = M2_B.running_time;
            
            if (M2_B.optimality_gap == 0) {
                benders_gap = "-";
            }
            else if (M2_B.optimality_gap == -1) {
                benders_gap = "unb";
            }
            else {benders_gap = to_string(M2_B.optimality_gap);}
            
            benders_cuts = M2_B.sep.cut_count;

            // if MIP and benders have same objective
            if (abs(MIP_result[0][0] - benders_result[0]) <= M2_B.sep.epsilon) {objective_correct="true";}

            else {
                // diagnostics
                objective_correct = "false";
            }

            cout << "hello" << endl;
            for (int a=0; a<G.m+1; ++a){
                // include objective
                MIP_policy +=  to_string(MIP_result[1][a]) + ",";
                benders_policy += to_string(benders_result[a]) + ",";

                if (policy_agree){
                    if (MIP_result[1][a] != benders_result[a]) {
                        policy_agree = false;
                    }
                }
            }
            
            cout << "INSTANCE RESULTS: " << endl;
            cout << "Graph, Followers: " << graph_name << ", " << followers << endl;
            cout << "MIP obj: " << MIP_result[0][0] << endl;
            cout << "Benders obj: " << benders_result[0] << endl;
            
            cout << "MIP running_time: " << MIP_time << endl;
            cout << "MIP gap: " << MIP_gap << endl;
            cout << "Benders running_time: " << benders_time << endl;
            cout << "Benders gap: " << benders_gap << endl;
            cout << "Benders cuts: " << benders_cuts << endl;
            
            if (objective_correct=="true") {cout << "OBJECTIVES MATCH" << endl;}
            else {cout << "OBJECTIVES DON'T MATCH" << endl;}
            
            cout << "--------------------------------------------------------------" << endl;
            cout << "--------------------------------------------------------------" << endl;
            cout << "--------------------------------------------------------------" << endl;
            cout << "--------------------------------------------------------------" << endl;
            cout << "--------------------------------------------------------------" << endl;
            cout << "--------------------------------------------------------------" << endl;
            cout << "--------------------------------------------------------------" << endl;
            cout << "--------------------------------------------------------------" << endl;

            // set latex string 
            // $n$ & Density & Followers & $r_0$ & MIP (s) & MIP Gap (\%) & Benders (s) & Benders Gap (\%) & Benders Cuts
            latex_line = to_string(n) + " & " + to_string(density) + " & " + to_string(followers) + " & " + to_string(r_0) + " & " + to_string(MIP_time) + " & " + (MIP_gap) + " & " + to_string(benders_time) + " & " + (benders_gap) + " & " + to_string(benders_cuts) + "\\\\" + "\n" + "\\hline" + "\n";

            if (policy_agree) {policy_agree_str="true";}
            else {policy_agree_str="false";}

            // data for each graph with interdiction policies to post process and analyze in python
            exp_logline = graph_name + "," + to_string(followers) + "," + objective_correct + "," + policy_agree_str + "," + MIP_policy + benders_policy + "\n";

            // write latex string to file, and a hline when necessary
            out << latex_line;
            out2 << exp_logline;
        }
    }       
}

int main()
{
    // TESTING MODELS
    int n=6;
    int p=4;
    int k=2;
    int r_0=2;
    const string filename = "dat/simplegraph3.txt";
    string test = "test";

    const LayerGraph G = LayerGraph(filename, n);
    M2ProblemInstance M2 = M2ProblemInstance(G, 30, 80, p, k, r_0, test, test); 
    
    auto final_solution = enumSolve(M2);
    // M2ModelLinear M_L = M2ModelLinear(&M2);
    // // M2Benders M_B = M2Benders(&M2);

    // vector<vector<float> > x_MIP = M_L.solve();

    // cout << "objective: " << x_MIP[0][0] << endl;

    // for (int w=1; w<k+1; ++w){
    //     cout << "policy " << w << ": ";

    //     for (int a=0; a<M2.m; ++a){
    //         cout << x_MIP[w][a];
    //     }
    //     cout << endl;
    // }

    // k=1;
    // M2ProblemInstance M22 = M2ProblemInstance(G, 30, 80, p, k, r_0, test, test); 
    // M2ModelLinear M_L2 = M2ModelLinear(&M2);
    // // M2Benders M_B = M2Benders(&M2);

    // vector<vector<float> > x_MIP2 = M_L2.solve();

    // cout << "objective: " << x_MIP2[0][0] << endl;

    // for (int w=1; w<k+1; ++w){
    //     cout << "policy " << w << ": ";

    //     for (int a=0; a<M22.m; ++a){
    //         cout << x_MIP2[w][a];
    //     }
    //     cout << endl;
    // }
    // vector<float> x_bend = M_B.solve();
    // float MIP_obj = x_MIP[0];
    // float bend_obj = x_bend[0];

    // x_MIP.erase(x_MIP.begin());
    // x_bend.erase(x_bend.begin());

    // float MIP_valid_obj = M2.validatePolicy(x_MIP);
    // float bend_valid_obj = M2.validatePolicy(x_bend);

    // vector<int> sp_temp;

    // for (int q=0; q<p; ++q){
    //     sp_temp=M2.Dijkstra(q);
    //     cout << "SP q=" << q << ": " << sp_temp[0] << endl;
    // }

    // cout << "mip obj: " << MIP_obj << endl;
    //cout << "benders obj: " << bend_obj << endl;
    //cout << "mip obj from dij: " << MIP_valid_obj << endl;
    //cout << "benders obj from dij: " << bend_valid_obj << endl;

    // FOR TESTING THE ENUMERATION CODE
    //vector<int> nums;
    //vector<vector<int> > result;
    //nums.push_back(0);
    //nums.push_back(1);
    //nums.push_back(2);
    //
    //result = enum_combs(3, nums);
    //for (int i=0; i<result.size(); ++i){
    //    cout << "\ncombination " << i << "\n";
    //    for (int j=0; j < result[i].size(); ++j){
    //        cout << "\n" << result[i][j] << "\n";
    //    }
    //}

    // FOR TESTING (SINGLE INSTANCE) M2 STUFF 
    // int n = 11;
    //float running_time;
    //std::vector<float> x_final_benders;
    //float x_final_linear;
    // const std::string filename = "dat/set1_08-24-21/set1_08-24-21_11_0.5.txt";
    // const LayerGraph G = LayerGraph(filename, n);
    // M2ProblemInstance M2 = M2ProblemInstance(G, 150, 160, 3, 2);
    
    // TESTING DIJSKTRA 
    // int n=6;
    // const string filename = "dat/simplegraph3.txt";
    // vector<vector<int> > empty_vector1;
    // vector<int> empty_vector2;
    // cout << "hello 0" << endl;

    // const LayerGraph G = LayerGraph(filename, n);
    // cout << "hello 1" << endl;

    // string instance_name = "simple";
    // string set_name = "simplegraph";
    // cout << "hello from main" << endl;
    // M2ProblemInstance M2 = M2ProblemInstance(G, 150, 160, 1, 0, instance_name, set_name); 
    // M2.printInstance();

    // vector<int> sp_result = M2.Dijkstra(0);

    // cout << "Path objective: " << sp_result[0] << endl;
    // for (int a=1; a<G.m+1; ++a){
    //     cout << sp_result[a] << endl;
    // }



    // COMPUTATIONAL EXPERIMENT FOR M2
    // string setname = "set1_09-17-21";
    // const string logfilename = "dat/" + setname + "/" + setname + ".log";

    // int num_instances;
    // string line;
    // ifstream myfile(logfilename);

    // int n_temp;
    // int r_0_temp;
    // vector<int> sizes; 
    // vector<string> graph_names;
    // vector<int> r_0s; 
    // vector<int> followers_set; 
    // 
    // if (myfile.is_open()){
    //     int line_counter = 0;

    //     while (getline(myfile, line)) {

    //         stringstream ss(line);
    //         string word;

    //         if (line_counter == 0) {
    //             // sizes
    //             // the size (n) will also determine the r_0 to keep it independent of graph density
    //             // for now we will set it as r_0 = floor(n * 0.5)
    //             while (ss >> word) {
    //                 n_temp = stoi(word);
    //                 sizes.push_back(n_temp);
    //                 r_0_temp = floor(0.5 * n_temp);
    //                 r_0s.push_back(r_0_temp);
    //             }
    //         }

    //         else if (line_counter == 1) {
    //             // graph names 
    //             while (ss >> word) {
    //                 graph_names.push_back(word);
    //             }
    //         }

    //         ++line_counter;
    //     }
    // }

    // // check the totals are correct
    // // cout << "Number of sizes: " << sizes.size() << endl;
    // // cout << "Number of names: " << graph_names.size() << endl;
    // // cout << "Number of r_0: " << r_0s.size() << endl;

    // // for (int i = 0; i < sizes.size(); ++i) {

    // //     cout << "Graph name: " << graph_names[i] << endl;
    // //     cout << "Graph size: " << sizes[i] << endl;

    // // }

    // string outfile = "dat/" + setname + "/outfile.txt";
    // string exp_logfile = "dat/" + setname + "/exp_logfile.csv";
    // followers_set.push_back(1);
    // followers_set.push_back(3);
    // followers_set.push_back(5);
    // followers_set.push_back(10);

    // comp_exp_M2(sizes, graph_names, r_0s, followers_set, setname, outfile, exp_logfile);

    


    // FOR TESTING MODELS

}
