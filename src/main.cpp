#include "../inc/M3.h"

void comp_exp_M2(vector<int>& sizes, vector<string>& graph_names, vector<int>& r_0s, vector<int>& followers_set, string& outfile) {
    /*
     * Function to run a computational experiment for M2
     * Sample latex line for SFFP report: 
     * $n$ & Density & Followers & $r_0$ & MIP (s) & MIP Gap (\%) & Benders (s) & Benders Gap (\%) & Benders Cuts
     *
     */ 

    // ----- Create File-Related Global Variables -----
    string graph_name;
    string outfile_name = outfile;

    // ----- Create Latex Variables ----- 
    int n;
    float density;
    int followers; 
    float MIP_time;
    float MIP_gap;
    float benders_time;
    float benders_gap;
    int benders_cuts;
    int r_0;

    // ----- Global Variables Related to Algorithm -----
    M2ProblemInstance M2;
    M2ModelLinear M2_L;
    M2Benders M2_B;
    int min;
    int max;

    cout << "in computational experiment now" << endl;

    // for each graph 
    for (int i = 0; i<graph_names.size(); ++i){
        graph_name = graph_names[i];
        cout << graph_name << endl;

        // for each number of followers (set followers, n, r_0)
        for (int j = 0; j<followers_set.size(); ++j){
            followers = followers_set[j];
            n = sizes[j];
            r_0 = r_0s[j];

            // read graph to create M2 instance, (set density)
            const LayerGraph G = LayerGraph(graph_name, n);
            M2 = M2ProblemInstance(G, min, max, followers, r_0);
            density = (float(n) / G.m );

            // solve with MIP, (set MIP stats)
            M2_L = M2ModelLinear(&M2);
            
            // solve with Benders, (set Benders stats)
            M2_B = M2Benders(&M2);
            
            // set latex string 
            
            // write latex string to file, and a hline when necessary
        }
    }       
    
}

int main()
{
    // FOR TESTING THE ENUMERATION CODE
    //vector<int> nums;
    //vector<vector<int>> result;
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

    // FOR TESTING M2 STUFF 
    int n = 11;
    //float running_time;
    //std::vector<float> x_final_benders;
    //float x_final_linear;
    const std::string filename = "dat/set1_08-24-21/set1_08-24-21_11_0.5.txt";
    const LayerGraph G = LayerGraph(filename, n);
    M2ProblemInstance M2 = M2ProblemInstance(G, 150, 160, 3, 2);

    // COMPUTATIONAL EXPERIMENT FOR M2
    const std::string logfilename = "dat/set1_08-24-21/set1_08-24-21.log";

    int num_instances;
    string line;
    ifstream myfile(logfilename);

    vector<int> sizes; 
    vector<string> graph_names;
    vector<int> r_0s; 
    vector<int> followers_set; 
    
    if (myfile.is_open()){

        const char *cline;
        int line_counter = 0;

        while (getline(myfile, line)) {

            stringstream ss(line);
            string word;

            // if (line_counter == 0) {
            //     // number of instances
            //     num_instances = stoi(line);
            // }

            if (line_counter == 0) {
                // sizes
                while (ss >> word) {
                    sizes.push_back(stoi(word));
                }
            }

            else if (line_counter == 1) {
                // graph names 
                while (ss >> word) {
                    graph_names.push_back(word);
                }
            }

            // else if (line_counter == 3) {
            //     // r_0
            //     while (ss >> word) {
            //         r_0s.push_back(stoi(word));
            //     }
            // }

            // else if (line_counter == 3) {
            //     // followers
            //     while (ss >> word) {
            //         followers_set.push_back(stoi(word));
            //     }
            // }

            ++line_counter;
        }
    }

    // check the totals are correct
    // cout << "Number of instances: " << num_instances << endl;
    cout << "Number of sizes: " << sizes.size() << endl;
    cout << "Number of names: " << graph_names.size() << endl;
    // cout << "Number of r_0: " << r_0s.size() << endl;
    // cout << "Number of followers: " << followers_set.size() << endl;

    
    for (int i = 0; i < sizes.size(); ++i) {

        cout << "Graph name: " << graph_names[i] << endl;
        cout << "Graph size: " << sizes[i] << endl;

    }

    string outfile = "exp1_08-24-21.txt";

    comp_exp_M2(sizes, graph_names, r_0s, followers_set, outfile);


    //M2ModelLinear M2_L = M2ModelLinear(M2);
    //M2_L.M2model->write("simplegraph1mip.lp");

    //M2Benders M2_Bend = M2Benders(M2);

    // M2ModelLinear (above) is the equivalent to (17)-(21) in Overleaf
    // There is no bilenear term in this model
    // If we want to experiment with the other formulation (6)-(10),
    // we could have a linear and a bilinear version for that one
    // currently the class M2ModelBilinear is unstable, it is (6)-(10)
    // without manual linearization (using gurobi quadratic constraints)

    // M2ModelBilinear M2_BL = M2ModelBilinear(M2);
    // M2_BL.M2model->write("bilinear.lp");
    // cout << "\n\n";
    // cout << "\nSolving Using LINEAR MODEL\n";
    // running_time = M2_L.solve();

    //cout << "\n\n";
    //cout << "\nSolving Using LINEAR MODEL\n";

    //x_final_linear = M2_L.solve();

    // printout for the linear model
    // cout << "\nObjective: " << x_final[0] << "\n";
    // for (int i = 1; i < n + 1; ++i)
    // {
    //     cout << "\nx_" << i - 1 << ": " << x_final[i] << "\n";
    // }

    //x_final_benders = M2_Bend.solve();
    //// printout for the benders model
    //cout << "\nObjective: " << x_final_benders[0] << "\n";
    //for (int a = 1; a < G.m + 1; ++a)
    //{
    //    cout << "\nx_" << a - 1 << ": " << x_final_benders[a] << "\n";
    //}

    //// cout << "\n\n\n\n\n";
    //// cout << "\nSolving Using BiLINEAR MODEL\n";
    //// running_time = M2_BL.solve();
    //delete M2_L.M2env;
    //delete M2_L.M2model;
    //// delete M2_BL.M2env;
    //// delete M2_BL.M2model;
    //delete M2_Bend.M2Bendersenv;
    //delete M2_Bend.M2Bendersmodel;

    //delete M2;
}
