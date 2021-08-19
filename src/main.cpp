#include "../inc/M3.h"

void comp_exp_M2(vector<string> graph_names, vector<int> sizes, vector<int> r_0s, vector<int> followers_set, string outfile) {
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

    
    // for each graph 
    for (int i = 0; i<graph_names.size(); ++i){
        graph_name = graph_names[i];

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
    int n = 3;
    //float running_time;
    //std::vector<float> x_final_benders;
    //float x_final_linear;
    const std::string filename = "set1_08-2-21_3_0.5.txt";
    const LayerGraph G = LayerGraph(filename, n);
    M2ProblemInstance *M2 = new M2ProblemInstance(G, 150, 160, 3, 2);

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
