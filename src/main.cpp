#include "../inc/M3.h"

// use this file to test and debug
// write good scripts in other file

int main()
{
    int n=10;
    int p=5;
    int k=2;
    int M = 500;
    int interdiction_cost = 10;
    int min = 30;
    int max = 80;
    int mean = 50;
    int stddev = 10;
    const string setname = "test-03_28_22-0";
    const string name = setname + "-" + to_string(n) + "_3";
    const string directory = "dat/" + setname + "/";
    const string filename = directory + name + ".txt";

    const LayerGraph G = LayerGraph(filename, n);
    int r_0=n/4;
    AdaptiveInstance m3 = AdaptiveInstance(p, k, r_0, G, directory, name); 
    m3.initCosts(interdiction_cost, mean, stddev, 1);

    auto k_solution = enumSolve(m3, G);
    printSolution(k_solution, "k solution");

    auto kprime_solution = extendByOne(k_solution, m3, G);
    printSolution(kprime_solution, "k+1 solution (heuristic)");

    m3.set_policies(k+1);
    auto kk_solution = enumSolve(m3, G);
    printSolution(kk_solution, "k+1 solution (exact)");


    // for (int a=0; a<m3.arcs; ++a) {
    //     cout << m3.interdiction_costs[a] << " ";
    // }
    // cout << endl;
    // cout << endl;

    // for (int q=0; q<m3.scenarios; ++q) {
    //     for (int a=0; a<m3.arcs; ++a) {
    //         cout << m3.arc_costs[q][a] << " ";
    //     }
    //     cout << endl;
    //     cout << endl;
    // }
    // G.printGraph(m3.arc_costs, m3.interdiction_costs);

    // for (int k=1; k<=max_k; ++k) {
    //     cout << "\n\nk: " << k << endl;
    //     m3.set_policies(k);
    //     SetPartitioningModel sp_model = SetPartitioningModel(M, m3);
    //     sp_model.configureModel(G, m3);
    //     string sp_filename = "sp_model_k_" + to_string(k) + ".lp";
    //     sp_model.m3_model->write(sp_filename);

    //     auto final_solution = enumSolve(m3, G);
    //     vector<vector<float> > x_spmodel = sp_model.solve();
    //     cout << "MIP OPTIMAL SOLUTION: " << x_spmodel[0][0] << endl;
    // }

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
}
