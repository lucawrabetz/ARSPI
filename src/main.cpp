#include "../inc/M3.h"

// use this file to test and debug
// write good scripts in other file

int main()
{
    int n=1202;
    int p=5;
    int k=1;
    int M = 500;
    int interdiction_cost = 5;
    int min = 30;
    int max = 80;
    int mean = 50;
    int stddev = 10;
    int a = 50;
    int b = 200;
    float fraction = 0.1;

    const string setname = "tests-06_20_22-0";
    const string name = setname + "-" + to_string(n) + "_2_3";
    const string directory = "dat/" + setname + "/";
    const string filename = directory + name + ".txt";

    vector<int> fake_vec = vector<int>();
    vector<vector<int>> fake_vec2 = vector<vector<int>>();

    const Graph G = Graph(filename, n);
    int r_0=1500;
    cout << "m: " << G.m << endl;
    AdaptiveInstance m3 = AdaptiveInstance(p, k, r_0, G, directory, name); 
    m3.initCosts(fraction, a, b, 3, G, true);
    
    cout << G.m << endl;

    // k = 1, enumeration
    // AdaptiveSolution enum_solution = enumSolve(m3, G);
    // double exact_objective_1 = enum_solution.worst_case_objective;
    // double enum_time_1 = enum_solution.most_recent_solution_time * 0.001;
    // enum_solution.logSolution(G, m3, "k=1 solution exact");

    // k = 1, sp
    SetPartitioningModel sp_model(M, m3);
    sp_model.configureModel(G, m3);
    sp_model.solve();
    AdaptiveSolution mip_solution = sp_model.current_solution;
    mip_solution.computeAllObjectives(G, m3);
    double mip_time_1 = mip_solution.most_recent_solution_time * 0.001;
    double exact_objective_1 = mip_solution.worst_case_objective;
    // mip_solution.logSolution(G, m3, "k=1 solution exact, sp model");

    string heuristic_objective_1 = "/";
    string heuristic_time_1 = "/";
    string enum_time_1 = "/";
    string enum_time_2 = "/";
    string enum_time_3 = "/";
    string enum_time_4 = "/";
    string enum_time_5 = "/";

    // k = 2
    mip_solution.extendByOne(m3, G);
    mip_solution.computeAllObjectives(G, m3);
    double heuristic_objective_2 = mip_solution.worst_case_objective;
    double heuristic_time_2 = mip_solution.most_recent_solution_time * 0.001;
    // enum_solution.logSolution(G, m3, "");

    // k = 3
    mip_solution.extendByOne(m3, G);
    mip_solution.computeAllObjectives(G, m3);
    double heuristic_objective_3 = mip_solution.worst_case_objective;
    double heuristic_time_3 = mip_solution.most_recent_solution_time * 0.001;
    // enum_solution.logSolution(G, m3, "");
    
    // k = 4
    mip_solution.extendByOne(m3, G);
    mip_solution.computeAllObjectives(G, m3);
    double heuristic_objective_4 = mip_solution.worst_case_objective;
    double heuristic_time_4 = mip_solution.most_recent_solution_time * 0.001;
    // enum_solution.logSolution(G, m3, "");
    

    // k = 5
    mip_solution.extendByOne(m3, G);
    mip_solution.computeAllObjectives(G, m3);
    double heuristic_objective_5 = mip_solution.worst_case_objective;
    double heuristic_time_5 = mip_solution.most_recent_solution_time * 0.001;
    // enum_solution.logSolution(G, m3, "");
    
    m3.set_policies(2);
    // AdaptiveSolution enum_solution2 = enumSolve(m3, G);
    // double enum_time_2 = enum_solution2.most_recent_solution_time * 0.001;
    SetPartitioningModel sp_model2(M, m3);
    sp_model2.configureModel(G, m3);
    sp_model2.solve();
    AdaptiveSolution mip_solution2 = sp_model2.current_solution;
    mip_solution2.computeAllObjectives(G, m3);
    double exact_objective_2 = mip_solution2.worst_case_objective;
    double mip_time_2 = mip_solution2.most_recent_solution_time * 0.001;

    m3.set_policies(3);
    // AdaptiveSolution enum_solution2 = enumSolve(m3, G);
    // double enum_time_2 = enum_solution2.most_recent_solution_time * 0.001;
    SetPartitioningModel sp_model3(M, m3);
    sp_model3.configureModel(G, m3);
    sp_model3.solve();
    AdaptiveSolution mip_solution3 = sp_model3.current_solution;
    mip_solution3.computeAllObjectives(G, m3);
    double exact_objective_3 = mip_solution3.worst_case_objective;
    double mip_time_3 = mip_solution3.most_recent_solution_time * 0.001;

    m3.set_policies(4);
    // AdaptiveSolution enum_solution2 = enumSolve(m3, G);
    // double enum_time_2 = enum_solution2.most_recent_solution_time * 0.001;
    SetPartitioningModel sp_model4(M, m3);
    sp_model4.configureModel(G, m3);
    sp_model4.solve();
    AdaptiveSolution mip_solution4 = sp_model4.current_solution;
    mip_solution4.computeAllObjectives(G, m3);
    double exact_objective_4 = mip_solution4.worst_case_objective;
    double mip_time_4 = mip_solution4.most_recent_solution_time * 0.001;

    m3.set_policies(5);
    // AdaptiveSolution enum_solution2 = enumSolve(m3, G);
    // double enum_time_2 = enum_solution2.most_recent_solution_time * 0.001;
    SetPartitioningModel sp_model5(M, m3);
    sp_model5.configureModel(G, m3);
    sp_model5.solve();
    AdaptiveSolution mip_solution5 = sp_model5.current_solution;
    mip_solution5.computeAllObjectives(G, m3);
    double exact_objective_5 = mip_solution5.worst_case_objective;
    double mip_time_5 = mip_solution5.most_recent_solution_time * 0.001;

    // k = 1 print line
    cout << exact_objective_1 << " & " << heuristic_objective_1 << " & " << mip_time_1 << " & " << enum_time_1 << " & " << heuristic_time_1 << endl;
    // k = 2 print line
    cout << exact_objective_2 << " & " << heuristic_objective_2 << " & " << mip_time_2 << " & " << enum_time_2 << " & " << heuristic_time_2 << endl;
    // k = 3 print line
    cout << exact_objective_3 << " & " << heuristic_objective_3 << " & " << mip_time_3 << " & " << enum_time_3 << " & " << heuristic_time_3 << endl;

    // k = 4 print line
    cout << exact_objective_4 << " & " << heuristic_objective_4 << " & " << mip_time_4 << " & " << enum_time_4 << " & " << heuristic_time_4 << endl;
    // k = 5 print line
    cout << exact_objective_5 << " & " << heuristic_objective_5 << " & " << mip_time_5 << " & " << enum_time_5 << " & " << heuristic_time_5 << endl;




    // k = 2 exact
    // m3.set_policies(2);
    // SetPartitioningModel sp_model2(M, m3);
    // sp_model2.configureModel(G, m3);
    // sp_model2.solve();
    // AdaptiveSolution mip_solution2 = sp_model2.current_solution;
    // mip_solution2.computeAllObjectives(G, m3);
    // mip_solution2.logSolution(G, m3, "k = 2 solution exact sp model");

    // k = 3 exact
    // m3.set_policies(3);
    // SetPartitioningModel sp_model3(M, m3);
    // sp_model3.configureModel(G, m3);
    // sp_model3.solve();
    // AdaptiveSolution mip_solution3 = sp_model3.current_solution;
    // mip_solution3.computeAllObjectives(G, m3);
    // mip_solution3.logSolution(G, m3, "k = 3 solution exact sp model");

    // // k = 4
    // enum_solution.extendByOne(m3, G);
    // enum_solution.logSolution("");

    // // k = 5
    // enum_solution.extendByOne(m3, G);
    // enum_solution.logSolution("");

    // k = 2
    // mip_solution.extendByOne(m3, G);
    // mip_solution.logSolution("k=2 after heuristic on sp solution");

    // m3.set_policies(1);
    // k = 1
    // AdaptiveSolution solution2 = enumSolve(m3, G);
    // solution2.logSolution("k=1 solution exact");

    // k = 2
    // solution2.extendByOne(m3, G, false);
    // solution2.logSolution("k=2 after heuristic (using enum as subroutine)");
    // m3.set_policies(2);
    // auto k2exact_solution = enumSolve(m3, G);
    // printSolution(k2exact_solution, "k=2 solution (exact)");

    // // k = 3
    // auto k3_solution = extendByOne(k2exact_solution, m3, G);
    // printSolution(k3_solution, "k=3 solution (heuristic, using exact k=2 as input)");

    // m3.set_policies(3);
    // auto k3exact_solution = enumSolve(m3, G);
    // printSolution(k3exact_solution, "k=3 solution (exact)");

    // k = 4
    // auto k4_solution = extendByOne(k3exact_solution, m3, G);
    // printSolution(k4_solution, "k=4 solution (heuristic, using exact k=3 as input)");

    // m3.set_policies(4);
    // auto k4exact_solution = enumSolve(m3, G);
    // printSolution(k4exact_solution, "k=4 solution (exact)");

    // // k = 5
    // auto k5_solution = extendByOne(k4exact_solution, m3, G);
    // printSolution(k5_solution, "k=5 solution (heuristic, using exact k=4 as input)");

    // m3.set_policies(5);
    // auto k5exact_solution = enumSolve(m3, G);
    // printSolution(k5exact_solution, "k=5 solution (exact)");










    // int n=20;
    // int p=3;
    // int k=1;
    // int M = 500;
    // int interdiction_cost = 10;
    // int min = 30;
    // int max = 80;
    // int mean = 50;
    // int stddev = 10;
    // int a = 50;
    // int b = 100;
    // const string setname = "test-04_04_22-0";
    // const string name = setname + "-" + to_string(n) + "_3";
    // const string directory = "dat/" + setname + "/";
    // const string filename = directory + name + ".txt";

    // const Graph G = Graph(filename, n);
    // int r_0=n*0.3;
    // AdaptiveInstance m3 = AdaptiveInstance(p, k, r_0, G, directory, name); 
    // m3.initCosts();

    // auto k_solution = enumSolve(m3, G);
    // printSolution(k_solution, "k solution");

    // auto kprime_solution = extendByOne(k_solution, m3, G);
    // printSolution(kprime_solution, "k+1 solution (heuristic)");

    // m3.set_policies(k+1);
    // auto kk_solution = enumSolve(m3, G);
    // printSolution(kk_solution, "k+1 solution (exact)");
}
