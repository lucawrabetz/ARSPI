#include "../inc/M2.h"

using std::cout;

int main()
{
    float running_time;
    const std::string filename = "./dat/simplegraph1.txt";
    const LayerGraph G = LayerGraph(filename, 2);
    M2ProblemInstance *M2 = new M2ProblemInstance(G, 150, 160, 3, 2);
    M2ModelLinear M2_L = M2ModelLinear(M2);
    M2_L.M2model->write("simplegraph1.lp");

    // M2ModelLinear (above) is the equivalent to (17)-(21) in Overleaf
    // There is no bilenear term in this model
    // If we want to experiment with the other formulation (6)-(10),
    // we could have a linear and a bilinear version for that one
    // currently the class M2ModelBilinear is unstable, it is (6)-(10)
    // without manual linearization (using gurobi quadratic constraints)

    // M2ModelBilinear M2_BL = M2ModelBilinear(M2);
    // M2_BL.M2model->write("bilinear.lp");

    cout << "\n\n";
    cout << "\nSolving Using LINEAR MODEL\n";
    running_time = M2_L.solve();

    // cout << "\n\n\n\n\n";
    // cout << "\nSolving Using BiLINEAR MODEL\n";
    // running_time = M2_BL.solve();
    delete M2_L.M2env;
    delete M2_L.M2model;
    // delete M2_BL.M2env;
    // delete M2_BL.M2model;
    delete M2;
}