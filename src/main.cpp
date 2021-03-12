#include "../inc/M2.h"

int main()
{
    const std::string filename = "./dat/simplegraph.txt";
    const LayerGraph G = LayerGraph(filename, 6);
    const M2ProblemInstance M2 = M2ProblemInstance(G, 150, 160, 3, 2);
    M2ModelLinear M2_L = M2ModelLinear(M2);
    M2ModelBilinear M2_BL = M2ModelBilinear(M2);
    M2_L.solve();
    M2_BL.solve();
}