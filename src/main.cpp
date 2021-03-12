#include "../inc/M2.h"

int main()
{
    const std::string filename = "./dat/simplegraph.txt";
    const LayerGraph G = LayerGraph(filename, 6);
    M2ProblemLinear M2_L = M2ProblemLinear(G, 150, 160, 3, 2);
    M2ProblemBilinear M2_BL = M2ProblemBilinear(G, 150, 160, 3, 2);
    M2_L.solve();
    M2_BL.solve();
}