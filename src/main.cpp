#include "../inc/M2.h"

int main()
{
    const std::string filename = "./dat/simplegraph.txt";
    const LayerGraph G = LayerGraph(filename, 6);
    M2ProblemBilinear M2 = M2ProblemBilinear(G, 150, 160, 3, 2);
    M2.solve();
}