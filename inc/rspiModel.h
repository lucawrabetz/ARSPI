#pragma once
#include "graph.h"
#include "/Library/gurobi902/mac64/include/gurobi_c++.h"

// rspiModel - MIP to solve a single policy robust SPI instance (MIP (7))

class rspiModel
{
private:
    GRBEnv *covenv;
    GRBModel *covmodel;

public:
    void test_function()
    {
        std::cout << "hello world\n";
    }
};
