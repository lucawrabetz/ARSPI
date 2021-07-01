#pragma once
#include "../inc/M2.h"

vector<vector<int>> enum_combs(int a, vector<int>& nums);
void enum_partitions(int k, vector<int>& p, vector<vector<int>>& so_far, vector<vector<vector<int>>>* result);

class M3ProblemInstance
{
    // Full Instance of an M3 problem (to be instantiated, just needs an M2 problem, 
    // and the integer k
public:
    M2ProblemInstance M2;
    int k;

    M3ProblemInstance();
    M3ProblemInstance(const M2ProblemInstance &the_M2, int the_k);

};


class M3MIP
{
};
