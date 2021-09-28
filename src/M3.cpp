#include "../inc/M3.h"

vector<vector<int>> enum_combs(int a, vector<int>& nums){
    // INPUTS:
    // a - size of each combination
    // nums - nums to pick from
    int n = nums.size();
    vector<vector<int>> result;
    vector<int> temp;
    
    // base case - a == nums.size() (nums is the only combination)
    if (a == n) {
        temp = nums;
        result.push_back(temp);
        return result;
    }
}

void enum_partitions(int k, vector<int>& p, vector<vector<int>>& so_far, vector<vector<vector<int>>>* result){
    // INPUTS:  
    // int k - the number of subsets in every partition
    // vector<int> p - the full universe of indexes to choose from (at the top of the recursion tree, 0,1, ... , p)
    // vector<vector<int>> so_far - a partition we are fixed to in this branch of the tree - empty at top level call 
    // vector<vector<vector<int>>>* result - the final vector of partitions, we are working in place 

    // base case - k == 1

    // if k > 1, recursive function
}
