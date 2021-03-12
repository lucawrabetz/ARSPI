#include "../inc/rspi.h"

rspiModel::rspiModel()
{
    /* 
        Default constructor
    */
    n, m, k, t = 0;
}

rspiModel::rspiModel(const std::string &filename)
{
    /* 
        Read the graph from a file - construct RSPI class
        Construct base gurobi model
    */
}

void rspiModel::update()
{
    /*
        Add main constraints to the gurobi model (depend on subset of cost vectors)
    */
}

void rspiModel::solve(std::vector<int> &data)
{
    /*
        Solve the robust model for each subset of costs in the k-partition
        - If the minimum out of the k solutions is better than our current v, update v, and update x_bar
        - If not, do nothing
        - In the bigger scheme of this package, this function runs at the bottom of the partitioning recursion
    */
    a = data.data();
    next = 0;
    z = -1;

    for (int j = 0; j < k; j++)
    // there are k sets of scenarios in this partition
    {
        // update current_scenarios to reflect new set, based on data, which defines the partition split
        current_scenarios.clear();
        current_scenarios.push_back(next);
        next++;
        while (*a == 0)
        {
            current_scenarios.push_back(next);
            ++a;
            next++;
        }

        update(); // we need update the gurobi model with this subset of scenarios
        model->optimize();

        if (z < 0)
        {
            z = model->get(GRB_DoubleAttr_ObjVal);
        }
        else if (model->get(GRB_DoubleAttr_ObjVal) < z)
        {
            z = model->get(GRB_DoubleAttr_ObjVal);
        }
    }

    if (z > v)
    {
        v = z;
    }
}