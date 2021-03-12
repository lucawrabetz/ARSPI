#include "../inc/rspi.h"

void subsetRecur(std::vector<int> &arr, int l, int k, int index, std::vector<int> &data, int i);

void subsetLoop(int l, int k)
{
    // A temporary vector to store current subset
    std::vector<int> data;
    // Original values (0->q-1)
    std::vector<int> arr;

    for (int i = 0; i < l; i++)
    {
        if (i < k)
        {
            data.push_back(0);
        }

        arr.push_back(i);
    }

    // Pass each subset of size k to rspi
    subsetRecur(arr, l, k, 0, data, 0);
}

/* arr[] ---> Input Array 
n	 ---> Size of input array 
r	 ---> Size of a combination to be printed 
index ---> Current index in data[] 
data[] ---> Temporary array to store current combination 
i	 ---> index of current element in arr[]	 */
void subsetRecur(std::vector<int> &arr, int l, int k, int index, std::vector<int> &data, int i)
{
    // Current subset is ready, run rspi.solve with it
    if (index == k)
    {
        // return candidate subset instead of printing
        // RUN THE INSIDE OF THE ALGORITHM LOOP HERE SOMEHOW::
        // data: the k indexes that make up the subset
        // rspi.solve(data)

        return;
    }

    // When no more elements are there to put in data[]
    if (i >= l)
    {
        return;
    }

    // recursions
    // index + 1, i + 1
    data[index] = arr[i];
    subsetRecur(arr, l, k, index + 1, data, i + 1);

    // index, i + 1
    subsetRecur(arr, l, k, index, data, i + 1);
}

int main()
{
    int k = 3;
    int l = 5;
    subsetLoop(l, k);
    return 0;
    // graph = directedGraph()

    // abstraction of ARSPI object

    // v = -infinity
    // X = k vectors, final policies

    // partition U
}