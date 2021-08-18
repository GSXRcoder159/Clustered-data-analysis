using namespace std;

// A structure storing information about a pair of clusters
struct Pair
{
    int l1; //the index of the cluster from layer 1 in the vector of clusters form layer 1
    int l2; //the index of the cluster from layer 2 in the vector of clusters form layer 2
    int dsize;  //the size difference between the two clusters

    //a constructor
    Pair(int l1 = 0, int l2 = 0, int dsize = 0);
};
