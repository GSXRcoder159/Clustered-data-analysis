#include <vector>

using namespace std;

//an object representing a mathematical vector
struct Vect
{
    std::vector<int> v; //a vector with x, y and optionally even z coordinate
    int x;  //x-coordinate
    int y;  //y-coordinate
    int z;  //z-coordinate

    //a constructor
    Vect(std::vector<int> a);

    //get the dot product with another vector
    float dot_product(std::vector<int> b);

    //get the cross product with another vector
    float cross_product(std::vector<int> b);

    //get the angle between this and another vector
    float get_angle(std::vector<int> b);
};
