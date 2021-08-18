#include <vector>

using namespace std;

//an object representing a point in 2D space
struct Point
{
    float x;    //x-coordinate
    float y;    //y-coordinate

    //An operator for sorting Points by their y-coordinate in ascending order
    bool operator<(const Point &other) const
    {
        if (y == other.y)
        {
            return x > other.x;
        }
        else
        {
            return y < other.y;
        }
    }

    //a constructor
    Point(int x = 0, int y = 0);

    /**
     *  Checks for the mutual orientation of there pixels.
     *
     *  @param[in]  a  first Pixel.
     *  @param[in] b   second Pixel
     *  @param[in] c  third Pixel
     *  @return  (int) Value based on the orientation: 0 = colinear, 1 = clockwise,
     *    2 = counterclockwise.
     */
    int orient(Point a, Point b, Point c);

    /**
     *  Get the convex hull of a cluster. 
     * 
     *  @param[in]  points    a vector containing all Pixel from a given cluster.
     *  @param[in] size  size of the vector containing all Pixel from a given cluster.
     *  @return  (vector<Point>) a vector of the Pixel forming the convex hull
     *  of the given cluster.
     */
    void convexHull(Point points[], int m);

    /**
     *  Get a position vector a given point. 
     * 
     *  @param[in]  a    a Point object.
     *  @return  (vector<int>) a position vector.
     */
    std::vector<int> get_vector(Point a);
};