#include <string>
#include <iostream>
#include <vector>

using namespace std;

//an object representing a line with its atributes
struct line
{
    int x;  //x-coordinate
    int y;  //y-coordinate
    float angle = 0; //angle with x-axis

    //a constructor
    line(int x, int y, float angle);
    /**
     *  Calculate the y-coordinate.
     * 
     *  @param[in]  x   the x-coordinate.
     *  @return  (int) y-coordinate.
     */
    int get_y(int x);
    
    /**
     *  Calculate the x-coordinate.
     * 
     *  @param[in]  y   the y-coordinate.
     *  @return  (int) x-coordinate.
     */
    int get_x(int y);

    /**
     *  Calculate the x-coordinate of an intercept between two lines.
     * 
     *  @param[in]  l1   first line object.
     *  @param[in]  l2   second line object.
     *  @return  (vector<int>) x-coordinate of the intercept.
     */
    vector<int> get_x_intercept(line l1, line l2);
    
    /**
     *  Calculate the y-coordinate of an intercept between two lines.
     * 
     *  @param[in]  l1   first line object.
     *  @param[in]  l2   second line object.
     *  @return  (vector<int>) y-coordinate of the intercept.
     */
    vector<int> get_y_intercept(line l1, line l2);
};
