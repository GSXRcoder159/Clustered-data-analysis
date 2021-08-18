/**
    Utef: Cluster analysis
    @file clusters.cpp
    @author Adam Kuča
    @version 1.1 4/10/20 
*/

#include "cluster.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <bits/stdc++.h>
#include <algorithm>

#define llu long long int
#define PI 3.14159265

using namespace std;

// Boint is an additional structure created solely to help calculate the convex hull of a data set
struct Boint
{

    llu x, y;

    // Operator for sorting points/boints ascendingly based on the x-coordinate
    bool operator<(Boint p)
    {
        return x < p.x || (x == p.x && y < p.y);
    }
};

// Returns cross product of three points
llu cross_product(Boint O, Boint A, Boint B)
{
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// Calculate the convex hull of a vector of points
vector<Boint> convex_hull(vector<Boint> A)
{
    int n = A.size(), k = 0;

    if (n <= 3)
        return A;

    vector<Boint> ans(2 * n);

    // Sort points lexicographically
    sort(A.begin(), A.end());

    // Build lower hull
    for (int i = 0; i < n; ++i)
    {

        // If the point at K-1 position is not a part
        // of hull as vector from ans[k-2] to ans[k-1]
        // and ans[k-2] to A[i] has a clockwise turn
        while (k >= 2 && cross_product(ans[k - 2],
                                       ans[k - 1], A[i]) <= 0)
            k--;
        ans[k++] = A[i];
    }

    // Build upper hull
    for (size_t i = n - 1, t = k + 1; i > 0; --i)
    {

        // If the point at K-1 position is not a part
        // of hull as vector from ans[k-2] to ans[k-1]
        // and ans[k-2] to A[i] has a clockwise turn
        while (k >= t && cross_product(ans[k - 2],
                                       ans[k - 1], A[i - 1]) <= 0)
            k--;
        ans[k++] = A[i - 1];
    }

    // Resize the array to desired size
    ans.resize(k - 1);

    return ans;
}

// Retruns the distance between two points
double dist(Boint a, Boint b)
{
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

// Calculate the perimiter of a vector of points
double perimeter(vector<Boint> ans)
{
    double perimeter = 0.0;

    // Find the distance between adjacent points
    for (int i = 0; i < ans.size() - 1; i++)
    {
        perimeter += dist(ans[i], ans[i + 1]);
    }

    // Add the distance between first and last point
    perimeter += dist(ans[0], ans[ans.size() - 1]);

    return perimeter;
}

// Cluster constructor
Cluster::Cluster()
{
}

/**
   *  Get the total energy of a given cluster.
   * 
   *  Adds together values of energies of pixels in a given cluster, then returns this value.
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> l1_clusters].
   *  @return  (float) total energy of a cluster.
  */
float Cluster::get_total_energy(int n)
{
    float totalEnergy = 0; //sum of the energies of pixels

    for (int i = 0; i < l1_clusters[n].size(); i++)
    {
        totalEnergy += l1_clusters[n][i].e;
    }

    return totalEnergy;
}

/**
   *  Get the minimum time in the given cluster. - 
   * 
   *  [Sorts times of all pixels in an ascending order, then return the first(smallest) value.]
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (double) minimum time in a cluster.
   * 
   *  
  */
double Cluster::get_minimum_time(int n)
{
    vector<double> t; //a vector with times of all pixels within a cluster

    for (int i = 0; i < l1_clusters[n].size(); i++)
    {
        t.push_back(l1_clusters[n][i].t);
    }

    std::sort(t.begin(), t.end());

    return t[0];
}

/**
   *  Get the number of elements in the given cluster.
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (int) total number of elements in a cluster.
  */
int Cluster::get_size(int n)
{
    return l1_clusters[n].size();
}

/**
   * Get the pixel with the smallest value of time in a given cluster. - 
   * 
   * [Gets the index of the pixel with the smallest time and returns the Pixel object stored under
   * that index.]
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  a Pixel object representing the pixel with the smallest time.
  */
Pixel Cluster::get_pixel_with_minimum_time(int n)
{
    int i = get_minimum_time_index(n, "l1_clusters"); //index of the pixel with the smallest time

    return l1_clusters[n][i];
}

/**
   *  Get the pixel with the greatest value of time in a given cluster. - 
   * 
   *  [Gets the index of the pixel with the greatest time and returns the Pixel object stored under
   *  that index.]
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  a Pixel object representing the pixel with the greatest time.
  */
Pixel Cluster::get_pixel_with_maximum_time(int n)
{
    int i = get_maximum_time_index(n); //index of the pixel with the greatest time

    return l1_clusters[n][i];
}

/**
   *  Get the time difference between the smallest and the greatest time in a give cluster. - 
   * 
   *  [Sorts a vector of time values in an ascending order and returns the difference between the
   *  greatest and the smallest value.]
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (double) value of the time difference.
  */
double Cluster::get_delta_time(int n)
{
    std::vector<double> t; //a vector with times of all pixels within a cluster

    for (int i = 0; i < l1_clusters[n].size(); i++)
    {
        t.push_back(l1_clusters[n][i].t);
    }

    sort(t.begin(), t.end());

    return (t[t.size() - 1] - t[0]);
}

/**
   *  Get the index of the pixel with the smallest time within a give cluster. - 
   * 
   *  [Sorts a vector of time values in an ascending order and searches for the pixel with the smallest 
   *  time, then it returns this index]
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @param[in] vect   a string containing either "l1_clusters" - for a cluster stored in l1_clusters, or
   *  "l2_clusters" for a cluster sotred in l2_clusters.
   *  @return  (int) index of the pixel with the smallest time.
  */
int Cluster::get_minimum_time_index(int n, string vect)
{
    vector<double> t; //a vector with times of all pixels within a cluster
    int index = 0;    //the index of the pixel with the smallest time
    if (vect == "l1_clusters")
    {
        for (int i = 0; i < l1_clusters[n].size(); i++)
        {
            t.push_back(l1_clusters[n][i].t);
        }

        std::sort(t.begin(), t.end());

        for (int i = 0; i < l1_clusters[n].size(); i++)
        {
            if (l1_clusters[n][i].t == t[0])
            {
                index = i;
                break;
            }
        }
    }
    else if (vect == "l2_clusters")
    {
        for (int i = 0; i < l2_clusters[n].size(); i++)
        {
            t.push_back(l2_clusters[n][i].t);
        }

        std::sort(t.begin(), t.end());

        for (int i = 0; i < l2_clusters[n].size(); i++)
        {
            if (l2_clusters[n][i].t == t[0])
            {
                index = i;
                break;
            }
        }
    }
    return index;
}

/**
   *  Get the index of the pixel with the largest time within a give cluster. - 
   * 
   *  [Sorts a vector of time values in an ascending order and searches for the pixel with the largest 
   *  time, then it returns this index]
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (int) index of the pixel with the largest time.
  */
int Cluster::get_maximum_time_index(int n)
{
    vector<double> t; //a vector with times of all pixels within a cluster
    int index = 0;    //the index of the pixel with the greatest time

    for (int i = 0; i < l1_clusters[n].size(); i++)
    {
        t.push_back(l1_clusters[n][i].t);
    }

    std::sort(t.begin(), t.end());

    for (int i = 0; i < l1_clusters[n].size(); i++)
    {
        if (l1_clusters[n][i].t == t[t.size() - 1])
        {
            index = i;
            break;
        }
    }

    return index;
}

/**
   *  Get the width of a given cluster. - 
   * 
   *  [Sorts a vector of x-coordinates in an ascending order, then finds the difference between the
   *   greatest x-coordinate and the smallest x-coordinate]
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (int) width of the cluster.
  */
int Cluster::get_x_dimension(int n)
{
    std::vector<int> x; //a vector with the x-coordinates of all vectors in a cluster

    for (int i = 0; i < l1_clusters[n].size(); i++)
    {
        x.push_back(l1_clusters[n][i].x);
    }

    sort(x.begin(), x.end());

    return int(x[x.size() - 1] - x[0]);
}

/**
   *  Get the height of a given cluster. - 
   * 
   *  [Sorts a vector of y-coordinates in an ascending order, then finds the difference between the
   *   greatest x-coordinate and the smallest y-coordinate]
   * 
   *  @param[in]  n   index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (int) height of the cluster.
  */
int Cluster::get_y_dimension(int n)
{
    std::vector<int> y; //a vector with the y-coordinates of all vectors in a cluster

    for (int i = 0; i < l1_clusters[n].size(); i++)
    {
        y.push_back(l1_clusters[n][i].y);
    }

    sort(y.begin(), y.end());

    return int(y[y.size() - 1] - y[0]);
}

/**
   *  Output data into a given file. - 
   * 
   *  [Based on a given function name, this function outputs data returned by that function into 
   *  a file located under a given path]
   * 
   *  @param[in]  path   a string containing an absolute path to the output file.
   *  @param[in]  function   a string containing the name of the function to be outputted 
   *   ["get_size", "get_delta_time",
   *  "get_total_energy", "get_minimum_time", "get_x_pixel_with_minimum_time", 
   *  "get_y_pixel_with_minimum_time", "get_x_dimension", "get_y_dimension", "get_dimension", 
   *  "get_pairs", "get_convex_hull", "get_delta_time_between_layers"].
   *  @return  Nothing.
  */
void Cluster::output(std::string path, string function)
{
    ofstream outdata; //output stream

    outdata.open(path, std::ofstream::out | std::ofstream::trunc); //open file
    outdata.close();                                               //close the file, therefore also delete all contents

    outdata.open(path); //open the file again

    outdata.clear(); //clear the output stream

    if (function == "get_size" || function == "get_delta_time")
    {
        for (int i = 0; i < l1_clusters.size(); i++)
        {
            outdata << to_string(get_size(i)) + " " + to_string(get_delta_time(i)) << endl;
        }
    }
    else
    {
        if (function == "get_total_energy")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                outdata << to_string(get_total_energy(i)) << endl;
            }
        }
        else if (function == "get_minimum_time")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                outdata << to_string(get_minimum_time(i)) << endl;
            }
        }
        else if (function == "get_x_pixel_with_minimum_time")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                outdata << to_string(get_pixel_with_minimum_time(i).x) << endl;
            }
        }
        else if (function == "get_y_pixel_with_minimum_time")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                outdata << to_string(get_pixel_with_minimum_time(i).y) << endl;
            }
        }
        else if (function == "get_x_pixel_with_maximum_time")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                outdata << to_string(get_pixel_with_maximum_time(i).x) << endl;
            }
        }
        else if (function == "get_y_pixel_with_maximum_time")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                outdata << to_string(get_pixel_with_maximum_time(i).y) << endl;
            }
        }
        else if (function == "get_x_dimension")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                outdata << to_string(get_x_dimension(i)) << endl;
            }
        }
        else if (function == "get_y_dimension")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                outdata << to_string(get_y_dimension(i)) << endl;
            }
        }
        else if (function == "get_dimensions")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                for (int j = 0; j < l1_clusters[i].size(); j++)
                {
                    outdata << to_string(l1_clusters[i][j].x) + " " + to_string(l1_clusters[i][j].y) + "\t";
                }
                outdata << "" << endl;
            }
        }
        else if (function == "get_pairs")
        {
            for (int i = 0; i < l1_clusters.size(); i++)
            {
                get_pairs(i);
                if (!pairs[i].empty())
                {
                    for (int j = 0; j < pairs[i].size(); j++)
                    {
                        for (int k = 0; k < l1_clusters[i].size(); k++)
                        {
                            outdata << to_string(l1_clusters[i][k].x) + " " + to_string(l1_clusters[i][k].y) + " " + to_string(l1_clusters[i][k].t) << endl;
                        }
                        outdata << '$' << endl;
                        for (int k = 0; k < l2_clusters[pairs[i][j]].size(); k++)
                        {
                            outdata << to_string(l2_clusters[pairs[i][j]][k].x) + " " + to_string(l2_clusters[pairs[i][j]][k].y) + " " + to_string(l2_clusters[pairs[i][j]][k].t) << endl;
                        }
                        outdata << '#' << endl;
                    }
                }
            }
        }
        else if (function == "get_convex_hull")
        {
            vector<Point> border; //a vector containing the points making up the border of the convex hull

            for (int i = 0; i < l1_clusters.size(); i++)
            {
                border = get_convex_hull(l1_clusters[i], l1_clusters[i].size());

                if (!border.empty())
                {
                    for (int j = 0; j < border.size(); j++)
                    {
                        outdata << to_string(border[j].x) + " " + to_string(border[j].y) << endl;
                    }
                }
                else
                {
                    for (int j = 0; j < l1_clusters[i].size(); j++)
                    {
                        outdata << to_string(l1_clusters[i][j].x) + " " + to_string(l1_clusters[i][j].y) << endl;
                    }
                }

                outdata << '#' << endl;
            }
        }
        else if (function == "get_delta_time_between_layers")
        {
            vector<Pair> dt = get_delta_size_between_layers(); //a vector containing the size differences between
                                                               //layers

            for (int i = 0; i < dt.size(); i++)
            {
                outdata << to_string(dt[i].l1) + " " + to_string(dt[i].l2) + " " + to_string(dt[i].dsize) << endl;
            }
        }
    }

    outdata.close(); //close the output stream
}

/**
   *  Output data into a given file, for export to MS excel. - 
   * 
   *  [This function exports data from the original clustered file into a new text file,
   *  which is more easily legible by MS excel]
   * 
   *  @param[in]  path   a string containing an absolute path to the output file.
   *  @return  Nothing.
  */
void Cluster::exp(std::string path)
{
    ofstream outdata;

    outdata.open(path, std::ofstream::out | std::ofstream::trunc);
    outdata.close();

    outdata.open(path);

    outdata.clear();

    outdata << "Energy,";
    for (int i = 0; i < l1_clusters.size(); i++)
    {
        for (int j = 0; j < l1_clusters[i].size(); j++)
        {
            outdata << ' ' + to_string(l1_clusters[i][j].e);
        }
        outdata << '\t';
    }
    outdata << '\n';

    outdata << "Time,";
    for (int i = 0; i < l1_clusters.size(); i++)
    {
        for (int j = 0; j < l1_clusters[i].size(); j++)
        {
            outdata << ' ' + to_string(l1_clusters[i][j].t);
        }
        outdata << '\t';
    }
    outdata << '\n';

    outdata << "X,";
    for (int i = 0; i < l1_clusters.size(); i++)
    {
        for (int j = 0; j < l1_clusters[i].size(); j++)
        {
            outdata << ' ' + to_string(l1_clusters[i][j].x);
        }
        outdata << '\t';
    }
    outdata << '\n';

    outdata << "Y,";
    for (int i = 0; i < l1_clusters.size(); i++)
    {
        for (int j = 0; j < l1_clusters[i].size(); j++)
        {
            outdata << ' ' + to_string(l1_clusters[i][j].y);
        }
        outdata << '\t';
    }
    outdata << '\n';

    outdata.close();
}

/**
   *  Get the convex hull of a cluster. 
   * 
   *  @param[in]  points    a vector containing all Pixel from a given cluster.
   *  @param[in] size  size of the vector containing all Pixel from a given cluster.
   *  @return  (vector<Point>) a vector of the Pixel forming the convex hull
   *  of the given cluster.
  */
vector<Point> Cluster::get_convex_hull(vector<Pixel> points, int size)
{
    vector<Point> borderPoints;
    Point point;
    if (size < 3) //at least three points required
        return borderPoints;
    int n[size];
    for (int i = 0; i < size; i++)
        n[i] = -1;
    int l = 0; //initialize result.
    for (int i = 1; i < size; i++)
        if (points[i].x < points[l].x)
            l = i; //find left most point
    int p = l, q;
    int limit = 0;
    do
    {
        q = (p + 1) % size;
        for (int i = 0; i < size; i++)
            if (orient_points(points[p], points[i], points[q]) == 2)
                q = i;
        n[p] = q;
        p = q;

        if (limit > 10 * points.size())
        {
            return borderPoints;
        }
        limit++;
    } while (p != l);
    for (int i = 0; i < size; i++)
    {
        if (n[i] != -1)
        {
            point = {points[i].x, points[i].y};
            borderPoints.push_back(point);
        }
    }

    return borderPoints;
}

/**
   *  Checks for the mutual orientation of there pixels.
   *
   *  @param[in]  a  first Pixel.
   *  @param[in] b   second Pixel
   *  @param[in] c  third Pixel
   *  @return  (int) Value based on the orientation: 0 = colinear, 1 = clockwise,
   *    2 = counterclockwise.
  */
int Cluster::orient_points(Pixel a, Pixel b, Pixel c)
{
    int v = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);
    if (v == 0)
        return 0;           // colinear
    return (v > 0) ? 1 : 2; // clock or counterclock wise
}

/**
   *  Get the roundness of a cluster. - 
   * 
   *  [Finds the longest diagonal of the cluster, then it assumes a circle with the
   *  diameter as long as the longest diagonal of the cluster and returns the
   *  percentage of the pixels that are within this circle]
   *
   *  @param[in]  n index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (float) roundness of a cluster (1 - round, 0 - not round).
  */
float Cluster::get_roundness(int n)
{
    float longestDiagonal = get_longest_diagonal(n) + 1;
    float Area = pow(longestDiagonal, 2);
    float clusterArea = l1_clusters[n].size();

    if (Area == 0)
    {
        return 0;
    }

    return (clusterArea / Area);
}

/**
   *  Get coordinates of the longest diagonal of a cluster. - 
   * 
   *  [The function first of all, finds the convex hull of the given cluster of pixels, then it calculates
   *  the length of the diagonal between every two points on the convex hull and stores the coordinates 
   *  of the points on the longest diagonal. These points are then returned.]
   * 
   *  @param[in]  n index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (vector<Point>) vector containing two Points - starting point of the longest diagonal and
   *  ending point of the longest diagonal.
  */
vector<Point> Cluster::get_coordinates_of_longest_diagonal(int n)
{
    vector<Point> border = get_convex_hull(l1_clusters[n], l1_clusters[n].size());  //the convex hull of the given cluster

    int size = border.size();   // number of points in the convex hull of the cluster
    float diagonal = 0; //curent diagonal
    float longestDiagonal = 0;  //longest diagonal found
    int x_i, y_i, x_f, y_f = 0; //the coordinates of the points on the longest diagonal

    if (border.empty())
    {
        size = l1_clusters[n].size();
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {

                diagonal = sqrt(pow(l1_clusters[n][i].x - l1_clusters[n][j].x, 2) + pow(l1_clusters[n][i].y - l1_clusters[n][j].y, 2));
                if (diagonal > longestDiagonal)
                {
                    longestDiagonal = diagonal;
                    x_i = l1_clusters[n][i].x;
                    y_i = l1_clusters[n][i].y;
                    x_f = l1_clusters[n][j].x;
                    y_f = l1_clusters[n][j].y;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                diagonal = sqrt(pow(border[i].x - border[j].x, 2) + pow(border[i].y - border[i].y, 2));
                if (diagonal > longestDiagonal)
                {
                    longestDiagonal = diagonal;
                    x_i = border[i].x;
                    y_i = border[i].y;
                    x_f = border[j].x;
                    y_f = border[j].y;
                }
            }
        }
    }

    vector<Point> p;
    Point p1 = {x_i, y_i};
    Point p2 = {x_f, y_f};

    return p;
}

/**
   *  Get linearity of cluster of pixels. - 
   * 
   *  [Calculates the percentage of the pixels from a given cluster that lay on the longest diagonal of
   *  that cluster]
   * 
   *  @param[in]  n index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (float) percentage of pixels that lay on the longest diagonal.
  */
float Cluster::get_linearity(int n)
{
    int longestDiagonal = get_longest_diagonal(n) + 1;  //the longest diagonal of the given cluster
    float Area = longestDiagonal;   //the length of the longest diagonal (the number of pixels it includes)
    float clusterArea = l1_clusters[n].size();  //the number of pixels in the given cluster

    return (clusterArea / Area);
}

/**
   *  Get the longest diagonal of a cluster of pixels. - 
   * 
   *  [The function first of all, finds the convex hull of the given cluster of pixels, then it calculates
   *  the length of the diagonal between every two points on the convex hull, stores that value and returns
   *  it at the end.]
   * 
   *  @param[in]  n index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (int) Length of the longest diagonal in the given cluster.
  */
int Cluster::get_longest_diagonal(int n)
{
    float A_min = INFINITY;

    vector<Point> border = get_convex_hull(l1_clusters[n], l1_clusters[n].size());  //the convex hull of the given cluster

    int size = border.size();   // number of points in the convex hull of the cluster
    float diagonal = 0; //curent diagonal
    float longestDiagonal = 0;  //longest diagonal found

    if (border.empty())
    {
        size = l1_clusters[n].size();
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {

                diagonal = sqrt(pow(l1_clusters[n][i].x - l1_clusters[n][j].x, 2) + pow(l1_clusters[n][i].y - l1_clusters[n][j].y, 2));
                if (diagonal > longestDiagonal)
                {
                    longestDiagonal = diagonal;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                diagonal = sqrt(pow(border[i].x - border[j].x, 2) + pow(border[i].y - border[i].y, 2));
                if (diagonal > longestDiagonal)
                {
                    longestDiagonal = diagonal;
                }
            }
        }
    }

    return longestDiagonal;
}

/**
   *  Get the time difference between a cluster from layer 1 and layer 2. - 
   * 
   *  [Findes the difference between the size of the first cluster and the second cluster.]
   * 
   *  @return  ( vector<Pair> ) a vector containing Pairs (a structure storing the index of the
   *  cluster in the vector of clusters of the first layer, index of a cluster in the clusters of 
   *  the second layer, the difference in the sizes of the two clusters).
  */
vector<Pair> Cluster::get_delta_size_between_layers()
{
    if (pairs.empty())
    {
        for (int i = 0; i < l1_clusters.size(); i++)
        {
            get_pairs(i);
        }
    }

    vector<Pair> v;
    double dsize;
    Pair p;
    int a = 0;
    int b = 0;

    for (int i = 0; i < pairs.size(); i++)
    {
        if (!pairs[i].empty())
        {
            for (int j = 0; j < pairs[i].size(); j++)
            {
                p.l1 = i;
                p.l2 = pairs[i][j];
                a = l1_clusters[p.l1].size();
                b = l2_clusters[p.l2].size();
                p.dsize = abs(a - b);
                v.push_back(p);
            }
        }
    }
    return v;
}

/**
   *  Loads a clustered file into a given vector of clusters - 
   * 
   *  [Findes the difference between the size of the first cluster and the second cluster.]
   * 
   *  @param[in]  path the absolute path to the clustered file.
   *  @param[in]  set the number of clusters to be stored at once
   *  @param[in]  params the pixel atributes to be read from the source file [x-coordinate, y-coordinate, time, energy]
   *  @param[in]  second signaling whether the source file is for the first layer of the second layer
   *  @return  ( vector<vector<Pixel>> ) a vector containing the loaded Pixels.
  */
vector<vector<Pixel>> Cluster::load_data(std::string path, int set = 500, std::vector<bool> params = {false, false, false, false}, bool second = false)
{
    ifstream File(path);    //source file path
    string line;    //one line from the source file
    std::vector<Pixel> lines;   //all the lines of a single cluster
    string Word;    //one word/information
    vector<vector<Pixel>> c;    //the output vector

    char ch = '0';  //current character
    int i = 0;      //the value of i determines whether the loaded information is x-coordinate, y-coordinate, time or energy
    int X, Y = 0;   //x-coordinate, y-coordinate
    double T = 0;   //time
    float E = 0;    //energy

    if (second == true)
    {
        auto load_all = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            switch (i)
                            {
                            case 0:
                                printf("i is 0");
                                break;
                            case 1:
                                X = 255 - stoi(Word);
                                break;
                            case 2:
                                Y = stoi(Word);
                                break;
                            case 3:
                                T = stod(Word);
                                break;
                            case 4:
                                E = stof(Word);
                                break;
                            default:
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_x = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            X = 255 - stoi(Word);
                            break;

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_y = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            if (i == 2)
                            {
                                Y = stoi(Word);
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_t = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            if (i == 3)
                            {
                                T = stod(Word);
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_e = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            if (i == 4)
                            {
                                Y = stof(Word);
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_xy = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            if (i == 1)
                            {
                                X = 255 - stoi(Word);
                                break;
                            }
                            else if (i == 2)
                            {
                                Y = stoi(Word);
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        if (params[0] && params[1] && params[2] && params[3])
        {
            load_all();
        }
        else if (params[0] && !params[1] && !params[2] && !params[3])
        {
            load_x();
        }
        else if (params[1] && !params[0] && !params[2] && !params[3])
        {
            load_y();
        }
        else if (params[2] && !params[1] && !params[0] && !params[3])
        {
            load_t();
        }
        else if (params[3] && !params[1] && !params[2] && !params[0])
        {
            load_e();
        }
        else if (params[0] && params[1] && !params[2] && !params[3])
        {
            load_xy();
        }
        else
        {
            load_all();
        }
    }
    else
    {
        auto load_all = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            switch (i)
                            {
                            case 0:
                                printf("i is 0");
                                break;
                            case 1:
                                X = stoi(Word);
                                break;
                            case 2:
                                Y = stoi(Word);
                                break;
                            case 3:
                                T = stod(Word);
                                break;
                            case 4:
                                E = stof(Word);
                                break;
                            default:
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_x = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            X = stoi(Word);
                            break;

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_y = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            if (i == 2)
                            {
                                Y = stoi(Word);
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_t = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            if (i == 3)
                            {
                                T = stod(Word);
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_e = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            if (i == 4)
                            {
                                Y = stof(Word);
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        auto load_xy = [&]()
        {
            Word.clear();

            while (std::getline(File, line) && c.size() <= set)
            {
                if (line == "#")
                {
                    c.push_back(lines);
                    lines.clear();
                }
                else
                {
                    for (int j = 0; j < line.size(); j++)
                    {
                        ch = line[j];

                        if (ch == ' ')
                        {
                            i++;

                            if (i == 1)
                            {
                                X = stoi(Word);
                                break;
                            }
                            else if (i == 2)
                            {
                                Y = stoi(Word);
                                break;
                            }

                            Word.clear();
                        }
                        else
                        {
                            Word = Word + ch;
                        }
                    }
                    lines.push_back(Pixel(X, Y, T, E));
                }
                i = 0;
            }
            return c;
        };

        if (params[0] && params[1] && params[2] && params[3])
        {
            load_all();
        }
        else if (params[0] && !params[1] && !params[2] && !params[3])
        {
            load_x();
        }
        else if (params[1] && !params[0] && !params[2] && !params[3])
        {
            load_y();
        }
        else if (params[2] && !params[1] && !params[0] && !params[3])
        {
            load_t();
        }
        else if (params[3] && !params[1] && !params[2] && !params[0])
        {
            load_e();
        }
        else if (params[0] && params[1] && !params[2] && !params[3])
        {
            load_xy();
        }
        else
        {
            load_all();
        }
    }

    return c;
}

/**
   *  Get pairs of clusters form the first layer detector and the second layer detector. - 
   * 
   *  [For each cluster from the first layer the function checks all clusters from the second
   *  layer within a certain range and stores the indicies of those that pass all conditions into
   *  the return vector.]
   * 
   *  @param[in]  v_index the index of a cluster from layer 1 to find the pairs for.
   *  @return  Nothing.
  */
void Cluster::get_pairs(int v_index)
{
    double t1 = l1_clusters[v_index][get_minimum_time_index(v_index, "l1_clusters")].t; //the minimum time from a given cluster in the layer 1
    double t2 = 0;  //the minimum time of the coresponding cluster in layer 2
    int s1 = l1_clusters[v_index].size();   //the size of a given cluster in layer 1
    int s2 = 0; //the size of the coresponding cluster in layer 2
    const int increment = 500;  //the amount of clusters to be loaded at once
    bool ok = false;    //true for pairs that pass all conditions

    pairs.resize((l1_clusters.size()));
    std::vector<bool> params = {true, true, true, false};   //parameters for the loading of the data [x-coordinate, y-coordinate, time, energy]

    int start_index = 0;    //index of the cluster form the layer 2 to start looking for pairs

    l2_clusters = load_data("C:\\Users\\adamk\\OneDrive\\Plocha\\Work\\clustered_data\\office1_f8_px.txt", set, params, true);  //the vector containing all cluster from the layer 2

    if (!l2_clusters.empty())
    {
        for (int i = v_index; i >= 0; i--)
        {
            if (!pairs[i].empty())
            {
                start_index = pairs[i][pairs[i].size() - 1];
                break;
            }
        }

        for (int i = start_index; i < l2_clusters.size(); i++)
        {
            t2 = l2_clusters[i][get_minimum_time_index(i, "l2_clusters")].t;
            s2 = l2_clusters[i].size();

            if (abs(t1 - t2) < 200 && abs(s1 - s2) < 3)
            {
                ok = is_pair(v_index, i);
                if (ok)
                {
                    pairs[v_index].push_back(i);
                }
            }
            else if (t1 - t2 < -201)
            {
                break;
            }
            if (i >= set)
            {
                l2_clusters = load_data("C:\\Users\\adamk\\OneDrive\\Plocha\\Work\\clustered_data\\office1_f8_px.txt", set += increment, params, true);
            }
        }
    }
    else
    {
        printf("File not loaded!");
    }
}

/**
   *  Check if two clusters are pairs. - 
   * 
   *  [Checks if the angle between the two clusters is somewhat similar and if the sizedifference is smaller than 3]
   * 
   *  @param[in]  c1 an index of a cluster from layer 1.
   *  @param[in]  c2 an index of a cluster from layer 2.
   *  @return  true - if the two clusters are pairs; false - if the two clusters aren't pairs.
  */
bool Cluster::is_pair(int c1, int c2)
{
    vector<Pixel> cluster1; //vector of Pixels from first cluster (layer 1)
    vector<Pixel> cluster2; //vector of Pixels from second cluster (layer 2)

    cluster1 = l1_clusters[c1];
    cluster2 = l2_clusters[c2];

    std::sort(cluster1.begin(), cluster1.end());

    Point min_point;    //the Pixel with the lowest time (first detected something) [layer 1]
    Point max_point;    //the Pixel with the greatest time (last detected something) [layer 1]

    if (cluster1.size() >= 6)
    {
        min_point.x = (cluster1[cluster1.size() - 1].x + cluster1[cluster1.size() - 2].x + cluster1[cluster1.size() - 3].x) / 3;
        min_point.y = (cluster1[cluster1.size() - 1].y + cluster1[cluster1.size() - 2].y + cluster1[cluster1.size() - 3].y) / 3;
        max_point.x = (cluster1[0].x + cluster1[1].x + cluster1[2].x) / 3;
        max_point.y = (cluster1[0].y + cluster1[1].y + cluster1[2].y) / 3;
    }
    else if (cluster1.size() >= 4)
    {
        min_point.x = (cluster1[cluster1.size() - 1].x + cluster1[cluster1.size() - 2].x) / 2;
        min_point.y = (cluster1[cluster1.size() - 1].y + cluster1[cluster1.size() - 2].y) / 2;
        max_point.x = (cluster1[0].x + cluster1[1].x) / 2;
        max_point.y = (cluster1[0].y + cluster1[1].y) / 2;
    }
    else
    {
        min_point.x = cluster1[cluster1.size() - 1].x;
        min_point.y = cluster1[cluster1.size() - 1].y;
        max_point.x = cluster1[0].x;
        max_point.y = cluster1[0].y;
    }

    Vect direction = min_point.get_vector(max_point);   //A vector connecting the maxPoint and minPoint [layer 1]

    std::sort(cluster2.begin(), cluster2.end());

    Point min_point2;   //the Pixel with the lowest time (first detected something) [layer 2]
    Point max_point2;   //the Pixel with the greatest time (last detected something) [layer2]

    if (cluster2.size() >= 6)
    {
        min_point2.x = (cluster2[cluster2.size() - 1].x + cluster2[cluster2.size() - 2].x + cluster2[cluster2.size() - 3].x) / 3;
        min_point2.y = (cluster2[cluster2.size() - 1].y + cluster2[cluster2.size() - 2].y + cluster2[cluster2.size() - 3].y) / 3;
        max_point2.x = (cluster2[0].x + cluster2[1].x + cluster2[2].x) / 3;
        max_point2.y = (cluster2[0].y + cluster2[1].y + cluster2[2].y) / 3;
    }
    else if (cluster2.size() >= 4)
    {
        min_point2.x = (cluster2[cluster2.size() - 1].x + cluster2[cluster2.size() - 2].x) / 2;
        min_point2.y = (cluster2[cluster2.size() - 1].y + cluster2[cluster2.size() - 2].y) / 2;
        max_point2.x = (cluster2[0].x + cluster2[1].x) / 2;
        max_point2.y = (cluster2[0].y + cluster2[1].y) / 2;
    }
    else
    {
        min_point2.x = cluster2[cluster2.size() - 1].x;
        min_point2.y = cluster2[cluster2.size() - 1].y;
        max_point2.x = cluster2[0].x;
        max_point2.y = cluster2[0].y;
    }

    Vect direction2 = min_point2.get_vector(max_point2);    //A vector connecting the maxPoint and minPoint [layer 2]

    float dot_product = direction.dot_product(direction2.v);    //Dot product of the two vectors

    float angle = direction.get_angle(direction2.v);    //The angle between the two angles

    if (angle != angle)
    {
        min_point.x = cluster1[cluster1.size() - 1].x;
        min_point.y = cluster1[cluster1.size() - 1].y;
        max_point.x = cluster1[0].x;
        max_point.y = cluster1[0].y;

        min_point2.x = cluster2[cluster2.size() - 1].x;
        min_point2.y = cluster2[cluster2.size() - 1].y;
        max_point2.x = cluster2[0].x;
        max_point2.y = cluster2[0].y;

        direction = min_point.get_vector(max_point);
        direction2 = min_point2.get_vector(max_point2);
        dot_product = direction.dot_product(direction2.v);
        angle = direction.get_angle(direction2.v);
    }

    if (angle > PI / 6 || angle < -PI / 6)
    {
        return false;
    }

    if (cluster1.size() == 1 || cluster2.size() == 1)
    {
        if (cluster1.size() < cluster2.size())
        {
            for (int i = 0; i < cluster1.size(); i++)
            {
                if (sqrt(pow(cluster1[i].x, 2) + pow(cluster2[i].x, 2)) > 10 || sqrt(pow(cluster1[i].y, 2) + pow(cluster2[i].y, 2)) > 10)
                {
                    return false;
                    break;
                }
            }
            return true;
        }
        else
        {
            for (int i = 0; i < cluster1.size(); i++)
            {
                if (sqrt(pow(cluster1[i].x, 2) + pow(cluster2[i].x, 2)) > 10 || sqrt(pow(cluster1[i].y, 2) + pow(cluster2[i].y, 2)) > 10)
                {
                    return false;
                    break;
                }
            }
            return true;
        }
    }

    for (int i = 0; i < cluster2.size(); i++)
    {
        if (direction.x < 0 && cluster2[i].x > max_point.x)
        {
            return false;
            break;
        }
        else if (direction.x > 0 && cluster2[i].x < max_point.x)
        {
            return false;
            break;
        }

        if (direction.y < 0 && cluster2[i].y > max_point.y)
        {
            return false;
            break;
        }
        else if (direction.y > 0 && cluster2[i].y < max_point.y)
        {
            return false;
            break;
        }
    }

    return true;
}

/**
   *  Get the angle of impact of a particle. - 
   * 
   *  [Finds the vector from the pixel that impacted first to the pixel that impacted last.
   *   Afterwards, takes the dot product with the xy-plane and thus finds the impact angle.]
   * 
   *  @param[in]  c1_min a pixel from a cluster with the smallest time from layer 1.
   *  @param[in]  c1_max a pixel from a cluster with the greatest time from layer 1.
   *  @param[in]  c2_min a pixel from a cluster with the smallest time from layer 2.
   *  @param[in]  c2_max a pixel from a cluster with the greatest time from layer 2.
   *  @return  (float) angle of impact.
  */
float get_impact_angle(Pixel c1_min, Pixel c1_max, Pixel c2_min, Pixel c2_max)
{
    // 260 micrometers is the average size of pixel
    std::vector<int> temp = {(c1_max.x - c1_min.x) * 260, (c1_max.y - c1_min.y) * 260, 500}; //a vector containing the x size, y size and height

    Vect d1(temp);   //create a vector object from that vector

    std::vector<int> plane = {0, 0, 1}; //a vector representing the xy plane

    float angle1 = acos((d1.dot_product(plane)) / (sqrt(pow(d1.x, 2) + pow(d1.y, 2) + pow(d1.z, 2))));   //the angle is equal to the dot product / the magnitude

    temp = {(c2_max.x - c2_min.x) * 260, (c2_max.y - c2_min.y) * 260, 500};
    Vect d2(temp);
    float angle2 = acos((d2.dot_product(plane)) / (sqrt(pow(d2.x, 2) + pow(d2.y, 2) + pow(d2.z, 2))));
    
    float angle = (angle1 +  angle2)/2;

    return angle;
}

/**
   *  Gets the distance between two detector layers. - 
   * 
   *  [Finds the pixels with smallest and greatest times from both layers, then calls a function to compute the angle
   *  of impact, which it then uses to compute the distance between the two layers.]
   * 
   *  @param[in]  number index of the cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (float) distance between the two detector layers.
  */
float Cluster::get_distance(int number)
{
    float angle = 0;    //angle of impact of a particle
    float distance = -1;    //distance between the two detector layers
    if (!pairs.empty())
    {
        for (int i = 0; i < pairs[number].size(); i++)
        {
            std::vector<Pixel> t = l1_clusters[number];   //a vector with all the pixels

            std::sort(t.begin(), t.end());

            Pixel c1_min = t[t.size() - 1];
            Pixel c1_max = t[0];

            t = l2_clusters[pairs[number][i]];

            std::sort(t.begin(), t.end());

            Pixel c2_min = t[t.size() - 1];
            Pixel c2_max = t[0];

            angle = get_impact_angle(c1_min, c1_max, c2_min, c2_max);

            return (tan(angle) * sqrt((abs(c1_min.x - c2_max.x) * 260) ^ 2 + (abs(c1_min.y - c2_max.y) * 260) ^ 2));
        }
    }

    return distance;
}

/**
   *  Gets the angle of incidence of a particle. - 
   * 
   *  [Takes in two points - a and b. It first finds the gradient of the line between these two points in the first
   *   layer. Then it calculates the y-intercept to find the definition of the line connecting the two points. So,
   *   two lines are found, then their intercept is found, which represents a third point. Afterwards, all side of
   *   the triangle these there line form are found and thus also the desired angle. This angle is used to find the
   *   height (distance between the two layers).]
   * 
   *  @param[in]  a1 a Point representing the pixel with the lowest time of cluster 1.
   *  @param[in]  a2 a Point representing the pixel with the greatest time of cluster 1.
   *  @param[in]  b1 a Point representing the pixel with the lowest time of cluster 2.
   *  @param[in]  b2 a Point representing the pixel with the greatest time of cluster 2.
   *  @return  (float) the angle of incidence.
  */
float Cluster::get_angle_of_incidence(Point a1, Point a2, Point b1, Point b2)
{
    bool switcha, switchb = false;
    //If the x-coordinate of the second point is not greater than that of the first one, swap them.
    if (a2.x < a1.x)
    {
        int a = a1.x;
        a1.x = a2.x;
        a2.x = a;

        a = a1.y;
        a1.y = a2.y;
        a2.y = a;
        switcha = true;
    }

    if (b2.x < b1.x)
    {
        int b = b1.x;
        b1.x = b2.x;
        b2.x = b;

        b = b1.y;
        b1.y = b2.y;
        b2.y = b;
        switchb = true;
    }

    //Calculating the gradients m1, m2
    float m1 = (a2.y - a1.y) / (a2.x - a1.x);
    float m2 = (b2.y - b1.y) / (b2.x - b1.x);

    //Calculating the intercepts with y-axis
    float c1 = a2.y - m1 * a2.x;
    float c2 = b2.y - m2 * b2.x;

    //find the coordinates of the intercept
    Point P_i;
    P_i.x = (c2 - c1) / (m1 - m2);
    P_i.y = m1 * P_i.x + c1;

    Point P_a, P_b;
    if (switcha)
    {
        P_a.x = a2.x;
        P_a.y = a2.y;
    }
    else
    {
        P_a.x = a1.x;
        P_a.y = a1.y;
    }

    if (switchb)
    {
        P_b.x = b2.x;
        P_b.y = b2.y;
    }
    else
    {
        P_b.x = b1.x;
        P_b.y = b1.y;
    }

    float i = sqrt(pow(P_a.x - P_b.x, 2) + pow(P_a.y - P_b.y, 2));
    float a = sqrt(pow(P_i.x - P_b.x, 2) + pow(P_i.y - P_b.y, 2));
    float b = sqrt(pow(P_i.x - P_a.x, 2) + pow(P_i.y - P_a.y, 2));

    float u = (-pow(i, 2) + pow(a, 2) + pow(b, 2));
    float d = (2 * a * b);
    float f = u / d;

    float angle = acos(f);

    //angle ratio based on the base of a triangle
    float d1 = abs(P_i.x - P_a.x);
    float d2 = abs(P_i.x - P_b.x);
    float r1 = (d1) / (d1 + d2);
    float r2 = (d2) / (d1 + d2);
    float angle1 = r1 * angle;
    float angle2 = r2 * angle;

    float height1 = (d1) / tan(angle1);
    float height2 = (d2) / tan(angle2);

    return height1;
}

/**
   *  Gets the height between two detector layers. - 
   * 
   *  [Finds the pixels with smallest and greatest times for two pairs of clusters, which it then passes on to
   *   a function that calculates the angle of incidence from which this function calculates the height between
   *   the two detectors.]
   * 
   *  @param[in]  c1 index of a first cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @param[in]  c2 index of a second cluster in the vector of clusters @p [vector<vector<Pixel*>> clusters].
   *  @return  (float) distance height the two detector layers.
  */
float Cluster::get_height(int c1, int c2)
{
    vector<Pixel> v = l1_clusters[c1];
    Point c1_i1;    //cluster 1 layer 1
    Point c2_i2;    //cluster 1 layer 2
    Point c1_i2;    //cluster 2 layer 1
    Point c2_i1;    //cluster 2 layer 2

    //cluster 1 layer 1
    std::sort(v.begin(), v.end());
    c1_i1.x = v[0].x;
    c1_i1.y = v[0].y;

    //cluster 1 layer 2
    v.clear();
    v = l2_clusters[pairs[c1][0]];
    std::sort(v.begin(), v.end());
    c1_i2.x = v[v.size() - 1].x;
    c1_i2.y = v[v.size() - 1].y;

    //cluster 2 layer 1
    v.clear();
    v = l1_clusters[c2];
    std::sort(v.begin(), v.end());
    c2_i1.x = v[0].x;
    c2_i1.y = v[0].y;

    //cluster 1 layer 2
    v.clear();
    v = l2_clusters[pairs[c2][0]];
    std::sort(v.begin(), v.end());
    c2_i2.x = v[v.size() - 1].x;
    c2_i2.y = v[v.size() - 1].y;

    float height = get_angle_of_incidence(c1_i1, c1_i2, c2_i1, c2_i2);

    return height;
}
