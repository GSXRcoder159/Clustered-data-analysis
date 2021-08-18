#include <vector>
#include "pixel.h"  //Pixel
#include "point.h"  //Point
#include "pair.h"   //Pair
#include "line.h"   //Line
#include <iostream>
#include <string>
#include <list>
#include "vect.h"   //Vect

using namespace std;

struct Cluster
{
    std::vector<Pixel> lines;

    std::vector<vector<Pixel>> l1_clusters; //clusters from the first detector layer
    std::vector<vector<Pixel>> l2_clusters; //clusters form the second detector layer
    std::vector<vector<int>> pairs; //vector of indicies of pairs of clusters form layer1 and layer 2
    int set = 500;  //the maximum number of pixels to be loaded at once

    string Word;

    Cluster();

    void input(std::string path);

    float get_total_energy(int n);

    double get_minimum_time(int n);

    int get_size(int n);

    Pixel get_pixel_with_minimum_time(int n);

    Pixel get_pixel_with_maximum_time(int n);

    double get_delta_time(int n);

    int get_minimum_time_index(int n, string vector);

    int get_maximum_time_index(int n);


    int get_x_dimension(int n);

    int get_y_dimension(int n);

    float get_roundness(int n);

    float get_linearity(int n);

    int find_longest_diagonal(vector<Point> &points, int n, Point pix);

    int find_longest_diagonal_2(vector<Point> &points, int n, Point pix, std::vector<int> v, vector<Point> &used);

    vector<Pair> get_delta_size_between_layers();

    void output(std::string path, std::string function);

    void exp(std::string path);

    vector<Point> get_convex_hull(vector<Pixel> points, int m);

    int orient_points(Pixel a, Pixel b, Pixel c);

    int get_longest_diagonal(int n);

    vector<Point> get_coordinates_of_longest_diagonal(int n);

    vector<vector<Pixel>> load_data(std::string path, int set, std::vector<bool> params, bool second);

    void get_pairs(int v_index);

    bool is_pair(int c1, int c2);

    float get_distance(int number);

    float get_angle_of_incidence(Point a1, Point a2, Point b1, Point b2);

    float get_height(int c1, int c2);
};