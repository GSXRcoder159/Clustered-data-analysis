#include "point.h"
#include <iostream>
#include <vector>

using namespace std;
#define INF 100000

Point::Point(int x, int y)
{
    this->x = x;
    this->y = y;
}

int Point::orient(Point a, Point b, Point c)
{
    int v = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);
    if (v == 0)
        return 0;           // colinear
    return (v > 0) ? 1 : 2; // clock or counterclock wise
}

void Point::convexHull(Point points[], int m)
{
    if (m < 3) //at least three points required
        return;
    int n[m];
    for (int i = 0; i < m; i++)
        n[i] = -1;
    int l = 0; //initialize result.
    for (int i = 1; i < m; i++)
        if (points[i].x < points[l].x)
            l = i; //find left most point
    int p = l, q;
    do
    {
        q = (p + 1) % m;
        for (int i = 0; i < m; i++)
            if (orient(points[p], points[i], points[q]) == 2)
                q = i;
        n[p] = q;
        p = q;
    } while (p != l);
    for (int i = 0; i < m; i++)
    {
        if (n[i] != -1)
            cout << "(" << points[i].x << ", " << points[i].y << ")\n";
    }
}

std::vector<int> Point::get_vector(Point a)
{
    vector<int> v;

    v.push_back(a.x - this->x);
    v.push_back(a.y - this->y);

    return v;
}