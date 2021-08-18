#include "line.h"
#include <math.h>
#include <vector>

#define PI 3.14159265

line::line(int x, int y, float angle)
{
    this->x = x;
    this->y = y;
    this->angle = angle;
}

int line::get_y(int x)
{
    return (tan(angle * PI / 180.0) * (x - this->x) + this->y);
}

int line::get_x(int y)
{
    return (y / (tan(angle * PI / 180.0)) - this->y / (tan(angle * PI / 180.0)) + this->x);
}

vector<int> line::get_x_intercept(line l1, line l2)
{
    vector<int> v;
    bool intercept = false;
    float x = l1.x;
    float y = 0;
    float dx = 0.1;

    for (int i = 0; i < 10000; i++)
    {
        if (l1.get_y(x) == l2.get_y(x))
        {
            y = l1.get_y(x);
            intercept = true;
            break;
        }
        x -= dx;
    }

    if (y != 0)
    {
        v.push_back((int)x);
        v.push_back((int)y);
    }

    return v;
}

vector<int> line::get_y_intercept(line l1, line l2)
{
    vector<int> v;
    bool intercept = false;
    float x = 0;
    float y = l1.y;
    float dy = 0.1;

    for (int i = 0; i < 10000; i++)
    {
        if (l1.get_x(y) == l2.get_x(y))
        {
            x = l1.get_x(y);
            intercept = true;
            break;
        }
        y -= dy;
    }

    if (y != 0)
    {
        v.push_back((int)x);
        v.push_back((int)y);
    }

    return v;
}