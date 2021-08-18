#include "vect.h"
#include <vector>
#include <math.h>

using namespace std;

Vect::Vect(std::vector<int> a)
{
    this->v = a;
    this->x = a[0];
    this->y = a[1];
    this->z = a[2];
}

float Vect::dot_product(std::vector<int> b)
{
    int product = 0;
    for (int i = 0; i < this->v.size(); i++)
    {
        product = product + this->v[i] * b[i];
    }
    return product;
}

float Vect::cross_product(std::vector<int> b)
{
    return 1;
}

float Vect::get_angle(std::vector<int> b)
{
    float dot = dot_product(b);
    float mag_a = sqrt(pow(this->x, 2) + pow(this->y, 2));
    float mag_b = sqrt(pow(b[0], 2) + pow(b[1], 2));

    float angle = acos(dot / (mag_a * mag_b));

    return angle;
}