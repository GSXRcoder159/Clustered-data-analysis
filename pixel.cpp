#include "pixel.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

Pixel::Pixel(int x, int y, double t, float e)
{
    this->x = x;
    this->y = y;
    this->t = t;
    this->e = e;
}
