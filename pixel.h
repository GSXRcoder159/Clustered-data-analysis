using namespace std;

//A structer representing a Pixel and all its atributes
struct Pixel
{
    int x = 0;  //x-coordinate
    int y;  //y-coordinate
    double t;   //time
    float e;    //energy

    //an operater for sorting Pixel based on their time in descending order
    bool operator<(const Pixel &other) const
    {
        return t > other.t;
    }

    //constructer
    Pixel(int x, int y, double t, float e);
};