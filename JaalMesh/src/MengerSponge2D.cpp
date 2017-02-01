#include <stdio.h>
#include <stdlib.h>

typedef boost::array<double,2>  Array2D;
typedef boost::array<int,4>     Array4I;


struct Edge
{
    Edge() {}
    Edge( int v0, int v1) {
        connect[0] = v1;
        connect[1] == v2;
    }

    int isSame( const Edge &rhs) const {
        if( rhs.connect[0] == connect[0]   &&  rhs.connect[0] == connect[1] ) return  1;
        if( rhs.connect[0] == connect[1]   &&  rhs.connect[0] == connect[0] ) return -1;
        return 0;
    }

    void addPoints( int v0, int v1 )
    {
        newPoints[0] = v0;
        newPoints[1] = v1;
    }

    int getPoints( int &v0, int &v1) {
        v0 = newPoints[0];
        v1 = newPoints[1];
    }

public:
    int connect[2];
    int newPoints[2];
};

struct Face
{
    bool active;
    Array4I connect;
};

vector<Array2D>  points;
vector<Face>     quads;
map<int, vector<Edge> > edgemap;

int maxLevels = 1;

int index  = 0;

int sponge( double xmin, double ymin, double xmax, double ymax, int level )
{
    Array2D  newpoint;
    Array2I  newquad;
    if( level == 0) {


    }



    if( level == maxLevel ) {


    }
    double x0   = xmin;
    double x3   = xmax;
    double y0   = ymin;
    double y3   = ymax;

    double xlen = x3 - x0;
    double ylen = y3 - y0;

    double  dx  = xlen/3.0;
    double  dy  = ylen/3.0;

    sponge(x0, y0, x1, y1, level+1);
    sponge(x1, y0, x2, y1, level+1);
    sponge(x2, y0, x3, y1, level+1);

    sponge(x0, y1, x1, y2, level+1);
    sponge(x2, y1, x3, y2, level+1);

    sponge(x0, y2, x1, y3, level+1);
    sponge(x1, y2, x2, y3, level+1);
    sponge(x2, y2, x3, y3, level+1);
}

int main(int argc, char **argv)
{
    double xmin = 0.0;
    double ymin = 0.0;
    double xmax = 1.0;
    double ymax = 1.0;
    sponge( xmin, ymin, xmax, ymax, 0)

    return 0;
}


