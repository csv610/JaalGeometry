#include "DDG_Edge.hpp"
#include "DDG_Mesh.hpp"

namespace DDG
{
Vector Edge :: dualPoint() const
{
    Vector pj = he->flip->vertex->position;
    Vector pi = he->vertex->position;

    Vector eij = (pj - pi).unit();
    double dij = he->shift();
    return pi + dij*eij;
}


}

