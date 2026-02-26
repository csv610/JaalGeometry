#include "ConvexHull.hpp"

#include "GeomPredicates.hpp"

void Split( const vector<Vertex*> &nodes, Vertex *v0, Vertex *v0, vector<Vertex*> &hullnodes)
{
    vector<Vertex*> lNodes, rNodes;

    Point3D pa = v0->getXYZCoords();
    Point3D pb = v1->getXYZCoords();
    for( int i = 0; i < nodes.size(); i++) {






    }

    JNodeSequence QuickHull2D( const vector<Vertex*>  &nodes)
    {
        //First two extreme points of the pointset..

        Point3D p3d = points[0]->getXYZCoords();

        double xmin = p3d[0];
        double ymin = p3d[1];
        double xmax = p3d[0];
        double ymax = p3d[1];

        Vertex *v0, *v1;

        for( int i = 0; i < nodes.size(); i++)
            nodes[i]->setVisitBit(0);

        v0->setVisitBit(1);
        v1->setVisitBit(1);
        Split(nodes, v0, v1, hullnodes);













    }
