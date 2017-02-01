
#ifdef CSV
#include "Tube.hpp"

void Tube :: setSegments( const JEdgeSequence &edges)
{
    clear();

    if( edges.empty() ) return;

    numEdges = edges.size();
    JNodePtr vtx;
    Point3D xyz;

    bool closed = 0;
    if( edges.front()->getNodeAt(0) == edges.back()->getNodeAt(1) ) {
        closed = 1;
    } else {
        closed = 0;
    }

    numPoints = numEdges + 3;
    points = new RowArray[numPoints];

    JEdgePtr edge;
    if( !closed ) {
        edge =  edges.front();
        vtx = edge->getNodeAt(1);
        xyz = vtx->getXYZCoords();
        points[0][0] = xyz[0];
        points[0][1] = xyz[1];
        points[0][2] = xyz[2];

        for( size_t i = 0; i < numEdges; i++) {
            vtx = edges[i]->getNodeAt(0);
            xyz = vtx->getXYZCoords();
            points[i+1][0] = xyz[0];
            points[i+1][1] = xyz[1];
            points[i+1][2] = xyz[2];
        }

        edge =  edges.back();
        vtx = edge->getNodeAt(1);
        xyz = vtx->getXYZCoords();
        points[numEdges+1][0] = xyz[0];
        points[numEdges+1][1] = xyz[1];
        points[numEdges+1][2] = xyz[2];

        edge =  edges.back();
        vtx = edge->getNodeAt(0);
        xyz = vtx->getXYZCoords();
        points[numEdges+2][0] = xyz[0];
        points[numEdges+2][1] = xyz[1];
        points[numEdges+2][2] = xyz[2];
    }

    if( closed ) {
        edge =  edges[numEdges-2];
        vtx  = edge->getNodeAt(0);
        xyz  = vtx->getXYZCoords();
        points[0][0] = xyz[0];
        points[0][1] = xyz[1];
        points[0][2] = xyz[2];

        edge =  edges[numEdges-1];
        vtx  = edge->getNodeAt(0);
        xyz  = vtx->getXYZCoords();
        points[1][0] = xyz[0];
        points[1][1] = xyz[1];
        points[1][2] = xyz[2];

        for( size_t i = 0; i < numEdges; i++) {
            vtx = edges[i]->getNodeAt(0);
            xyz = vtx->getXYZCoords();
            points[i+2][0] = xyz[0];
            points[i+2][1] = xyz[1];
            points[i+2][2] = xyz[2];
        }
        edge =  edges.front();
        vtx = edge->getNodeAt(1);
        xyz = vtx->getXYZCoords();
        points[numEdges+2][0] = xyz[0];
        points[numEdges+2][1] = xyz[1];
        points[numEdges+2][2] = xyz[2];
    }
}
////////////////////////////////////////////////////////////////////////////////

void Tube :: draw()
{
    glColor3f( color[0], color[1], color[2]  );

    glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, 120);
    glDisable( GL_CULL_FACE);

    glEnable( GL_COLOR_MATERIAL );
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glShadeModel(GL_SMOOTH );
    gleSetJoinStyle (TUBE_NORM_FACET | TUBE_NORM_EDGE | TUBE_JN_ANGLE | TUBE_JN_CAP);
//   gleSetJoinStyle (TUBE_JN_CAP);

    gleSetJoinStyle (TUBE_NORM_FACET);
    gleSetJoinStyle (TUBE_JN_ROUND );
    gleSetNumSides(numSides);
    glPushMatrix();
    glePolyCylinder (numPoints, points, nullptr, radius);
    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
#endif
