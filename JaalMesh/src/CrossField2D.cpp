#include "CrossField2D.hpp"

//////////////////////////////////////////////////////////////////////////////////
void JCrossField2D :: setMesh( const JMeshPtr &m)
{
    mesh = m;
}

//////////////////////////////////////////////////////////////////////////////////
void JCrossField2D :: genVecField( )
{
    // Idea from : A frontal Delaunay quad mesh generator using the Linf norm
    //             Remacle etc ..

    if( mesh == nullptr) return;
    mesh->pruneAll();

    vector<JEdgeSequence>   boundEdges;
    JNodeSequence   boundNodes;

    mesh->getTopology()->getBoundary(boundEdges);
    mesh->getTopology()->getBoundary(boundNodes);

    Vec3D normVec;
    normVec[0] = 0.0;
    normVec[1] = 0.0;
    normVec[2] = 0.0;

    for( const JNodePtr &vtx : boundNodes )
        vtx->setAttribute("Normal", normVec);

    int numloops = boundEdges.size();
    for( int i = 0; i < numloops; i++) {
        if( JEdgeGeometry::getOrientation( boundEdges[i] )  < 0.0)
            JEdgeTopology::reverse( boundEdges[i] );

        int numedges = boundEdges[i].size();
        for( size_t j = 0; j < numedges; j++) {
            const JNodePtr &v0 = boundEdges[i][j]->getNodeAt(0);
            const JNodePtr &v1 = boundEdges[i][j]->getNodeAt(1);
            const Point3D  &p0 = v0->getXYZCoords();
            const Point3D  &p1 = v1->getXYZCoords();
            double dx = p0[0] - p1[0];
            double dy = p1[1] - p0[1];
            v0->getAttribute("Normal", normVec);
            normVec[0] += dy;
            normVec[1] += dx;
            v0->setAttribute("Normal", normVec);

            v1->getAttribute("Normal", normVec);
            normVec[0] += dy;
            normVec[1] += dx;
            v1->setAttribute("Normal", normVec);
        }
    }

    for( const JNodePtr &vtx : boundNodes ) {
        vtx->getAttribute("Normal", normVec);
        double  dx = normVec[0];
        double  dy = normVec[1];
        double  dl  = sqrt(dx*dx + dy*dy);
        normVec[0] /= dl;
        normVec[1] /= dl;
        vtx->setAttribute("Normal", normVec);
    }

    map<JNodePtr, double> boundAngle;
    for( const JNodePtr &vtx : boundNodes ) {
        vtx->getAttribute("Normal", normVec);
        double dx = normVec[0];
        double dy = normVec[1];
        boundAngle[vtx] = atan2(dy,dx);
    }

    int numnodes = mesh->getSize(0);
    sin4t.resize(numnodes);
    cos4t.resize(numnodes);
    theta.resize(numnodes);

    JHarmonicField   hfield;
    hfield.setMesh(mesh);
    hfield.setFieldName("Dirichlet");

    // Solver for "cos(4t) " field ...
    double val;
    for( const JNodePtr &vtx : boundNodes ) {
        double t = boundAngle[vtx];
        val = cos(4*t);
        vtx->setAttribute("Dirichlet", val);
    }
    hfield.solveSystem();

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->getAttribute("Dirichlet", val);
        cos4t[i] = val;
    }

    // Solver for "sin(4t)" field ...
    for( const JNodePtr &vtx : boundNodes ) {
        double t = boundAngle[vtx];
        val = sin(4*t);
        vtx->setAttribute("Dirichlet", val);
    }
    hfield.solveSystem();
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->getAttribute("Dirichlet", val);
        sin4t[i] = val;
    }

    // Get Angle at internal nodes ...
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        theta[i] = 0.25*atan2(sin4t[i], cos4t[i]);
        vtx->setAttribute("CrossFieldAngle", theta[i] );
    }

    // Set the edge orientation ...
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        int  v1  = edge->getNodeAt(0)->getID();
        int  v2  = edge->getNodeAt(1)->getID();
        double t   = 0.5*(theta[v1] + theta[v2]);
        edge->setAttribute("CrossFieldAngle", t);
    }

    mesh->deleteNodeAttribute("Dirichlet");
}

//////////////////////////////////////////////////////////////////////////////////
