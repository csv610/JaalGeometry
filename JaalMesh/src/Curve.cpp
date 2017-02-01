#include "Curve.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
JCurve* JCurve :: readFromFile ( const string &filename )
{
    if( filename.empty() ) return NULL;

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() ) return NULL;

    JCurve *newCurve = new JCurve;
    string str;

    infile >> str;
    if( str != "#Nodes") {
        return NULL;
    }

    int numNodes;
    infile >> numNodes;

    Point3D p3d;
    for( int i = 0; i < numNodes; i++) {
        infile >> p3d[0] >> p3d[1] >> p3d[2];
        JNodePtr vtx = JNode::newObject();
        vtx->setXYZCoords(p3d);
        newCurve->addControlNode(vtx);
    }
    return newCurve;
}

///////////////////////////////////////////////////////////////////////////////
double JCurve :: getLength() const
{
    double len = 0.0;
    for( size_t i = 0; i < nodes.size()-1; i++)
        len += JNodeGeometry::getLength(nodes[i], nodes[i+1] );

    if( type == CLOSE_CURVE)
        len += JNodeGeometry::getLength( nodes.back(), nodes.front() );

    return len;
}

///////////////////////////////////////////////////////////////////////////////
void JCurve :: refineSegments(int n)
{
    if( n < 3) return;

    JNodeSequence newnodes, edgenodes;
    int nSize = nodes.size();
    for( int i = 0; i < nSize; i++) {
        JNodePtr v0 = nodes[i];
        JNodePtr v1 = nodes[(i+1)%nSize];
        JEdgeGeometry::generateLinearNodes(v0, v1, n, edgenodes);
        for( size_t j = 0; j < edgenodes.size()-1; j++)
            newnodes.push_back( edgenodes[j] );
    }
    nodes = newnodes;
}
///////////////////////////////////////////////////////////////////////////////

double JCurve :: getMinimumLocalDistance() const
{
    double dist, minDist = 0.0;

    if(ctrlnodes.size() < 2 ) return 0.0;
    minDist = JNodeGeometry::getLength( ctrlnodes[0], ctrlnodes[1] );

    int nSize = ctrlnodes.size();
    for( int  i = 1; i < nSize-1; i++) {
        dist  = JNodeGeometry::getLength( ctrlnodes[i], ctrlnodes[i+1] );
        minDist = min( minDist, dist );
    }

    if( type == CLOSE_CURVE ) {
        dist  = JNodeGeometry::getLength( ctrlnodes.back(), ctrlnodes.front() );
        minDist = min( minDist, dist );
    }

    return minDist;
}

///////////////////////////////////////////////////////////////////////////////
double JCurve :: getMinimumGlobalDistance() const
{
    double dist, minDist = 0.0;

    if(ctrlnodes.size() < 2 ) return 0.0;
    minDist = JNodeGeometry::getLength( ctrlnodes[0], ctrlnodes[1] );

    int nSize = ctrlnodes.size();
    for( int i = 0;   i < nSize; i++) {
        for( int j = i+1; j < nSize; j++) {
            dist  = JNodeGeometry::getLength( ctrlnodes[i], ctrlnodes[j] );
            minDist = min( minDist, dist );
        }
    }
    return minDist;
}
///////////////////////////////////////////////////////////////////////////////

int JCurve :: getTangentAt(const JNodePtr vtx, Vec3D &vec) const
{
    int pos = -1;

    size_t nSize = nodes.size();
    for( size_t i = 0; i < nSize; i++) {
        if( nodes[i] == vtx ) {
            pos = i;
        }
    }

    if( pos == -1) return 1;

    Point3D pprev, pnext, pcurr;

    if( pos == 0) {
        if( type == OPEN_CURVE) {
            pcurr = nodes[0]->getXYZCoords();
            pnext = nodes[1]->getXYZCoords();
            vec[0] = pnext[0] - pcurr[0];
            vec[1] = pnext[1] - pcurr[1];
            vec[2] = pnext[2] - pcurr[2];
        } else {
            pprev  = nodes[nSize-1]->getXYZCoords();
            pnext  = nodes[1]->getXYZCoords();
            vec[0] = 0.5*( pnext[0] - pprev[0]);
            vec[1] = 0.5*( pnext[1] - pprev[1]);
            vec[2] = 0.5*( pnext[2] - pprev[2]);
        }
        return 0;
    }

    if( size_t(pos) == nSize-1) {
        if( type == OPEN_CURVE) {
            pnext = nodes[nSize-1]->getXYZCoords();
            pcurr = nodes[nSize-2]->getXYZCoords();
            vec[0] = pnext[0] - pcurr[0];
            vec[1] = pnext[1] - pcurr[1];
            vec[2] = pnext[2] - pcurr[2];
        } else {
            pprev  = nodes[nSize-2]->getXYZCoords();
            pnext  = nodes[0]->getXYZCoords();
            vec[0] = 0.5*( pnext[0] - pprev[0]);
            vec[1] = 0.5*( pnext[1] - pprev[1]);
            vec[2] = 0.5*( pnext[2] - pprev[2]);
        }
        return 0;
    }

    pprev  = nodes[pos-1]->getXYZCoords();
    pnext  = nodes[pos+1]->getXYZCoords();
    vec[0] = 0.5*( pnext[0] - pprev[0]);
    vec[1] = 0.5*( pnext[1] - pprev[1]);
    vec[2] = 0.5*( pnext[2] - pprev[2]);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void JCurve :: deleteAll()
{
    cout << "EXIT" << endl;
    exit(0);
    /*
         size_t nSize = edges.size();
         for( size_t i = 0; i < nSize; i++)
              delete edges[i];
         edges.clear();

         nSize = nodes.size();
         for( size_t i = 0; i < nSize; i++)
              delete nodes[i];
         nodes.clear();
    */
}

