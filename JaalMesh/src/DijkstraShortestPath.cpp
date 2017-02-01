#include "StopWatch.hpp"
#include "DijkstraShortestPath.hpp"

///////////////////////////////////////////////////////////////////////////////

void JDijkstraShortestPath:: setMesh(JMeshPtr &m) {
        mesh = m;
        vdst = nullptr;
        filter = nullptr;
}

void JDijkstraShortestPath:: setMesh(JMeshPtr &m, const JMeshFilterPtr &f) {
        mesh = m;
        vdst = nullptr;
        filter = f;
}

bool JDijkstraShortestPath::atomicOp( const JNodePtr &vi)
{
    if (vi == vdst) return 0;

/*
    if( filter ) {
        if( vi != vsrc  ) {
            if( !filter->passThrough( vi ) ) {
                vdst = vi;
                return 0;
            }
        }
    }
*/

    NodeDist viAttrib, vjAttrib;
    vi->getAttribute("NodeDist", viAttrib);
    double di = viAttrib.first;

    JNodeSequence vneighs;
    JNode::getRelations( vi, vneighs );
    int nSize = vneighs.size();

    //   int offset = JMath::random_value((int)0, nSize-1);
    int offset = 0;

    for (int i = 0; i < nSize; i++) {
        JNodePtr vj = vneighs[(i+offset)%nSize];
//      double   dij = JNodeGeometry::getLength(vi,vj);
        double   dij = 1.0;

        vj->getAttribute("NodeDist", vjAttrib);
        double dj = vjAttrib.first;

        if (dj > di + dij) {
            vjAttrib.first   = di + dij;
            vjAttrib.second  = vi;
            nodesQ.push(vj);
            vj->setAttribute("NodeDist", vjAttrib);
        }
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void JDijkstraShortestPath::traceback()
{
   nodepath.clear();
   if( vdst == nullptr) return;

   NodeDist nodedist;
   JNodePtr currnode = vdst;
   while(1) {
        if( currnode == nullptr) break;
        nodepath.push_back(currnode);
        currnode->getAttribute("NodeDist", nodedist);
        currnode = nodedist.second;
   }
}
///////////////////////////////////////////////////////////////////////////////

void JDijkstraShortestPath :: initialize()
{
    if( mesh == nullptr) return;

    size_t numNodes = mesh->getSize(0);

    NodeDist nodedist;
    nodedist.first = std::numeric_limits<double>::max();
    nodedist.second= nullptr;

    for( size_t i = 0; i < numNodes; i++) {
         const JNodePtr &vtx = mesh->getNodeAt(i);
         vtx->setAttribute("NodeDist", nodedist);
    }

    while(!nodesQ.empty()) nodesQ.pop(); 

    nodedist.first = 0.0;
    nodedist.second= nullptr;
    for( const JNodePtr &vtx : sourceNodes) {
         vtx->setAttribute("NodeDist", nodedist);
         nodesQ.push(vtx);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JDijkstraShortestPath :: fastmarching()
{
   initialize();

   bool progress;
   while (!nodesQ.empty()) {
        JNodePtr currVertex = nodesQ.top();
        nodesQ.pop();
        progress = atomicOp(currVertex);
        if (!progress) break;
   }
}
///////////////////////////////////////////////////////////////////////////////
JNodeSequence  JDijkstraShortestPath :: getPath( const JNodePtr &s, const JNodePtr &d)
{
    sourceNodes.empty();
    sourceNodes.push_back(s);
    vdst = d;

    fastmarching();
    traceback();
    return nodepath;
}
///////////////////////////////////////////////////////////////////////////////
