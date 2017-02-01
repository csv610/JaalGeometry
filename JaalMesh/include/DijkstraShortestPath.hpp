#pragma once 

#include "Mesh.hpp"
#include "basic_math.hpp"
#include <queue>

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////

class  JDijkstraShortestPath 
{
     struct Compare
     {
        bool operator() (const JNodePtr &v0, const JNodePtr &v1) const
        {
             NodeDist d0, d1;
             v0->getAttribute("NodeDist", d0);
             v1->getAttribute("NodeDist", d1);
             return d0.first > d1.first;
        }
      };
    std::priority_queue<JNodePtr, JNodeSequence, Compare> nodesQ;

public:

    void setMesh(JMeshPtr &m);
    void setMesh(JMeshPtr &m, const JMeshFilterPtr &f);

    JNodeSequence getPath( const JNodePtr &vs, const JNodePtr &vd);

    void setDistance( const JNodePtr &vs);
    void setDistance( const JNodeSequence &vs);

private:
    // Input parameters ...
    JMeshPtr mesh;
    JNodePtr vdst;
    JMeshFilterPtr  filter;

    typedef std::pair<double,JNodePtr> NodeDist;

    JNodeSequence  sourceNodes;

    // Output data ...
    JNodeSequence  nodepath;

    double getCost(const JNodePtr &vi, const JNodePtr &vj ) const;

    bool  atomicOp( const JNodePtr &node);
    void  initialize();
    void  fastmarching();  // Fast Marching Style algorithm O(nlogn) with heap
    void  traceback();
};

