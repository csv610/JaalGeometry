#pragma once

#include "Mesh.hpp"
#include "DelaunayMesh.hpp"
#include "NearestNeighbours.hpp"

class JCircleCover
{
public:
    void setMesh( const JMeshPtr &m);

    vector<JCircle>  getCircles();
    JMeshPtr         getCoverMesh();

private:
    JMeshPtr mesh, covermesh;
    deque<JNodePtr> nodesQ;
    JNodeSequence farNodesQ, nodesSelected;
    JNodeSequence newnodes;
    JFaceSequence newfaces;

    struct NodeDist
    {
       NodeDist( const JNodePtr &v, double l) { vertex = v; dist = l; }
       JNodePtr  vertex;
       double    dist;
     
    };

    int getNewCircle(JCircle &c);
};
