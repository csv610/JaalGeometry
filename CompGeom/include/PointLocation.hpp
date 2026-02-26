#pragma once

#include "Mesh.hpp"

#include <ANN/ANN.h>

using namespace Jaal;

class JPointLocation
{
public:
    JPointLocation();

    void setMesh( const JMeshPtr &m);

    // Returns false, if given point is completely outside the mesh.
    bool   inDomain( const Point3D &p);

    JNodeSequence getNeighbors( const Point3D &p, int k);

    // Returns the face if the given point is completely inside it.
    // If the point is on the boundary of the edge, or at one of the
    // end point, it is also considered inside,

    JEdgePtr searchEdge( const Point3D &p);
    JFacePtr searchFace( const Point3D &p, bool include_boundary);
    JCellPtr searchCell( const Point3D &p, bool include_boundary);

private:
    JMeshPtr mesh;
    ANNpointArray     dataPts;
    ANNpoint          queryPoint;
    ANNidxArray       nnIdx;
    ANNdistArray      dists;
    ANNkd_tree       *kdTree;

    void preprocess();
    JNodePtr getSeedNode( const Point3D &qp);

    JFacePtr getNextFace( const JFacePtr &f, const JNodePtr &vtx);
    JFacePtr searchFace( const JFacePtr &f, const Point2D &p);

    JCellPtr getNextCell( const JCellPtr &c, const JNodePtr &vtx);
    JCellPtr searchCell( const JCellPtr &c, const Point3D &p);
};

