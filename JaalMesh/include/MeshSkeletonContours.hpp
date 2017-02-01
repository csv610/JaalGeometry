#pragma once

#include "MeshSkeleton.hpp"
#include "MeshSlicer.hpp"
#include "MeshBoolean.hpp"

struct JSlice
{
    void  clear()
    {
        boundContour.clear();
        mesh.reset();
    }

    bool  isEmpty() const {
        return boundContour.empty();
    }

    bool     active = 1;
    int      id;
    int      branchID;
    Point3D  center;
    JMeshPtr mesh;
    JEdgeSequence boundContour;

    void genMesh();
};

typedef boost::shared_ptr<JSlice> JSlicePtr;

class JMeshSkeletonContours
{
public:
    void setSkeleton( const JMeshSkeletonPtr &sk);

    void  clearAll();

    // Get all the contours in a given branch ..
    vector<JSlicePtr> getSlices(int branchID);

    // Get Slice center of the edge 
    JSlicePtr  getSlice( const JEdgePtr &edge);

    // Get the  first slice in a given branch. A valid slice must not intersect with
    // Any of the slice in other branches..
    JSlicePtr  getFirstSlice( int branch);

    // Get the last slice in a given branch. A valid slice must not intersect with
    // Any of the slice in other branches..
    JSlicePtr  getLastSlice( int  branch);

    // Get all the slices in the model...
    vector<JSlicePtr> getAllSlices();

    JMeshPtr  getContourMesh();

private:
    boost::scoped_ptr<JMeshSlicer>  meshSlicer;
//  boost::scoped_ptr<JMeshBoolean> meshBoolean;
    JMeshSkeletonPtr  meshSkel;
    vector<JSlicePtr> allSlices;

    bool anyIntersection( const JSlicePtr &s);


    int  getClosestContour(const vector<JEdgeSequence> &contours,
                           const Point3D &center);
};
