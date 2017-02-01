#include "MeshSkeletonContours.hpp"

///////////////////////////////////////////////////////////////////////////////
void JSlice :: genMesh()
{
    JNodePtr v0 = JNode::newObject();
    v0->setXYZCoords(center);

    JNodeSequence nodes;
    JMeshTopology::getEntitySet( boundContour, nodes);
    nodes.push_back(v0);
    mesh = JMesh::newObject();
    mesh->addObjects(nodes);

    mesh->addObjects( boundContour);

    JFaceSequence faces;
    for( const JEdgePtr &e : boundContour) {
        const JNodePtr &v1 = e->getNodeAt(0);
        const JNodePtr &v2 = e->getNodeAt(1);
        JFacePtr tri = JTriangle::newObject(v0,v1,v2);
        faces.push_back(tri);
        mesh->addObject(tri);
    }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonContours :: setSkeleton( const JMeshSkeletonPtr &sk)
{
    meshSkel = sk;
    if( meshSkel == nullptr) return;
    JMeshPtr mesh = meshSkel->getMesh();  // Cut with the original mesh ...
    meshSlicer.reset( new JMeshSlicer);
    meshSlicer->setMesh(mesh);
}
///////////////////////////////////////////////////////////////////////////////

int JMeshSkeletonContours :: getClosestContour(const vector<JEdgeSequence> &contours,
        const Point3D &center)
{
    if( contours.size() == 1) return 0;

    JNodeSequence nodes;

    double minDist = std::numeric_limits<double>::max();
    int contourID;

    for( size_t i = 0; i < contours.size(); i++) {
        JMeshTopology::getEntitySet(contours[i], nodes);
        for( const JNodePtr &vtx : nodes) {
            const Point3D &xyz = vtx->getXYZCoords();
            double dist = JMath::length(xyz, center);
            if( dist < minDist) {
                minDist = dist;
                contourID = i;
            }
        }
    }
    return contourID;
}

///////////////////////////////////////////////////////////////////////////////
JSlicePtr JMeshSkeletonContours :: getSlice(const JEdgePtr &edge)
{
    JSlicePtr oldSlice;
    int err = edge->getAttribute("Slice", oldSlice);
    if( !err ) {
        if( oldSlice->active) return oldSlice;
        return nullptr;
    }

    JSlicePtr newSlice;
    Vec3D normal;
    Point3D passThru;
    vector<JEdgeSequence> contours;

    JNodeSequence nodes;

    const Point3D  &p0   = edge->getNodeAt(0)->getXYZCoords();
    const Point3D  &p1   = edge->getNodeAt(1)->getXYZCoords();
    JMath::unit_vector(p1, p0, normal);
    passThru[0]  = 0.5*(p0[0] + p1[0]);
    passThru[1]  = 0.5*(p0[1] + p1[1]);
    passThru[2]  = 0.5*(p0[2] + p1[2]);
    contours = meshSlicer->getContours(normal, passThru);
    int cid = getClosestContour( contours, passThru);
    newSlice.reset( new JSlice );
    newSlice->center   = passThru;
    newSlice->boundContour = contours[cid];
    newSlice->genMesh();
    edge->setAttribute("Slice", newSlice);
    return newSlice;
}
///////////////////////////////////////////////////////////////////////////////

JSlicePtr JMeshSkeletonContours :: getFirstSlice(int branchID )
{
    JMeshSkeletonBranchPtr branch = meshSkel->getBranch(branchID);
    if( branch == nullptr) return nullptr;
    JEdgeSequence branchEdges = branch->skelEdges;

    int numEdges = branchEdges.size();
    for( int i = 0; i < numEdges; i++) {
        JSlicePtr newSlice = getSlice( branchEdges[i] );
        if( newSlice ) {
            newSlice->id       = i;
            newSlice->branchID = branchID;
            return newSlice;
        }
    }
    return nullptr;
}
///////////////////////////////////////////////////////////////////////////////
JSlicePtr JMeshSkeletonContours :: getLastSlice(int branchID )
{
    JMeshSkeletonBranchPtr branch = meshSkel->getBranch(branchID);
    if( branch == nullptr) return nullptr;

    JEdgeSequence branchEdges = branch->skelEdges;
    int numEdges = branchEdges.size();

    for( int i = 0; i < numEdges; i++) {
        JSlicePtr newSlice = getSlice( branchEdges[numEdges-1-i] );
        if( newSlice ) {
            newSlice->id       = numEdges-i-1;
            newSlice->branchID = branchID;
            return newSlice;
        }
    }
    return nullptr;
}

///////////////////////////////////////////////////////////////////////////////
vector<JSlicePtr> JMeshSkeletonContours :: getSlices(int branchID )
{
    vector<JSlicePtr> slices;
    JMeshSkeletonBranchPtr branch = meshSkel->getBranch(branchID);
    if( branch == nullptr) return slices;

    JEdgeSequence branchEdges = branch->skelEdges;

    JSlicePtr newSlice;

    size_t numEdges = branchEdges.size();
    for( int i = 0; i < numEdges; i++) {
        newSlice = getSlice( branchEdges[i] );
        if( newSlice ) {
            newSlice->id       = i;
            newSlice->branchID = branchID;
            slices.push_back(newSlice);
        }
    }
    return slices;
}
///////////////////////////////////////////////////////////////////////////////

vector<JSlicePtr> JMeshSkeletonContours :: getAllSlices()
{
    vector<JSlicePtr> branchSlices;
    allSlices.clear();

    if( meshSkel == nullptr) return branchSlices;

    size_t numBranches = meshSkel->getNumBranches();
    if( numBranches < 2) return branchSlices;

    vector<JMeshPtr>  brachMesh(numBranches);

    for( size_t i = 0; i <  numBranches; i++) {
        branchSlices = getSlices( i );
        boost::copy( branchSlices, back_inserter( allSlices) );
    }

    // Now check for intersection of slices. If any slice in a  particular
    // branch intersects with any slice in other branch, then  discard this
    // slice. "meshA" is the particular slice and "meshB" contains all other
    // slices.
    JMeshBoolean meshBool;

    int op = JMeshBoolean::MESH_INTERSECTION;

    JMeshPtr meshA, meshB;
    meshA = JMesh::newObject();
    for( size_t i = 0; i <  numBranches; i++) {
        meshA->clearAll();
        for( auto slice : allSlices) {
            if( slice->branchID != i )
                meshA->addObject( slice->mesh );
        }
        if( meshA->getSize(2) ) {
            meshBool.setMesh1(meshA);

            for( auto &slice : allSlices) {
                if( slice->branchID == i ) {
                    meshB =  slice->mesh;
                    if( meshB->getSize(2) ) {
                        meshBool.setMesh2( meshB);
                        JMeshPtr commMesh = meshBool.getMesh(op);
                        if(commMesh->getSize(0)) slice->active = 0;
                    }
                }
            }
        }
    }

    branchSlices = allSlices;
    allSlices.clear();

    for( auto slice : branchSlices) {
        if( !slice->active)  allSlices.push_back(slice);
    }
    return allSlices;
}
///////////////////////////////////////////////////////////////////////////////

/*
vector<JEdgeSequence> JMeshSkeletonContours :: getBranchContours( int branchID)
{
    if( allSlices.empty() ) getSlices();

    vector<JEdgeSequence> branchContours;
    for( auto slice  : allSlices) {
        if( slice->branchID == branchID ) {
            if( !slice->isEmpty()  )
                branchContours.push_back(slice->boundContour);
        }
    }
    return branchContours;
}
///////////////////////////////////////////////////////////////////////////////

vector<JEdgeSequence> JMeshSkeletonContours :: getContours()
{
    if( allSlices.empty() ) getSlices();

    vector<JEdgeSequence> allContours;
    for( auto slice  : allSlices) {
        if( !slice->isEmpty()  )
            allContours.push_back(slice->boundContour);
    }
    return allContours;
}

JMeshPtr JMeshSkeletonContours :: getDisks()
{
    if( allSlices.empty() ) getSlices();

    JMeshPtr disks = JMesh::newObject();
    for( auto &slice : allSlices) {
        if( slice->mesh ) disks->addObject(slice->mesh);
    }
    return disks;
}
*/
///////////////////////////////////////////////////////////////////////////////

