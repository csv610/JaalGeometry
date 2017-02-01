#include "MeshOptBoundaryLayer.hpp"

////////////////////////////////////////////////////////////////////////////////

void JMeshOptBoundaryLayer :: setNormals()
{
    if( mesh == nullptr) return;

    JEdgeSequence boundedges;
    Vec3D vec;

    mesh->getTopology()->getBoundary(boundedges);
    for( const JEdgePtr &edge: boundedges) {
        const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
        const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();
        double dx = p1[0] - p0[0];
        double dy = p1[1] - p0[1];
        double dl = sqrt(dx*dx + dy*dy);
        vec[0] = -dy/dl;
        vec[1] =  dx/dl;
        vec[2] =  0.0;
        edge->setAttribute("ContourNormal", vec);
    }

    vector<JEdgeSequence> boundloops;
    mesh->getTopology()->getBoundary(boundloops);

    JNodeSequence boundnodes;
    for( size_t i = 0; i < boundloops.size(); i++) {
        JEdgeTopology::getChainNodes( boundloops[i], boundnodes);
        int nnodes = boundnodes.size();
        for( int j = 0; j < nnodes; j++) {
            const Point3D &p0 = boundnodes[(j+2)%nnodes]->getXYZCoords();
            const Point3D &p1 = boundnodes[(j+0)%nnodes]->getXYZCoords();
            double dx = p1[0] - p0[0];
            double dy = p1[1] - p0[1];
            double dl = sqrt(dx*dx + dy*dy);
            vec[0] = -dy/dl;
            vec[1] =  dx/dl;
            vec[2] =  0.0;
            boundnodes[(j+1)%nnodes]->setAttribute("ContourNormal", vec);
        }
    }
    normals_available = 1;
}

////////////////////////////////////////////////////////////////////////////////

void JMeshOptBoundaryLayer :: update()
{
    if( mesh == nullptr) return;
    if( !normals_available) setNormals();

    mesh->buildRelations(0,0);

    JNodeSequence boundnodes, vneighs;
    mesh->getTopology()->getBoundary(boundnodes);

    JNodeSet vSet;
    int   loopid = 1;
    for( const JNodePtr &bnode: boundnodes) {
        JNode::getRelations( bnode, vneighs);
        const Point3D &xyz = bnode->getXYZCoords();
        bnode->setAttribute("Constraint", loopid);
        for( const JNodePtr &vtx : vneighs)
            vSet.insert(vtx);
    }

    loopid = 2;
    Vec3D vec;
    JNodePtr vBound;
    Point3D ptail, phead;
    for( const JNodePtr &vi: vSet) {
        if( !vi->isBoundary() ) {
            JNode::getRelations( vi, vneighs);
            int nCount = 0;
            vBound     = nullptr;
            for( const JNodePtr &vj : vneighs) {
                if( vj->isBoundary() ) {
                    nCount++;
                    vBound = vj;
                }
            }
            if( nCount == 1 ) {
                vBound->getAttribute("ContourNormal", vec);
                ptail  = vBound->getXYZCoords();
                phead[0] = ptail[0] - offset*vec[0];
                phead[1] = ptail[1] - offset*vec[1];
                phead[2] = ptail[2] - offset*vec[2];
                vi->setAttribute("PrevPos", phead);

                phead[0] = ptail[0] - offset*vec[0];
                phead[1] = ptail[1] - offset*vec[1];
                phead[2] = ptail[2] - offset*vec[2];
//              vi->setXYZCoords(phead);
//              vi->setAttribute("Constraint", loopid);
                vi->setAttribute("TargetPos", phead);
            }
        }
    }
    /*
        JLaplaceMeshSmoother lap;
        lap.setMesh(mesh);
        lap.smoothAll();
    */
    cout << " Locally injective mapping" << endl;

    lim.setMesh(mesh);
    lim.solve();

}
////////////////////////////////////////////////////////////////////////////////
