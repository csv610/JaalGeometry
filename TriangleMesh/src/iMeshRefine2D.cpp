#include "MeshRefine2D.h"

////////////////////////////////////////////////////////////////////////////////

int MeshRefine2D :: initialize()
{
    int err, result, namelen;
    numfacesRefined = 0;

    insertedNodes.clear();
    insertedFaces.clear();

    const char *tag1 = "VERTEX_ON_EDGE";
    namelen = strlen( tag1 );
    iMesh_createTag(mesh, tag1, 1, iBase_ENTITY_HANDLE, &vertex_on_edge_tag, &result, namelen);

    const char *tag2 = "RemoveFace";
    namelen = strlen( tag2 );
    iMesh_createTag(mesh, tag2, 1, iBase_INTEGER,  &remove_tag, &result, namelen);

    const char *tag3 = "Boundary";
    namelen = strlen( tag3 );
    iMesh_createTag(mesh, tag3, 1, iBase_INTEGER,  &boundary_tag, &result, namelen);

    iMesh_getTagHandle(mesh, "GLOBAL_ID", &globalID_tag, &err, strlen("GLOBAL_ID") );

    iMesh_getRootSet(mesh, &rootSet, &err);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int MeshRefine2D :: finalize()
{
    int err;

    SimpleArray<iBase_EntityHandle> faces;
    iMesh_getEntities( mesh, rootSet, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                       ARRAY_INOUT(faces), &err);
    assert( !err );

    int val;
    for( size_t i = 0; i < faces.size(); i++) {
        iMesh_getIntData( mesh, faces[i], remove_tag, &val, &err);
        if( err == 0 && val == 1 ) {
            iMesh_rmvTag(mesh, faces[i], remove_tag, &err);
            iMesh_deleteEnt( mesh, faces[i], &err);
        }
    }
    iMesh_destroyTag(mesh, remove_tag, 1, &err);

    SimpleArray<iBase_EntityHandle> nodes;
    iMesh_getEntities( mesh, rootSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                       ARRAY_INOUT(nodes), &err);

    for( size_t i = 0; i < nodes.size(); i++)
        iMesh_setIntData(mesh, nodes[i], globalID_tag, i, &err);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool
MeshRefine2D::searchEdge( iBase_EntityHandle v1, iBase_EntityHandle v2,
                          iBase_EntityHandle &edgehandle ) const
{
    int err;
    iBase_EntityHandle vmin = min(v1,v2);
    iBase_EntityHandle vmax = max(v1,v2);

    SimpleArray<iBase_EntityHandle> edges;
    iMesh_getEntAdj(mesh, vmin, iBase_EDGE, ARRAY_INOUT(edges), &err);

    SimpleArray<iBase_EntityHandle> edgenodes;
    for( size_t i = 0; i < edges.size(); i++) {
        iMesh_getEntAdj(mesh, edges[i], iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);
        if( edgenodes[0] == vmin && edgenodes[1] == vmax )  {
            edgehandle = edges[i];
            return 1;
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool
MeshRefine2D:: allow_edge_refinement( iBase_EntityHandle &edgehandle ) const
{
    int err, bound_val;
    iMesh_getIntData( mesh, edgehandle, boundary_tag, &bound_val, &err);

    // If err == 0, means edge is boundary
    if( !err && !boundary_split_flag ) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

iBase_EntityHandle
MeshRefine2D::create_new_edge( iBase_EntityHandle v1, iBase_EntityHandle v2)
{
    int status, err;
    iBase_EntityHandle edgehandle;
    static vector<iBase_EntityHandle> pnodes(2);

    pnodes[0] = std::min(v1,v2);
    pnodes[1] = std::max(v1,v2);
    iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &pnodes[0], 2, &edgehandle, &status, &err);

    return edgehandle;
}

///////////////////////////////////////////////////////////////////////////////

Point3D
MeshRefine2D:: edge_centroid( iBase_EntityHandle v1, iBase_EntityHandle v2) const
{
    int  err;
    double x, y, z;
    Point3D p3d;

    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;

    iMesh_getVtxCoord(mesh, v1, &x, &y, &z, &err);
    p3d[0] += x;
    p3d[1] += y;
    p3d[2] += z;

    iMesh_getVtxCoord(mesh, v2, &x, &y, &z, &err);
    p3d[0] += x;
    p3d[1] += y;
    p3d[2] += z;

    p3d[0] *= 0.5;
    p3d[1] *= 0.5;
    p3d[2] *= 0.5;

    return p3d;
}

///////////////////////////////////////////////////////////////////////////////

double
MeshRefine2D:: edge_length( iBase_EntityHandle v1, iBase_EntityHandle v2) const
{
    int  err;

    double x0, y0, z0;
    iMesh_getVtxCoord(mesh, v1, &x0, &y0, &z0, &err);

    double x1, y1, z1;
    iMesh_getVtxCoord(mesh, v2, &x1, &y1, &z1, &err);

    double dx = x1-x0;
    double dy = y1-y0;
    double dz = z1-z0;

    double len = sqrt(dx*dx + dy*dy + dz*dz );

    return len;
}

///////////////////////////////////////////////////////////////////////////////

Point3D
MeshRefine2D:: face_centroid( iBase_EntityHandle facehandle) const
{
    int  err;
    double  x, y, z;

    SimpleArray<iBase_EntityHandle> facenodes;
    iMesh_getEntAdj(mesh, facehandle, iBase_VERTEX, ARRAY_INOUT(facenodes), &err);

    size_t numNodes = facenodes.size();
    assert( numNodes > 2 );

    Point3D p3d;
    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;
    for(size_t i =  0; i < numNodes; i++) {
        iMesh_getVtxCoord(mesh, facenodes[i], &x, &y, &z, &err);
        p3d[0] += x;
        p3d[1] += y;
        p3d[2] += z;
    }

    p3d[0] /= (double)numNodes;
    p3d[1] /= (double)numNodes;
    p3d[2] /= (double)numNodes;

    return p3d;
}
///////////////////////////////////////////////////////////////////////////////

double
MeshRefine2D:: face_aspect_ratio( iBase_EntityHandle facehandle) const
{
    int  err;
    double  x, y, z;

    SimpleArray<iBase_EntityHandle> facenodes;
    iMesh_getEntAdj(mesh, facehandle, iBase_VERTEX, ARRAY_INOUT(facenodes), &err);

    size_t numNodes = facenodes.size();
    assert( numNodes > 2 );

    vector<double> elen;
    elen.resize( numNodes );

    Point3D p3d;
    for(size_t i =  0; i < numNodes; i++) {
        iBase_EntityHandle  v0 = facenodes[i];
        iBase_EntityHandle  v1 = facenodes[(i+1)%numNodes];
        elen[i] = edge_length(v0,v1);
    }

    double minlen = *boost::min_element( elen );
    double maxlen = *boost::max_element( elen );

    return minlen/maxlen;
}
///////////////////////////////////////////////////////////////////////////////

int
MeshRefine2D::setVertexOnEdge( iBase_EntityHandle v1, iBase_EntityHandle v2,
                               iBase_EntityHandle &vmid )
{
    int err;
    iBase_EntityHandle edgehandle;

    // Case I: Edge doesn't exist then create and set the tag value
    if( searchEdge(v1,v2,edgehandle) == 0)
    {
        edgehandle = create_new_edge(v1,v2);
        if( allow_edge_refinement(edgehandle) )
        {
            Point3D p3d = edge_centroid(v1,v2);
            iMesh_createVtx(mesh, p3d[0], p3d[1], p3d[2], &vmid, &err);
            iMesh_setEHData(mesh, edgehandle, vertex_on_edge_tag, vmid, &err);
            return 0;
        }
        return 1;
    }

    //Edge exist but not the tag value.
    if( allow_edge_refinement( edgehandle) )
    {
        iMesh_getEHData(mesh, edgehandle, vertex_on_edge_tag, &vmid, &err);
        if( err ) {
            Point3D p3d = edge_centroid(v1,v2);
            iMesh_createVtx(mesh, p3d[0], p3d[1], p3d[2], &vmid, &err);
            iMesh_setEHData(mesh, edgehandle, vertex_on_edge_tag, vmid, &err);
        }
        return 0;
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int MeshRefine2D::getVertexOnEdge( iBase_EntityHandle v1, iBase_EntityHandle v2,
                                   iBase_EntityHandle &hangVertex ) const
{
    int err;
    iBase_EntityHandle edgehandle;

    searchEdge( v1, v2, edgehandle) ;
    iMesh_getEHData(mesh, edgehandle, vertex_on_edge_tag, &hangVertex, &err);

    if( err ) hangVertex = 0;

    return err;
}

///////////////////////////////////////////////////////////////////////////////

int Sqrt3Refine2D :: execute()
{
    /*
       CentroidRefine2D refine(mesh);
       EdgeFlip eflip(mesh);

       for ( int itime = 0; itime < numIterations; itime++) {
           refine.execute();
           eflip.execute();
       }
    */

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int LongestEdgeRefine2D :: initialize()
{
    MeshRefine2D::initialize();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int LongestEdgeRefine2D :: atomicOp( const iBase_EntityHandle &oldface)
{
    int err;
    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);

    iBase_EntityHandle  v1, v2, vmid;
    vector<double> elen(3);

    double maxlen = 0.0;
    for( int i = 0; i < 3; i++) {
        v1      = eConnect[(i+1)%3];
        v2      = eConnect[(i+2)%3];
        elen[i] = edge_length( v1, v2);
        maxlen  = std::max(elen[i], maxlen);
    }

    for( int i = 0; i < 3; i++) {
        if( elen[i]/maxlen > 0.90 )
        {
            v1 = eConnect[(i+1)%3];
            v2 = eConnect[(i+2)%3];
            setVertexOnEdge(v1,v2, vmid);
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int LongestEdgeRefine2D::execute()
{
    int err;
    initialize();

    SimpleArray<iBase_EntityHandle> faces;
    iMesh_getEntities( mesh, rootSet, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                       ARRAY_INOUT(faces), &err);

    size_t ncount = 0;
    for( size_t i = 0; i < faces.size(); i++) {
        double ratio  = face_aspect_ratio( faces[i] );
        if( ratio < cutOffAspectRatio) {
            err = atomicOp( faces[i] );
            if( !err ) ncount++;
        }
    }

    if( ncount ) {
        ConsistencyRefine2D consistency;
        consistency.setMesh(mesh);
        consistency.execute();
        finalize();
    } else
        cout << "Warning: No Edge was refined " << endl;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////


int ConsistencyRefine2D :: initialize()
{
    hangingVertex.resize(3);
    edge0.set(0);
    edge1.set(1);
    edge2.set(2);

    insertedNodes.clear();
    insertedFaces.clear();

    MeshRefine2D::initialize();

    return 0;
}

//#############################################################################

int  ConsistencyRefine2D :: execute()
{
    initialize();

    makeConsistent();

    finalize();

    return 0;
}

//#############################################################################

void ConsistencyRefine2D :: subDivideQuad2Tri( const vector<iBase_EntityHandle> &connect)
{
    int err, status;
    assert( connect.size() == 4 );
    //********************************************************************
    // Subdivide a quadrilateral cell into two triangles. We can choose
    // either quadrilateral diagonal for spliiting, but we choose to
    // select quadrilateral which gives, maximum of minimum aspect
    // ratio of the resulting two triangles ....
    //*******************************************************************
    double diagonal[2];

    iBase_EntityHandle v0 = connect[0];
    iBase_EntityHandle v1 = connect[1];
    iBase_EntityHandle v2 = connect[2];
    iBase_EntityHandle v3 = connect[3];

    diagonal[0]  = edge_length( v0, v3 );
    diagonal[1]  = edge_length( v1, v2 );

    iBase_EntityHandle facehandle;
    vector<iBase_EntityHandle> pnodes(3);

    if( diagonal[0] < diagonal[1] ) {
        pnodes[0] = v0;
        pnodes[1] = v1;
        pnodes[2] = v3;
        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back( facehandle );

        pnodes[0] = v0;
        pnodes[1] = v3;
        pnodes[2] = v2;
        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back( facehandle );
    } else {
        pnodes[0] = v0;
        pnodes[1] = v1;
        pnodes[2] = v2;
        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back( facehandle );

        pnodes[0] = v1;
        pnodes[1] = v3;
        pnodes[2] = v2;

        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back( facehandle );
    }
}

//#############################################################################

void ConsistencyRefine2D :: makeConsistent1( const iBase_EntityHandle &oldface)
{
    int err, status;
    //------------------------------------------------------------------
    // When only one edge is inconsistent, we can direcly join the
    // hanging node to the opposite node of the triangle. Therefore
    // one additional triangle is generated with this case.
    //------------------------------------------------------------------
    SimpleArray<iBase_EntityHandle> pnodes;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(pnodes), &err);

    if( pnodes.size() != 3 ) return;

    iBase_EntityHandle n1 = pnodes[0];
    iBase_EntityHandle n2 = pnodes[1];
    iBase_EntityHandle n3 = pnodes[2];

    int val = 1;
    iMesh_setIntData(mesh, oldface, remove_tag, val, &err);

    iBase_EntityHandle facehandle;

    // If the only hanging vertex lies of the edge0. i.e. opposite to the
    // vertex 0 of the existing triangle.
    if( bitvec == edge0) {
        pnodes[0] = n1;
        pnodes[1] = n2;
        pnodes[2] = hangingVertex[0];

        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back(facehandle);

        // Create a new triangle ....
        pnodes[0] = n3;
        pnodes[1] = n1;
        pnodes[2] = hangingVertex[0];
        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back(facehandle);
        return;
    }

    // If the only hanging vertex lies of the edge1. i.e. opposite to the
    // vertex 1 of the existing triangle.
    if( bitvec == edge1)
    {
        // Replace the existing triangle ....
        pnodes[0] = n2;
        pnodes[1] = n3;
        pnodes[2] = hangingVertex[1];
        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back(facehandle);

        // Create a new triangle ....
        pnodes[0] = n2;
        pnodes[1] = hangingVertex[1];
        pnodes[2] = n1;

        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back(facehandle);

        return;
    }

    // If the only hanging vertex lies of the edge2. i.e. opposite to the
    // vertex 2 of the existing triangle.
    if( bitvec == edge2) {
        // Replace the existing triangle ....
        pnodes[0] = n3;
        pnodes[1] = n1;
        pnodes[2] = hangingVertex[2];

        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back(facehandle);

        // Create a new triangle ....
        pnodes[0] = n3;
        pnodes[1] = hangingVertex[2];
        pnodes[2] = n2;

        iMesh_createEnt(mesh, iMesh_TRIANGLE, &pnodes[0], 3, &facehandle, &status, &err);
        insertedFaces.push_back(facehandle);

        return;
    }

}

//#############################################################################

void ConsistencyRefine2D :: refineEdge0(const iBase_EntityHandle &oldface)
{
    int err, status;
    iBase_EntityHandle facehandle;

    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);
    assert( eConnect.size() == 3 );

    iBase_EntityHandle n1 = eConnect[0];
    iBase_EntityHandle n2 = eConnect[1];
    iBase_EntityHandle n3 = eConnect[2];

    static vector<iBase_EntityHandle> tnodes(3);
    tnodes[0] = n1;
    tnodes[1] = hangingVertex[2];
    tnodes[2] = hangingVertex[1];
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tnodes[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back(facehandle);

    // One Quadrilateral is created, divide into 2 triangles.
    static vector<iBase_EntityHandle> qnodes(4);
    qnodes[0] =  hangingVertex[2];
    qnodes[1] =  n2;
    qnodes[2] =  hangingVertex[1];
    qnodes[3] =  n3;

    subDivideQuad2Tri(qnodes );
}

//#############################################################################

void ConsistencyRefine2D :: refineEdge1(const iBase_EntityHandle &oldface)
{
    int err, status;
    iBase_EntityHandle facehandle;

    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);
    assert( eConnect.size() == 3 );

    iBase_EntityHandle n1 = eConnect[0];
    iBase_EntityHandle n2 = eConnect[1];
    iBase_EntityHandle n3 = eConnect[2];

    static vector<iBase_EntityHandle> tnodes(3);
    tnodes[0] = n2;
    tnodes[1] = hangingVertex[0];
    tnodes[2] = hangingVertex[2];

    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tnodes[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    // One Quadrilateral is created, divide into 2 triangles.
    static vector<iBase_EntityHandle> qnodes(4);
    qnodes[0] =  n1;
    qnodes[1] =  hangingVertex[2];
    qnodes[2] =  n3;
    qnodes[3] =  hangingVertex[0];

    subDivideQuad2Tri( qnodes );
}


//#############################################################################

void ConsistencyRefine2D :: refineEdge2(const iBase_EntityHandle &oldface)
{
    int err, status;
    iBase_EntityHandle facehandle;

    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);
    assert( eConnect.size() == 3 );

    iBase_EntityHandle n1 = eConnect[0];
    iBase_EntityHandle n2 = eConnect[1];
    iBase_EntityHandle n3 = eConnect[2];

    static vector<iBase_EntityHandle> tnodes(3);
    tnodes[0] = n3;
    tnodes[1] = hangingVertex[1];
    tnodes[2] = hangingVertex[0];
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tnodes[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    static vector<iBase_EntityHandle> qnodes(4);
    qnodes[0] =  n1;
    qnodes[1] =  n2;
    qnodes[2] =  hangingVertex[1];
    qnodes[3] =  hangingVertex[0];

    subDivideQuad2Tri( qnodes );
}


//#############################################################################

void ConsistencyRefine2D :: makeConsistent2( const iBase_EntityHandle &oldface)
{
    //--------------------------------------------------------------------
    // When there are two edges which are inconsistent, then we create
    // one triangle and one quadrilateral. This quadrilateral is further
    // divided into 2 triangle, which produces better aspect ratio.
    // Therefore, three triangles are generated in this procedure.
    //--------------------------------------------------------------------
    // Find out which edge is consistent ...
    bitvec.flip();

    int err, val = 1;
    iMesh_setIntData(mesh, oldface, remove_tag, val, &err);

    if( bitvec == edge0) {
        refineEdge0(oldface);
        return;
    }

    if( bitvec == edge1) {
        refineEdge1(oldface);
        return;
    }

    if( bitvec == edge2) {
        refineEdge2(oldface);
        return;
    }

}

//#############################################################################

void ConsistencyRefine2D :: makeConsistent3( const iBase_EntityHandle &oldface)
{
    int err, status;

    int val = 1;
    iMesh_setIntData(mesh, oldface, remove_tag, val, &err);

    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);
    assert( eConnect.size() == 3 );

    iBase_EntityHandle n1 = eConnect[0];
    iBase_EntityHandle n2 = eConnect[1];
    iBase_EntityHandle n3 = eConnect[2];

    iBase_EntityHandle facehandle;

    static vector<iBase_EntityHandle> tnodes(3);
    // First Triangle  Using 0(old), 1,2(new): Replace existing triangle
    tnodes[0] = n1;
    tnodes[1] = hangingVertex[2];
    tnodes[2] = hangingVertex[1];
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tnodes[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back(facehandle);

    // Second Triangle  Using 1(old), 2,0(new): Create new triangle
    tnodes[0] = n2;
    tnodes[1] = hangingVertex[0];
    tnodes[2] = hangingVertex[2];
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tnodes[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back(facehandle);

    // Second Triangle  Using 2(old), 1,0(new): Create new triangle
    tnodes[0] = n3;
    tnodes[1] = hangingVertex[1];
    tnodes[2] = hangingVertex[0];
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tnodes[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back(facehandle);

    //  All new only : Create new triangle
    tnodes[0] = hangingVertex[0];
    tnodes[1] = hangingVertex[1];
    tnodes[2] = hangingVertex[2];
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tnodes[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back(facehandle);
}

//#############################################################################

void ConsistencyRefine2D :: checkFaceConsistency( const iBase_EntityHandle &oldface )
{
    int err;
    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);

    iBase_EntityHandle v1,v2;

    bitvec.reset();
    for( int i = 0; i < 3; i++) {
        v1  = eConnect[ (i+1)%3 ];
        v2  = eConnect[ (i+2)%3 ];
        getVertexOnEdge( v1, v2, hangingVertex[i] );
        if( hangingVertex[i] ) {
            bitvec.set(i);
        }
    }

}

//#############################################################################

int ConsistencyRefine2D :: atomicOp(const iBase_EntityHandle &oldface)
{
    int err, val;
    iMesh_getIntData(mesh, oldface, remove_tag, &val, &err);

    if( err == 0 && val == 1 ) return 1;

    checkFaceConsistency( oldface );

    switch( bitvec.count() )
    {
    case 1:
        numfacesRefined++;
        makeConsistent1( oldface );
        break;
    case 2:
        numfacesRefined++;
        makeConsistent2( oldface );
        break;
    case 3:
        numfacesRefined++;
        makeConsistent3( oldface );
        break;
    }
    return 0;
}

//#############################################################################

void ConsistencyRefine2D :: makeConsistent()
{
    //**********************************************************************
    // The previous step, will leave some hanging nodes. So adjacent cells
    // will be forced to refined. All the hanging nodes are attached
    // directly to the opposite node of the triangle, thus creating two
    // triangle. This may produce some bad triangle, which could be
    // improved using edge swapping algorithm. In this process, only new
    // edges are created and  no new vertices are introduced.
    //**********************************************************************

    int err;
    SimpleArray<iBase_EntityHandle> faces;
    iMesh_getEntities( mesh, rootSet, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                       ARRAY_INOUT(faces), &err);

    for( size_t i = 0; i < faces.size(); i++)
        atomicOp( faces[i] );
}


//#############################################################################

int CentroidRefine2D::refine_tri(const iBase_EntityHandle &oldface)
{
    int err, status;
    iBase_EntityHandle facehandle, vcenter;

    Point3D p3d = face_centroid(oldface);
    iMesh_createVtx(mesh, p3d[0], p3d[1], p3d[2], &vcenter, &err);

    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);

    assert( eConnect.size() == 3 );

    iBase_EntityHandle v1 = eConnect[0];
    iBase_EntityHandle v2 = eConnect[1];
    iBase_EntityHandle v3 = eConnect[2];

    static vector<iBase_EntityHandle> tconn(3);

    tconn[0] = vcenter;
    tconn[1] = v1;
    tconn[2] = v2;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    tconn[0] = vcenter;
    tconn[1] = v2;
    tconn[2] = v3;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    tconn[0] = vcenter;
    tconn[1] = v3;
    tconn[2] = v1;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    int val = 1;
    iMesh_setIntData(mesh, oldface, remove_tag, val, &err);

    return 0;
}


///////////////////////////////////////////////////////////////////////////////

int  CentroidRefine2D::refine_quad(const iBase_EntityHandle &oldface)
{
    int status, err;
    iBase_EntityHandle  facehandle, vcenter;

    Point3D p3d = face_centroid(oldface);
    iMesh_createVtx(mesh, p3d[0], p3d[1], p3d[2], &vcenter, &err);

    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);

    iBase_EntityHandle v0 = eConnect[0];
    iBase_EntityHandle v1 = eConnect[1];
    iBase_EntityHandle v2 = eConnect[2];
    iBase_EntityHandle v3 = eConnect[3];

    vector<iBase_EntityHandle> tconn(3);

    tconn[0] = vcenter;
    tconn[1] = v0;
    tconn[2] = v1;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    tconn[0] = vcenter;
    tconn[1] = v1;
    tconn[2] = v3;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    tconn[0] = vcenter;
    tconn[1] = v3;
    tconn[2] = v2;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    tconn[0] = vcenter;
    tconn[1] = v2;
    tconn[2] = v0;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    int val = 1;
    iMesh_setIntData(mesh, oldface, remove_tag, val, &err);

    return 0;
}


///////////////////////////////////////////////////////////////////////////

int CentroidRefine2D::atomicOp(const iBase_EntityHandle &oldface)
{
    int err, face_topo;
    iMesh_getEntTopo(mesh, oldface, &face_topo, &err);

    switch( face_topo )
    {
    case iMesh_TRIANGLE:
        refine_tri(oldface);
        break;
    case iMesh_QUADRILATERAL:
        refine_quad(oldface);
        break;
    default:
        cout << "Warning: Element not supported for refinement " << endl;
        break;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////

int CentroidRefine2D::initialize()
{
    int err, result;

    const char *tag = "RemoveFace";
    int namelen = strlen( tag );
    iMesh_createTag(mesh, tag, 1, iBase_INTEGER,  &remove_tag, &result, namelen);

    iMesh_getTagHandle(mesh, "GLOBAL_ID", &globalID_tag, &err, strlen("GLOBAL_ID") );

    iMesh_getRootSet(mesh, &rootSet, &err);

    insertedFaces.clear();
}

///////////////////////////////////////////////////////////////////////////


int CentroidRefine2D::execute()
{
    int err;

    initialize();

    SimpleArray<iBase_EntityHandle> faces;
    iMesh_getEntities( mesh, rootSet, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                       ARRAY_INOUT(faces), &err);

    for( int i = 0; i < faces.size(); i++)
        atomicOp( faces[i]  );

    finalize();

    return 0;
}
///////////////////////////////////////////////////////////////////////////

int ObtuseRefine2D :: initialize()
{
    return 0;
};

///////////////////////////////////////////////////////////////////////////

int ObtuseRefine2D :: atomicOp( const iBase_EntityHandle &facehandle)
{
    int err;
    SimpleArray<iBase_EntityHandle> facenodes;
    iMesh_getEntAdj(mesh, facehandle, iBase_VERTEX, ARRAY_INOUT(facenodes), &err);

    Point3D pv0, pv1, pv2;
    Point3D vec1, vec2;

    iBase_EntityHandle vmid;
    double angle;
    for( int i = 0; i < 3; i++) {
        getVertexCoords( facenodes[(i+0)%3], pv0);
        getVertexCoords( facenodes[(i+1)%3], pv1);
        getVertexCoords( facenodes[(i+2)%3], pv2);
        vec1 = JMath::create_vector( pv2, pv0);
        vec2 = JMath::create_vector( pv1, pv0);
        angle = JMath::getVectorAngle(vec1, vec2);
        if( angle > cutoffAngle) {
            setVertexOnEdge(facenodes[(i+1)%2], facenodes[(i+2)%2], vmid );
            return 0;
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////
int ObtuseRefine2D :: execute()
{
    int err;

    initialize();

    SimpleArray<iBase_EntityHandle> faces;
    iMesh_getEntities( mesh, rootSet, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                       ARRAY_INOUT(faces), &err);

    size_t ncount = 0;
    for( size_t i = 0; i < faces.size(); i++) {
        err = atomicOp( faces[i]  );
        if( !err ) ncount++;
    }

    ConsistencyRefine2D refine;
    if( ncount ) {
        refine.setMesh(mesh);
        refine.execute();
        finalize();
    } else
        cout << "Warning: No triangle was refined " << endl;

    return 0;
}

//////////////////////////////////////////////////////////////////////////////

int Refine2D14 :: initialize()
{
    MeshRefine2D::initialize();
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int Refine2D14::refine_quad(const iBase_EntityHandle &oldface)
{
    int status, err;
    iBase_EntityHandle  facehandle, vcenter;

    Point3D p3d = face_centroid(oldface);
    iMesh_createVtx(mesh, p3d[0], p3d[1], p3d[2], &vcenter, &err);

    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);
    assert( eConnect.size() == 4 );

    iBase_EntityHandle v0 = eConnect[0];
    iBase_EntityHandle v1 = eConnect[1];
    iBase_EntityHandle v2 = eConnect[2];
    iBase_EntityHandle v3 = eConnect[3];

    iBase_EntityHandle v01, v23, v02, v13;
    err = setVertexOnEdge(v0,v1, v01);
    err = setVertexOnEdge(v2,v3, v23);
    err = setVertexOnEdge(v0,v2, v02);
    err = setVertexOnEdge(v1,v3, v13);

    static vector<iBase_EntityHandle> qconn(4);

    qconn[0] = v0;
    qconn[1] = v01;
    qconn[2] = v02;
    qconn[3] = vcenter;
    iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &qconn[0], 4, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    qconn[0] = v01;
    qconn[1] = v1;
    qconn[2] = vcenter;
    qconn[3] = v13;
    iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &qconn[0], 4, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    qconn[0] = v02;
    qconn[1] = vcenter;
    qconn[2] = v2;
    qconn[3] = v23;
    iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &qconn[0], 4, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    qconn[0] = vcenter;
    qconn[1] = v13;
    qconn[2] = v23;
    qconn[3] = v3;
    iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &qconn[0], 4, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    int val = 1;
    iMesh_setIntData(mesh, oldface, remove_tag, val, &err);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Refine2D14:: refine_tri( const iBase_EntityHandle &oldface)
{
    int err, status;
    iBase_EntityHandle  facehandle;

    SimpleArray<iBase_EntityHandle> eConnect;
    iMesh_getEntAdj(mesh, oldface, iBase_VERTEX, ARRAY_INOUT(eConnect), &err);
    assert( eConnect.size() == 3 );

    iBase_EntityHandle v0 = eConnect[0];
    iBase_EntityHandle v1 = eConnect[1];
    iBase_EntityHandle v2 = eConnect[2];

    iBase_EntityHandle v01, v12, v20;

    err = setVertexOnEdge(v0,v1, v01);
    err = setVertexOnEdge(v1,v2, v12);
    err = setVertexOnEdge(v2,v0, v20);

    vector<iBase_EntityHandle> tconn(3);

    tconn[0] = v0;
    tconn[1] = v01;
    tconn[2] = v20;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    tconn[0] = v01;
    tconn[1] = v1;
    tconn[2] = v12;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    tconn[0] = v12;
    tconn[1] = v2;
    tconn[2] = v20;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    tconn[0] = v01;
    tconn[1] = v12;
    tconn[2] = v20;
    iMesh_createEnt(mesh, iMesh_TRIANGLE, &tconn[0], 3, &facehandle, &status, &err);
    insertedFaces.push_back( facehandle );

    int val = 1;
    iMesh_setIntData(mesh, oldface, remove_tag, val, &err);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Refine2D14::atomicOp(const iBase_EntityHandle &oldface)
{
    int err, face_topo;
    iMesh_getEntTopo(mesh, oldface, &face_topo, &err);

    switch( face_topo )
    {
    case iMesh_TRIANGLE:
        refine_tri(oldface);
        break;
    case iMesh_QUADRILATERAL:
        refine_quad(oldface);
        break;
    default:
        cout << "Warning: Element not supported for refinement " << endl;
        break;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Refine2D14 ::execute()
{
    int err;
    initialize();

    SimpleArray<iBase_EntityHandle> faces;
    iMesh_getEntities( mesh, rootSet, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                       ARRAY_INOUT(faces), &err);

    for( size_t i = 0; i < faces.size(); i++)
        atomicOp( faces[i]  );

    ConsistencyRefine2D refine;

    refine.setMesh(mesh);
    refine.execute();

    MeshRefine2D::finalize();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int GradeRefine2D :: initialize()
{
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int GradeRefine2D :: atomicOp( const iBase_EntityHandle &apexVertex)
{
    /*
        SimpleArray<iBase_EntityHandle> vneighs;
        iMesh_getEntAdj(mesh, apexVertex, iBase_VERTEX, ARRAY_INOUT(vneighs), &err);

        size_t numNeighs = vneighs.size();

        if( numNeighs == 0 ) return 0;

        vector<double> elen;
        elen.resize( numNeighs);

        for( int i = 0; i < numNeighs; i++)
             elen[i] = length( vertex, vneighs[i] );

        sort( elen.begin(), elen.end() );

        double median_value = elen[numNeighs/2];

        for( int i = 0; i < numNeighs; i++) {
            if( elen[i] > 0.90*median_value)
                setVertexOnEdge( vertex, vneighs[i]);
        }
    */
    return 1;
}
////////////////////////////////////////////////////////////////////////////////

int GradeRefine2D :: finalize()
{
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int GradeRefine2D :: execute()
{
    /*
        SimpleArray<iBase_EntityHandle> nodes;
        iMesh_getEntities( mesh, rootSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                           ARRAY_INOUT(nodes), &err);

        initialize();

        for( int i = 0; i < nodes.size(); i++)
             atomicOp( nodes[i] );

        ConsistencyRefine2D refine(mesh, edgemap);
        refine.execute();
    */

    finalize();
}

///////////////////////////////////////////////////////////////////////////////

#ifdef TEST_MESHREFINE
int main(int argc, char **argv)
{
    iMesh_Instance mesh = read_off_file( "model.off");

//  LongestEdgeRefine2D  meshrefine;
    Refine2D14  meshrefine;
    meshrefine.setBoundarySplitFlag(0);
    meshrefine.setMesh( mesh );
    meshrefine.execute();

    write_off_file( mesh, "refine.off");

    return 0;
}
#endif

