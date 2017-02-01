#include "MeshRefine.hpp"

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////
void JMeshRefiner2D :: makeConsistent()
{
    /*
        if( inConsistentSet.empty() ) return;
        JFaceSequence faceSeq;
        copy( inConsistentSet.begin(), inConsistentSet.end(), back_inserter(faceSeq));

        ConsistencyRefiner2D refiner;
        refiner.setMesh(mesh);
        refiner.refine(faceSeq);
    */
}

////////////////////////////////////////////////////////////////////////////////
void JMeshRefiner2D :: upgradeInconsistent()
{
    if( mesh == nullptr) return;
    size_t numfaces = mesh->getSize(2);

    JNodePtr edgenode;
    JFacePtr poly;
    JNodeSequence polynodes;

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            bool hangingnode = 0;
            for( int j = 0; j < face->getSize(0); j++) {
                const JEdgePtr &edge = face->getEdgeAt(i);
                if( edge->hasAttribute("Steiner"))  hangingnode = 1;
            }
            if( hangingnode ) {
                polynodes.clear();
                for( int j = 0; j < face->getSize(0); j++) {
                    polynodes.push_back( face->getNodeAt(i) );
                    const JEdgePtr &edge = face->getEdgeAt(i);
                    int err = edge->getAttribute("Steiner", edgenode);
                    if( !err) polynodes.push_back(edgenode);
                }
                if( polynodes.size() == 4)
                    poly = JQuadrilateral::newObject();
                else
                    poly = JPolygon::newObject();
                poly->setNodes( polynodes);
                face->setStatus(JMeshEntity::REMOVE);
                outmesh->addObject( poly );
            }
        }
    }

}
////////////////////////////////////////////////////////////////////////////////

int JMeshRefiner2D :: initialize()
{
    newnodes.clear();
    newedges.clear();
    newfaces.clear();
    newcells.clear();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JMeshRefiner2D :: finalize()
{
    if( outmesh ) {
        outmesh->deleteEdgeAttribute("Steiner");
        outmesh->deleteFaceAttribute("Steiner");
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

JNodePtr JMeshRefiner2D ::setNodeOnEdge( const JEdgePtr &edge ) const
{
    JNodePtr v;
    if( !allow_edge_refinement(edge) ) return v;

    if( edge->hasAttribute("Steiner") ) {
        edge->getAttribute("Steiner", v);
        return v;
    }

    const JNodePtr &v0 = edge->getNodeAt(0);
    const JNodePtr &v1 = edge->getNodeAt(1);

    Point3D p3d = JNodeGeometry::getMidPoint(v0, v1);
    v = JNode::newObject();
    v->setXYZCoords( p3d );
    edge->setAttribute("Steiner", v);
    outmesh->addObject(v);
    return v;
}

///////////////////////////////////////////////////////////////////////////////

JNodePtr JMeshRefiner2D ::getNodeOnEdge( const JEdgePtr &edge ) const
{
    JNodePtr v;
    if( edge->hasAttribute("Steiner") ) {
        edge->getAttribute("Steiner", v);
    }
    return v;
}

///////////////////////////////////////////////////////////////////////////////

bool
JMeshRefiner2D ::allow_edge_refinement( const JEdgePtr &edge) const
{
    if( edge->hasAttribute("Boundary") && boundary_split_flag == 0) return 0;
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRefiner2D::append( const JNodePtr &vertex)
{
    if(outmesh) outmesh->addObject( vertex );
    newnodes.push_back( vertex );
}

///////////////////////////////////////////////////////////////////////////////

JFacePtr JMeshRefiner2D::append( const JNodePtr &v0, const JNodePtr &v1 , const JNodePtr &v2)
{
    JFacePtr face = JTriangle::newObject( v0, v1, v2);
    if(outmesh) outmesh->addObject(face);
    newfaces.push_back( face );
    return face;
}

///////////////////////////////////////////////////////////////////////////////

JFacePtr JMeshRefiner2D::append( const JNodePtr &v0, const JNodePtr &v1 , const JNodePtr &v2, const JNodePtr &v3)
{
    JFacePtr face = JQuadrilateral::newObject(v0,v1,v2,v3);
    if(outmesh) outmesh->addObject(face);
    newfaces.push_back( face );
    return face;
}

///////////////////////////////////////////////////////////////////////////////

int JLongestEdgeRefiner2D :: atomicOp( const JFacePtr &oldface)
{
    vector<double> elen(3);

    double maxlen = 0.0;
    for( int i = 0; i < 3; i++) {
        const JEdgePtr &edge = oldface->getEdgeAt(i);
        elen[i] = JEdgeGeometry::getLength(edge);
        maxlen  = std::max(elen[i], maxlen);
    }

    for( int i = 0; i < 3; i++) {
        if( elen[i]/maxlen > 0.90 ) {
            const JEdgePtr &edge = oldface->getEdgeAt(i);
            setNodeOnEdge( edge );
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JLongestEdgeRefiner2D::refineAll()
{
    if( mesh == nullptr ) return 1;

    size_t numfaces = mesh->getSize(2);

    size_t ncount = 0;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        double ratio  = JFaceGeometry::getAspectRatio(face);
        if( ratio < cutOffAspectRatio) {
            int err = atomicOp( face );
            if( !err ) ncount++;
        }
    }

    /*
        if( ncount ) {
            ConsistencyRefiner2D refiner;
            refiner.setMesh(mesh);
            refiner.refineAll();
        }
    */

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void JConsistencyRefiner2D :: refineAll()
{
    if( mesh == nullptr ) return;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces ; i++)
        atomicOp( mesh->getFaceAt(i) );
}

///////////////////////////////////////////////////////////////////////////////

void JConsistencyRefiner2D :: refine( JFaceSequence &faces2refine)
{
    size_t numfaces = faces2refine.size();
    for( size_t i = 0; i < numfaces ; i++)
        atomicOp( faces2refine[i] );
}

//#############################################################################

void JConsistencyRefiner2D :: quad2tris( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3)
{
    //********************************************************************
    // Subdivide a quadrilateral cell into two triangles. We can choose
    // either quadrilateral diagonal for spliiting, but we choose to
    // select quadrilateral which gives, maximum of minimum aspect
    // ratio of the resulting two triangles ....
    //*******************************************************************
    double diagonal[2];

    diagonal[0]  = JNodeGeometry::getLength( v0, v3 );
    diagonal[1]  = JNodeGeometry::getLength( v1, v2 );

    if( diagonal[0] < diagonal[1] ) {
        append( v0, v1, v3);
        append( v0, v3, v2);
    } else {
        append( v0, v1, v2);
        append( v1, v3, v2);
    }
}

//#############################################################################

void JConsistencyRefiner2D :: makeConsistent1( const JFacePtr &oldface)
{
    //------------------------------------------------------------------
    // When only one edge is inconsistent, we can direcly join the
    // hanging node to the opposite node of the triangle. Therefore
    // one additional triangle is generated with this case.
    //------------------------------------------------------------------
    if( oldface->getSize(0) != 3 ) return;

    const JNodePtr &n0 = oldface->getNodeAt(0);
    const JNodePtr &n1 = oldface->getNodeAt(1);
    const JNodePtr &n2 = oldface->getNodeAt(2);

    oldface->setStatus( JMeshEntity::REMOVE );

    if( edgeNode[0] ) {
        append(n1, n2, edgeNode[0] );
        append(n2, n0, edgeNode[0] );
        return;
    }


    if( edgeNode[1] ) {
        append(n0, n1, edgeNode[1] );
        append(n2, n0, edgeNode[1] );
        return;
    }

    if( edgeNode[2] ) {
        append(n0, n1, edgeNode[2] );
        append(n1, n2, edgeNode[2] );
        return;
    }

}

//#############################################################################

void JConsistencyRefiner2D :: makeConsistent2( const JFacePtr &oldface)
{
    //--------------------------------------------------------------------
    // When there are two edges which are inconsistent, then we create
    // one triangle and one quadrilateral. This quadrilateral is further
    // divided into 2 triangle, which produces better aspect ratio.
    // Therefore, three triangles are generated in this procedure.
    //--------------------------------------------------------------------
    const JNodePtr &n0 = oldface->getNodeAt(0);
    const JNodePtr &n1 = oldface->getNodeAt(1);
    const JNodePtr &n2 = oldface->getNodeAt(2);

    if( edgeNode[0] != nullptr && edgeNode[1] != nullptr && edgeNode[2] == nullptr ) {
        quad2tris(n0, edgeNode[0], edgeNode[1], n2);
        append(n1, edgeNode[1], edgeNode[0]);
    }

    if( edgeNode[0] == nullptr && edgeNode[1] != nullptr && edgeNode[2] != nullptr ) {
        quad2tris(n0, n1, edgeNode[1], edgeNode[2]);
        append(n2, edgeNode[2], edgeNode[1]);
    }

    if( edgeNode[0] != nullptr && edgeNode[1] == nullptr && edgeNode[2] != nullptr ) {
        quad2tris(n1, n2, edgeNode[2], edgeNode[0]);
        append(n0, edgeNode[0], edgeNode[2]);
    }

    oldface->setStatus( JMeshEntity::REMOVE );
}

//#############################################################################

void JConsistencyRefiner2D :: makeConsistent3( const JFacePtr &oldface)
{
    const JNodePtr &n0 = oldface->getNodeAt(0);
    const JNodePtr &n1 = oldface->getNodeAt(1);
    const JNodePtr &n2 = oldface->getNodeAt(2);

    append( n0, edgeNode[0], edgeNode[2] );
    append( n1, edgeNode[1], edgeNode[0] );
    append( n2, edgeNode[2], edgeNode[1] );
    append( edgeNode[0], edgeNode[1], edgeNode[2] );

    oldface->setStatus( JMeshEntity::REMOVE );
}

//#############################################################################

int JConsistencyRefiner2D :: atomicOp( const JFacePtr &oldface)
{
    if( !oldface->isActive() )  return 1;

    edgeNode[0] = nullptr;
    edgeNode[1] = nullptr;
    edgeNode[2] = nullptr;

    int nCount = 0;
    for( int i = 0; i < 3; i++) {
        const JEdgePtr &edge = oldface->getEdgeAt(i);
        if( edge->hasAttribute("Steiner") ) {
            edge->getAttribute("Steiner", edgeNode[i] );
            nCount++;
        }
    }

    switch( nCount ) {
    case 1:
        makeConsistent1( oldface );
        break;
    case 2:
        makeConsistent2( oldface );
        break;
    case 3:
        makeConsistent3( oldface );
        break;
    }

    if( nCount) {
        numfacesRefined++;
        oldface->setStatus( JMeshEntity::REMOVE );
    }

    return 0;
}

//#############################################################################

int JCentroidRefiner2D::refine_tri(const JFacePtr &oldface)
{
    JNodePtr vcenter = JNode::newObject();

    Point3D pc;
    oldface->getAvgXYZ( pc );
    vcenter->setXYZCoords( pc );

    const JNodePtr &v1 = oldface->getNodeAt( 0 );
    const JNodePtr &v2 = oldface->getNodeAt( 1 );
    const JNodePtr &v3 = oldface->getNodeAt( 2 );

    append( vcenter );
    append( vcenter, v1, v2);
    append( vcenter, v2, v3);
    append( vcenter, v3, v1);

    oldface->setStatus( JMeshEntity::REMOVE );
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JCentroidRefiner2D::refine_quad(const JFacePtr &oldface)
{
    JNodePtr vcenter = JNode::newObject();

    Point3D pc;
    oldface->getAvgXYZ( pc );
    vcenter->setXYZCoords( pc );

    const JNodePtr &v0 = oldface->getNodeAt( 0 );
    const JNodePtr &v1 = oldface->getNodeAt( 1 );
    const JNodePtr &v2 = oldface->getNodeAt( 2 );
    const JNodePtr &v3 = oldface->getNodeAt( 3 );

    append( vcenter, v0, v1);
    append( vcenter, v1, v3);
    append( vcenter, v3, v2);
    append( vcenter, v2, v0);

    oldface->setStatus( JMeshEntity::REMOVE );

    return 0;
}

///////////////////////////////////////////////////////////////////////////

int JCentroidRefiner2D::atomicOp(const JFacePtr &oldface)
{
    if( !oldface->isActive() ) return 1;

    if( oldface->getSize(0) == 3 ) return refine_tri( oldface );
    if( oldface->getSize(0) == 4 ) return refine_quad( oldface );

    cout << "Warning: Element not supported for refinement " << endl;

    return 1;
}

///////////////////////////////////////////////////////////////////////////

void JCentroidRefiner2D::refineAll(int n )
{
    if( mesh == nullptr ) return;

    newnodes.clear();
    newedges.clear();
    newfaces.clear();
    newcells.clear();

    for( int j = 0; j < numIters; j++) {
        size_t nSize = mesh->getSize(2);
        for( size_t i = 0; i < nSize; i++)
            atomicOp( mesh->getFaceAt(i)  );
    }

}

///////////////////////////////////////////////////////////////////////////

int JObtuseRefiner2D :: atomicOp( const JFacePtr &oldface)
{
    assert( oldface != nullptr );
    /*
    Point3D pv0, pv1, pv2;
    Point3D vec1, vec2;

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
    */

    return 1;
}

///////////////////////////////////////////////////////////////////////////
int JObtuseRefiner2D :: refineAll()
{
    if ( mesh == nullptr ) return 1;

    size_t numfaces = mesh->getSize(2);

    size_t ncount = 0;
    for( size_t i = 0; i < numfaces; i++) {
        int err = atomicOp( mesh->getFaceAt(i) );
        if( !err ) ncount++;
    }

    /*
        if( ncount ) {
            ConsistencyRefiner2D refiner;
            refiner.setMesh(mesh);
            refiner.refineAll();
        }
    */
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

int JRefiner2D14::refine_quad(const JFacePtr &oldface)
{
    const JNodePtr &v0 = oldface->getNodeAt( 0 );
    const JNodePtr &v1 = oldface->getNodeAt( 1 );
    const JNodePtr &v2 = oldface->getNodeAt( 2 );
    const JNodePtr &v3 = oldface->getNodeAt( 3 );

    JNodePtr edgenodes[4];
    for( int i = 0; i < 4; i++) {
        const JEdgePtr &edge  = oldface->getEdgeAt(i);
        edgenodes[i] = setNodeOnEdge( edge );
        assert( edgenodes[i] );
    }

    JNodePtr v01  = edgenodes[0];
    JNodePtr v12  = edgenodes[1];
    JNodePtr v23  = edgenodes[2];
    JNodePtr v30  = edgenodes[3];

    Point3D pc;
    oldface->getAvgXYZ( pc );

    JNodePtr vc = JNode::newObject();
    vc->setXYZCoords( pc );
    outmesh->addObject(vc);

    append( v0, v01, vc, v30);
    append( v1, v12, vc, v01);
    append( v2, v23, vc, v12);
    append( v3, v30, vc, v23 );

    oldface->setStatus( JMeshEntity::REMOVE );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JRefiner2D14:: refine_tri( const JFacePtr &oldface)
{
    JNodePtr v0 = oldface->getNodeAt( 0 );
    JNodePtr v1 = oldface->getNodeAt( 1 );
    JNodePtr v2 = oldface->getNodeAt( 2 );

    JNodePtr edgenodes[3];
    for( int i = 0; i < 3; i++) {
        const JEdgePtr &edge   = oldface->getEdgeAt(i);
        edgenodes[i] = setNodeOnEdge(edge);
    }

    JNodePtr v01 = edgenodes[0];
    JNodePtr v12 = edgenodes[1];
    JNodePtr v20 = edgenodes[2];

    append( v0, v01, v20);
    append( v1, v12, v01);
    append( v2, v20, v12);
    append( v01, v12, v20);

    oldface->setStatus( JMeshEntity::REMOVE );

    return 0;

}

///////////////////////////////////////////////////////////////////////////////

int JRefiner2D14::atomicOp(const JFacePtr &oldface)
{
    if( !oldface->isActive() )  return 1;

    int numnodes = oldface->getSize(0);

    if( numnodes == 3 ) return refine_tri( oldface );
    if( numnodes == 4 ) return refine_quad( oldface );

    if( numnodes > 4)
        cout << "Warning: Element not supported for refinement " << endl;
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int JRefiner2D14 :: refineAll()
{
    if( mesh == nullptr ) return 1;

    initialize();

    size_t numfaces = mesh->getSize(2);

    size_t ncount = 0;
    for( size_t i = 0; i < numfaces; i++) {
        int err = atomicOp( mesh->getFaceAt(i)  );
        if( !err ) ncount++;
    }

    finalize(); // Remove edges/faces temporary attributes ....

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGradeRefine2D :: initialize()
{
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGradeRefine2D :: atomicOp( const JNodePtr &)
{
    /*
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


int JGradeRefine2D :: execute()
{
    initialize();

    size_t numnodes = mesh->getSize(0);

    size_t ncount = 0;
    for( size_t i = 0; i < numnodes; i++) {
        int err = atomicOp( mesh->getNodeAt(i) );
        if( !err ) ncount++;
    }

    /*
        if( ncount ) {
            ConsistencyRefine2D refine(mesh, edgemap);
            refine.execute();
            finalize();
        }
    */

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

