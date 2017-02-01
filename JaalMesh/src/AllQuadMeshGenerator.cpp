#include "AllQuadMeshGenerator.hpp"

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////
void AllQuadMeshGenerator :: setMesh( const JMeshPtr &m)
{
    mesh    = m;
    if( mesh == nullptr) return;

    /*
        triMesh = m;
        size_t numfaces = mesh->getSize(2);

        JNodeSequence conn(3);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                int nn = face->getSize(0);
                if( nn >3) {
                    for( int j = 0; j < nn-2; j++) {
                        conn[0] = face->getNodeAt(0);
                        conn[1] = face->getNodeAt((j+1)%nn);
                        conn[2] = face->getNodeAt((j+2)%nn);
                        JFacePtr tri = Triangle::newObject(conn);
                        triMesh->addObject(tri);
                    }
                    face->setStatus(JMeshEntity::REMOVE);
                }
            }
        }
        triMesh->pruneFaces();
        triMesh->getTopology()->searchBoundary();
    */
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator :: getBaseQuadMesh(int algo)
{
    if( algo == 0) return getBaseQuads_Tri2Quads();
    if( algo == 1) return getBaseQuads_TriMatch();
    if( algo == 2) return getBaseQuads_Skeleton();
    if( algo == 4) return getBaseQuads_EdgeMatch();
    return nullptr;
}
////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator :: getBaseTriMesh(bool resample)
{
    if( mesh == nullptr) return nullptr;

    JDelaunayMesh2D del;

    if( !resample )  {
        del.setMesh(mesh);
        triMesh = del.getSimpleMesh();
        return triMesh;
    }

    JMeshPtr halfmesh = JMesh::newObject();

    vector<JEdgeSequence> boundloops;
    mesh->getTopology()->getBoundary(boundloops);

    JNodeSequence boundnodes;
    int numloops = boundloops.size();
    int indx = 0;
    Point3D pmid;

    skippedNodes.clear();
    midPoints.clear();

    for( int i = 0; i < numloops; i++) {
        JEdgeTopology::getChainNodes(boundloops[i], boundnodes);
        int nnodes = boundnodes.size();
        if( nnodes%2 ) {
            cout << "Error: Number of edges is not even " << endl;
            return nullptr;
        }
        for( int j = 0; j < nnodes/2; j++) {
            halfmesh->addObject(boundnodes[2*j] );
            skippedNodes.push_back( boundnodes[2*j+1]);
            boundnodes[2*j]->setID(indx++);
        }
        for( int j = 0; j < nnodes/2; j++) {
            const JNodePtr &v1 = boundnodes[2*j];
            const JNodePtr &v2 = boundnodes[(2*j+2)%nnodes];
            pmid = JNodeGeometry::getMidPoint(v1,v2);
            midPoints.push_back(pmid);
            JEdgePtr edge = JEdge::newObject(v1,v2);
            halfmesh->addObject( edge );
        }
    }

    del.setMesh(halfmesh);
    triMesh = del.getSimpleMesh();
    return triMesh;
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator :: getBaseQuads_Tri2Quads()
{
    getBaseTriMesh( 1 );
    getSimpleTris2Quads();

    JNodeSequence newBoundNodes;
    quadMesh->getTopology()->getBoundary( newBoundNodes);
    int nmid = midPoints.size();

    for( int i = 0; i < nmid; i++) {
        int found = 0;
        for( const JNodePtr &vtx : newBoundNodes) {
            const Point3D &p1 = vtx->getXYZCoords();
            double dist  = JMath::length(midPoints[i], p1);
            if( dist < 1.0E-10) {
                const Point3D &pold =  skippedNodes[i]->getXYZCoords();
                vtx->setXYZCoords( pold );
                found = 1;
                break;
            }
        }
        assert( found );
    }
    midPoints.clear();
    skippedNodes.clear();

    return quadMesh;
}
////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator :: getBaseQuads_TriMatch()
{
    getBaseTriMesh(0);
    JTriRefiner  refiner;
    refiner.setMesh(triMesh);
    refiner.refineAll(13);
    getEdmondsTrianglesMatching();
    JDoublet doublet;
    doublet.setMesh(quadMesh);
    doublet.removeAll();
    return quadMesh;
}
////////////////////////////////////////////////////////////////////////////////


JMeshPtr AllQuadMeshGenerator :: getBaseQuads_EdgeMatch()
{
    /*
        getBaseTriMesh(0);
        JEdmondGraphMatching matching;
        matching.setMesh(triMesh);
        quadMesh = matching.getEdgeMatching();
        return quadMesh;
    */
    return nullptr;
}
////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator :: getBaseQuads_Skeleton()
{
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator :: getTrianglesMatching(int algo)
{
    quadMesh.reset();
    switch(algo)
    {
    case GREEDY_MATCHING:
        getGreedyTrianglesMatching();
        break;
    case EDMONDS_MATCHING:
        getEdmondsTrianglesMatching();
        break;
    case BINARY_TREE_MATCHING:
        getBinaryTreeMatching();
        break;
    }

    if( quadMesh) {
       triMesh->deleteFaces();
       quadMesh->getTopology()->searchBoundary();
    }

    return quadMesh;
}

////////////////////////////////////////////////////////////////////////////////

int AllQuadMeshGenerator :: getGreedyTrianglesMatching()
{
    if( triMesh == nullptr) return 1;

    size_t numfaces = triMesh->getSize(2);

    JFaceSequence trifaces;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = triMesh->getFaceAt(i);
        if( f->isActive() &&  f->getSize(0) == 3) trifaces.push_back(f);
    }
    if( trifaces.empty() ) return 0;

    quadMesh = JMesh::newObject();

    JNodeSequence nodes = triMesh->getNodes();
    quadMesh->addObjects(nodes);

    // Some quads are already present ...
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = triMesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) == 4)
            quadMesh->addObject(face);
    }

    // Some needs to merge .
    JNodeSequence qnodes(4);
    JFaceSequence faceneighs;

    JEdgeSequence commonedges;
    for( size_t i = 0; i < trifaces.size()/2; i++) {
        const JFacePtr &tri0 = mesh->getFaceAt(2*i);
        const JFacePtr &tri1 = mesh->getFaceAt(2*i+1);
        JFace::getSharedEntities( tri0, tri1, commonedges);
        if( !commonedges.empty() ) {
            if( tri0->getSize(0) == 3 && tri1->getSize(0) == 3) {
                qnodes[1] = commonedges[0]->getNodeAt(0);
                qnodes[3] = commonedges[0]->getNodeAt(1);
                qnodes[0] = JTriangle::getOppositeNode(tri0, qnodes[1], qnodes[3]);
                qnodes[2] = JTriangle::getOppositeNode(tri1, qnodes[1], qnodes[3]);
                JFacePtr qface = JQuadrilateral::newObject( qnodes);
                quadMesh->addObject(qface);
                tri0->setStatus(JMeshEntity::REMOVE);
                tri1->setStatus(JMeshEntity::REMOVE);
                commonedges[0]->setStatus(JMeshEntity::REMOVE);
            }
        }
    }

    for( const JFacePtr &triface : trifaces) {
        if( triface->isActive() ) {
            for( int i = 0; i < 3; i++) {
                const JEdgePtr &edge = triface->getEdgeAt(i);
                JEdge::getRelations(edge, faceneighs);
                if( faceneighs.size() == 2) {
                    const JFacePtr &f0 = faceneighs[0];
                    const JFacePtr &f1 = faceneighs[1];
                    if( f0->getSize(0) == 3 && f1->getSize(0) == 3) {
                        qnodes[1] = edge->getNodeAt(0);
                        qnodes[3] = edge->getNodeAt(1);
                        qnodes[0] = JTriangle::getOppositeNode(f0, qnodes[1], qnodes[3]);
                        qnodes[2] = JTriangle::getOppositeNode(f1, qnodes[1], qnodes[3]);
                        JFacePtr qface = JQuadrilateral::newObject( qnodes);
                        quadMesh->addObject(qface);
                        f0->setStatus(JMeshEntity::REMOVE);
                        f1->setStatus(JMeshEntity::REMOVE);
                        edge->setStatus(JMeshEntity::REMOVE);
                        break;
                    }
                }
            }
        }
    }

// Some left out ...
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = triMesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) == 3)
            quadMesh->addObject(face);
    }

    quadMesh->getTopology()->collectEdges();
    triMesh->clearAll();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
JMeshPtr AllQuadMeshGenerator :: getCatmullClarkMesh()
{
    JMeshPtr quadMesh = JMesh::newObject();
    if( mesh == nullptr) return quadMesh;

    quadMesh->addObjects(mesh->getNodes());
    size_t numedges = mesh->getSize(1);

    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        JNodePtr  newnode = JEdgeGeometry::getMidNode(edge);
        quadMesh->addObject(newnode);
        edge->setAttribute("Steiner", newnode);
    }

    size_t numfaces = mesh->getSize(2);
    JNodeSequence edgenodes;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &oldface = mesh->getFaceAt(i);
        if( oldface->isActive() ) {
            JNodePtr newnode = JFaceGeometry::getCentroid(oldface);
            quadMesh->addObject(newnode);
            int nn = oldface->getSize(0);
            edgenodes.resize(nn);
            for( int j = 0; j < nn; j++) {
                const JEdgePtr &edge = oldface->getEdgeAt(j);
                edge->getAttribute("Steiner", edgenodes[j] );
            }
            for( int j = 0; j < nn; j++) {
                JNodePtr v0 = oldface->getNodeAt(j);
                JNodePtr v1 = edgenodes[j];
                JNodePtr v2 = newnode;
                JNodePtr v3 = edgenodes[(j + nn-1)%nn];
                JFacePtr newquad = JQuadrilateral::newObject(v0,v1,v2,v3);
                quadMesh->addObject(newquad);
            }
            oldface->setStatus(JMeshEntity::REMOVE);
        }
    }

    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        edge->setStatus(JMeshEntity::REMOVE);
    }
    quadMesh->enumerate(0);
    quadMesh->enumerate(0);
    quadMesh->enumerate(2);
    cout << "#Quads " <<quadMesh->getSize(0) << "  " << quadMesh->getSize(2) << endl;
    return quadMesh;
}
////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator :: getStructuredMesh( int *dim,  double *len, double *org, bool texCoord)
{
    JMeshPtr quadMesh;
    int nx = dim[0];
    int ny = dim[1];
    if( nx < 2 || ny < 2 ) return quadMesh;

    quadMesh = JMesh::newObject();

    double xlen = 1.0, ylen = 1.0;
    if( len ) {
        xlen = len[0];
        ylen = len[1];
    }

    double xorg = 0.0, yorg = 0.0;
    if( org ) {
        xorg = org[0];
        yorg = org[1];
    }

    double dx = xlen/(double)(nx - 1);
    double dy = ylen/(double)(ny - 1);

    double du = 1.0/(double)(nx-1);
    double dv = 1.0/(double)(ny-1);

    Point3D xyz;
    Point2D uv;

    int index = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            xyz[0] = xorg + i * dx;
            xyz[1] = yorg + j * dy;
            xyz[2] = 0.0;
            JNodePtr vnew = JNode::newObject();
            vnew->setXYZCoords(xyz);
            if( texCoord) {
                uv[0] = i*du;
                uv[1] = j*dv;
                vnew->setAttribute("UVCoords", uv);
            }
            vnew->setID(index++);
            quadMesh->addObject(vnew);
        }
    }

    JNodeSequence nodes(4);
    index = 0;
    JFacePtr newquad;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            int n0 = j * nx + i;
            int n1 = n0 + 1;
            int n2 = n1 + nx;
            int n3 = n0 + nx;
            nodes[0] = quadMesh->getNodeAt(n0);
            nodes[1] = quadMesh->getNodeAt(n1);
            nodes[2] = quadMesh->getNodeAt(n2);
            nodes[3] = quadMesh->getNodeAt(n3);
            newquad =  JQuadrilateral::newObject(nodes);
            newquad->setID(index++);
            quadMesh->addObject(newquad);
        }
    }

    int bound = 1;
    for( int i = 0; i < nx-1; i++) {
        const JNodePtr &v0 =  quadMesh->getNodeAt(i);
        const JNodePtr &v1 =  quadMesh->getNodeAt(i+1);
        JEdgePtr edge = JSimplex::getEdgeOf(v0,v1);
        edge->setAttribute("Boundary", bound);
    }

    bound = 2;
    for( int j = 0; j < ny-1; j++) {
        int offset  =  j*nx + nx-1;
        const JNodePtr &v0 =  quadMesh->getNodeAt(offset);
        const JNodePtr &v1 =  quadMesh->getNodeAt(offset+ nx);
        JEdgePtr edge = JSimplex::getEdgeOf(v0,v1);
        edge->setAttribute("Boundary", bound);
    }

    bound = 3;
    for( int i = 0; i < nx-1; i++) {
        int offset  = (ny-1)*nx + i;
        const JNodePtr &v0 =  quadMesh->getNodeAt(offset);
        const JNodePtr &v1 =  quadMesh->getNodeAt(offset+1);
        JEdgePtr edge = JSimplex::getEdgeOf(v0,v1);
        edge->setAttribute("Boundary", bound);
    }

    bound = 4;
    for( int j = 0; j < ny-1; j++) {
        int offset  =  j*nx;
        const JNodePtr &v0 =  quadMesh->getNodeAt(offset);
        const JNodePtr &v1 =  quadMesh->getNodeAt(offset+nx);
        JEdgePtr edge = JSimplex::getEdgeOf(v0,v1);
        edge->setAttribute("Boundary", bound);
    }

    return quadMesh;
}

////////////////////////////////////////////////////////////////////////////////
JMeshPtr AllQuadMeshGenerator :: getSquareHoles( int nx, int ny, double iLength, double oLength)
{
    assert( oLength > iLength );

    double len[2];
    len[0] = oLength*nx;
    len[1] = oLength*ny;
    int dim[2];
    dim[0] = nx+1;
    dim[1] = ny+1;

    JMeshPtr mesh = AllQuadMeshGenerator::getStructuredMesh( dim, len);

    size_t numFaces = mesh->getSize(2);
    Point3D org, xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        const JNodePtr &v00 = face->getNodeAt(0);
        const JNodePtr &v01 = face->getNodeAt(1);
        const JNodePtr &v02 = face->getNodeAt(2);
        const JNodePtr &v03 = face->getNodeAt(3);
        JFaceGeometry::getCentroid( v00, v01, v02, v03, org);
        JNodePtr v10 = JNode::newObject();
        JNodePtr v11 = JNode::newObject();
        JNodePtr v12 = JNode::newObject();
        JNodePtr v13 = JNode::newObject();
        xyz[0] = org[0] - 0.5*iLength;
        xyz[1] = org[1] - 0.5*iLength;
        v10->setXYZCoords(xyz);
        xyz[0] = org[0] + 0.5*iLength;
        xyz[1] = org[1] - 0.5*iLength;
        v11->setXYZCoords(xyz);
        xyz[0] = org[0] + 0.5*iLength;
        xyz[1] = org[1] + 0.5*iLength;
        v12->setXYZCoords(xyz);
        xyz[0] = org[0] - 0.5*iLength;
        xyz[1] = org[1] + 0.5*iLength;
        v13->setXYZCoords(xyz);
        JFacePtr f1 = JQuadrilateral::newObject(v00, v01, v11, v10);
        JFacePtr f2 = JQuadrilateral::newObject(v01, v02, v12, v11);
        JFacePtr f3 = JQuadrilateral::newObject(v02, v03, v13, v12);
        JFacePtr f4 = JQuadrilateral::newObject(v03, v00, v10, v13);
        mesh->addObject(v10);
        mesh->addObject(v11);
        mesh->addObject(v12);
        mesh->addObject(v13);

        mesh->addObject(f1);
        mesh->addObject(f2);
        mesh->addObject(f3);
        mesh->addObject(f4);
        face->setStatus( JMeshEntity::REMOVE);
    }

    mesh->pruneAll();

    return mesh;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator ::collapse_matched_triangles(const JMeshPtr &triMesh,
        const vector<JFacePair> &matching)
{
    if( triMesh == nullptr ) return nullptr;

    JMeshPtr quadMesh = JMesh::newObject();

    JNodeSequence tnodes = triMesh->getNodes();
    quadMesh->addObjects(tnodes );

    size_t nSize = matching.size();
    quadMesh->reserve( nSize, 2 );

    JEdgeSequence commedges;
    for (size_t i = 0; i < nSize; i++) {
        const JFacePtr &tri1 = matching[i].first;
        const JFacePtr &tri2 = matching[i].second;
        JFacePtr quad = JQuadrilateral::newObject(tri1, tri2);
        quadMesh->addObject(quad);
        JFace::getSharedEntities(tri1, tri2, commedges);
        tri1->setStatus(JMeshEntity::REMOVE);
        tri2->setStatus(JMeshEntity::REMOVE);
        assert( commedges.size() == 1);
        commedges[0]->setStatus( JMeshEntity::REMOVE);
    }

    nSize = triMesh->getSize(2);
    for( size_t i = 0; i < nSize; i++) {
        const JFacePtr &tri = triMesh->getFaceAt(i);
        if( tri->isActive() ) quadMesh->addObject(tri);
    }

    triMesh->clearAll();

    return quadMesh;
}

///////////////////////////////////////////////////////////////////////////////

int AllQuadMeshGenerator ::verify_matching(const JMeshPtr &mesh, const vector<JFacePair> &matching)
{
    // This is the Graph matching verification. It is useful to verify the
    // matching, if done by other algorithms.

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++) {
        JFacePtr f = mesh->getFaceAt(i);
        f->setVisitBit(0);
    }

    JEdgeSequence edges;
    size_t nsize = matching.size();
    for (size_t i = 0; i < nsize; i++) {
        JFacePtr f0 = matching[i].first;
        JFacePtr f1 = matching[i].second;
        JFace::getSharedEntities( f0, f1, edges);
        if( edges.size() != 1 ) {
            cout << "Error: Face matching error: faces not neighbors  " << endl;
            exit(0);
            return 1;
        }
        assert(!f0->getVisitBit());
        assert(!f1->getVisitBit());
        f0->setVisitBit(1);
        f1->setVisitBit(1);
    }

    for (size_t i = 0; i < numfaces; i++) {
        JFacePtr f = mesh->getFaceAt(i);
        assert(f->getVisitBit());
    }

    return 0;

}

///////////////////////////////////////////////////////////////////////////////

int AllQuadMeshGenerator :: refineTriangle(const JFacePtr &oldface)
{
    JNodeSequence newNodes;
    JEdgeSequence newEdges;
    JFaceSequence newFaces;

    if( oldface == nullptr) return 1;
    if( !oldface->isActive() )  return 2;
    if( oldface->getTypeID() != JFace::TRIANGLE) return 3;

    string name = "Steiner";

    JNodePtr vcenter;
    Point3D pc;
    int    err;

    if( oldface->hasAttribute(name) ) {
        oldface->getAttribute(name, vcenter);
    } else {
        vcenter = JNode::newObject();
        oldface->getAvgXYZ( pc );
        vcenter->setXYZCoords( pc );
        newNodes.push_back(vcenter);
    }

    const JNodePtr &v0 = oldface->getNodeAt( 0 );
    const JNodePtr &v1 = oldface->getNodeAt( 1 );
    const JNodePtr &v2 = oldface->getNodeAt( 2 );

    JEdgePtr  subedge1, subedge2;
    int tag;

    JNodeSequence edgenodes(3);
    JEdgeSequence triedges(3);

    for( int i = 0; i < 3; i++) {
        triedges[i] = oldface->getEdgeAt(i);
        if( !triedges[i]->hasAttribute(name) ) {
            const JNodePtr &v0 = oldface->getNodeAt(i);
            const JNodePtr &v1 = oldface->getNodeAt(i+1);
            edgenodes[i] = JNodeGeometry::getMidNode(v0,v1);
            triedges[i]->setAttribute(name, edgenodes[i]);
            newNodes.push_back(edgenodes[i]);
            subedge1 = JSimplex::getEdgeOf(v0, edgenodes[i], 1);
            subedge2 = JSimplex::getEdgeOf(v1, edgenodes[i], 1);
            err  = triedges[i]->getAttribute("Boundary", tag);
            if( !err) {
                subedge1->setAttribute("Boundary", tag);
                subedge2->setAttribute("Boundary", tag);
            }
            newEdges.push_back(subedge1);
            newEdges.push_back(subedge2);
        }
        triedges[i]->getAttribute(name, edgenodes[i]);
    }

    newFaces.resize(3);
    newFaces[0] = JQuadrilateral::newObject(vcenter, edgenodes[2], v0, edgenodes[0]);
    newFaces[1] = JQuadrilateral::newObject(vcenter, edgenodes[0], v1, edgenodes[1]);
    newFaces[2] = JQuadrilateral::newObject(vcenter, edgenodes[1], v2, edgenodes[2]);

    newEdges.push_back(JSimplex::getEdgeOf(vcenter, edgenodes[0], 1));
    newEdges.push_back(JSimplex::getEdgeOf(vcenter, edgenodes[1], 1));
    newEdges.push_back(JSimplex::getEdgeOf(vcenter, edgenodes[2], 1));

    // Old face is removed ...
    oldface->setStatus( JMeshEntity::REMOVE);

    for( int i = 0;  i < 3; i++) {
        if( triedges[i]->getNumRelations(2) == 0)
            triedges[i]->setStatus( JMeshEntity::REMOVE);
    }

    if( quadMesh ) {
        quadMesh->addObjects(newNodes );
        quadMesh->addObjects(newEdges );
        quadMesh->addObjects(newFaces );
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator:: getSimpleTris2Quads()
{
    triMesh = mesh;

    if( triMesh == nullptr ) {
        cout << "Warning : Null object passed " << endl;
        return quadMesh;
    }

    if( triMesh->getTopology()->getElementsType(2) != JFace::TRIANGLE ) {
        cout << "Warning : All triangle mesh required" << endl;
        return quadMesh;
    }

    if( triMesh->hasAttribute("Steiner", 1) )
        cout << "Warning: Edge has steiner attribute: It will be replaced " << endl;

    triMesh->deleteEdgeAttribute( "Steiner");

    size_t numNodes = triMesh->getSize(0);
    size_t numEdges = triMesh->getSize(1);
    size_t numFaces = triMesh->getSize(2);

    quadMesh = JMesh::newObject();

    quadMesh->reserve( numNodes + numEdges + numFaces, 0 );
    quadMesh->reserve( 3*numFaces, 2 );

    size_t index = numNodes;
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = triMesh->getNodeAt(i);
        if( vtx->isActive() ) quadMesh->addObject(vtx);
    }

    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &tri = triMesh->getFaceAt(i);
        if( tri->isActive() ) refineTriangle(tri);
    }

    triMesh->deleteEdgeAttribute( "Steiner");
    triMesh->pruneAll();

    quadMesh->enumerate(0);
    quadMesh->enumerate(1);
    quadMesh->enumerate(2);

    quadMesh->getTopology()->searchBoundary();
    return quadMesh;
}
////////////////////////////////////////////////////////////////////////////////
JMeshPtr AllQuadMeshGenerator:: getSimpleQuadMesh()
{
    /*
        if( mesh  == nullptr )
            return 1;
        }
        vector<JEdgeSequence>  boundary;
        mesh->getTopology()->getBoundary(boundary);

        int numloops = boundary.size();
        for(int i = 0; i < numloops ; i++) {
            if( boundary.size()%2 ) {
                cout << "One of the boundary loop has odd number of segments" << endl;
                return nullptr;
            }
        }
    */
}

////////////////////////////////////////////////////////////////////////////////
JMeshPtr JQuadExtractor :: getQuadMesh()
{
    if( mesh ==  nullptr) return nullptr;

    JMeshOFFExporter mexp;
    mexp.writeFile(mesh, "xyz.off");

    JMeshPtr uvmesh;
    mesh->getAttribute("TextureMesh", uvmesh);

    if( uvmesh == nullptr) return nullptr;

    size_t numnodes = uvmesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = uvmesh->getNodeAt(i);
        Point3D uv = vtx->getXYZCoords();
        uv[0] *= uvscale;
        uv[1] *= uvscale;
        vtx->setXYZCoords(uv);
    }

    mexp.writeFile(uvmesh, "uv.off");

    numnodes = uvmesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = uvmesh->getNodeAt(i);
        Point3D uv = vtx->getXYZCoords();
        uv[0] /= uvscale;
        uv[1] /= uvscale;
        vtx->setXYZCoords(uv);
    }

    int err = system("quadextract  xyz.off uv.off uvquads.off");
    if( err < 0) return nullptr;

    JMeshOFFImporter mimp;
    JMeshPtr uvQuads = mimp.readFile("uvquads.off");
    return uvQuads;
}
////////////////////////////////////////////////////////////////////////////////
size_t JQuadExtractor :: getNumIntegerUV()
{
    JMeshPtr uvmesh;
    mesh->getAttribute("TextureMesh", uvmesh);
    if( uvmesh  == nullptr) return 0;

    set<pair<int, int> >  aset;
    size_t numnodes = uvmesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = uvmesh->getNodeAt(i);
        const Point2D  &uv  = vtx->getXYCoords();
        double u = fabs(uv[0]);
        double v = fabs(uv[1]);
        if( fabs(floor(u) - u) == 0.0 &&  fabs(floor(v)-v) == 0.0) {
            aset.insert(make_pair( (int)u, (int)v ));
        }
    }
    return aset.size();
}
////////////////////////////////////////////////////////////////////////////////

int AllQuadMeshGenerator:: getBinaryTreeMatching()
{
    cout << "Exit" << endl;
    exit(0);
    /*
        Jaal::JBinaryTreeMatch t2quad;
        JMeshPtr quadMesh = t2quad.getQuadMesh(triMesh);
        return quadMesh;
    */
}

////////////////////////////////////////////////////////////////////////////////


JMeshPtr  square_holey_pattern()
{

    int dim[2] = {20,20};
    double len[2] = {13,13};
    double org[2] = {0.0,0.0};

    /*
        JMeshPtr mesh  = AllQuadMeshGenerator::getStructuredMesh( dim,  len, org);

        int numfaces = mesh->getSize(2);

        int nrandom = 20;
        for( int i = 0; i < nrandom; i++) {
            const JFacePtr &face = mesh->getRandomFace();
            face->setStatus( JMeshEntity::REMOVE);
        }
        mesh->pruneFaces();
        return mesh;
    i*/
}

int AllQuadMeshGenerator :: splitBoundQuad2QuadTriangle(const JEdgePtr &boundedge)
{
    /*
    JFaceSequence faceneighs;
    JEdge::getRelations( boundedge, faceneighs);
    if( faceneighs.size() != 1) {
        cout << "warning: Only boundary edge should be refined " << endl;
        return 1;
    }

    int pos = faceneighs[0]->getPosOf( boundedge);
    if( pos < 0) {
        cout << "warning: boundary edge does not belong to the face" << endl;
        return 1;
    }

    const JNodePtr v0 = faceneighs[0]->getNodeAt(pos+0);
    const JNodePtr v1 = faceneighs[0]->getNodeAt(pos+1);
    const JNodePtr v2 = faceneighs[0]->getNodeAt(pos+2);
    const JNodePtr v3 = faceneighs[0]->getNodeAt(pos+3);

    JNodePtr  v4 = JNodeGeometry::getMidNode(v0,v1);

    JFacePtr qface = Quadrilateral::newObject( v0,v4,v2,v3);
    JFacePtr tface = Triangle::newObject( v4,v1,v2);

    mesh->addObject(v4);
    JEdgeSequence newedges(3);
    newedges[0] = JEdge::newObject(v0,v4);
    newedges[1] = JEdge::newObject(v4,v1);
    newedges[2] = JEdge::newObject(v4,v2);
    mesh->addObjects(newedges);

    int bid;
    int err = boundedge->getAttribute("Boundary", bid);
    if( !err) {
       newedges[0]->setAttribute("Boundary", bid);
       newedges[1]->setAttribute("Boundary", bid);
    }

    mesh->addObject(qface);
    mesh->addObject(tface);

    faceneighs[0]->setStatus(JMeshEntity::REMOVE);
    boundedge->setStatus(JMeshEntity::REMOVE);
    */

    return 0;
}

////////////////////////////////////////////////////////////////////////////////


int AllQuadMeshGenerator :: makeEvenTriangulation()
{
    /*
        size_t numfaces = triMesh->getActiveSize(2);
        if( numfaces%2 == 0) return 0;

        JEdgePtr boundary;
        triMesh->getTopology()->searchBoundary();
        boundary = triMesh->getRandomEdge(1);

        JFaceSequence faceneighs;
        JEdge::getRelations(boundary, faceneighs);
        assert( faceneighs.size() == 1 );

        JFacePtr btri  = faceneighs[0];
        const JNodePtr &bv0 = boundary->getNodeAt(0);
        const JNodePtr &bv1 = boundary->getNodeAt(1);
        const JNodePtr &bv2 = Triangle::getOppositeNode(btri, bv0, bv1);

        Point3D p3d;
        JNodeGeometry::getMidPoint(bv0, bv1, p3d);

        JNodePtr bound = JNode::newObject();
        bound->setXYZCoords(p3d);
        triMesh->addObject(bound);

        JFacePtr tri0 = Triangle::newObject( bv0, bound, bv2 );
        triMesh->addObject(tri0);

        JFacePtr tri1 = Triangle::newObject( bound, bv1, bv2 );
        triMesh->addObject(tri1);

        JEdgePtr edge1 =  tri0->getEdgeAt(0);
        JEdgePtr edge2 =  tri1->getEdgeAt(0);

        int bmark;
        int err = boundary->getAttribute("Boundary", bmark);

        if( !err) {
            edge1->setAttribute("Boundary", bmark);
            edge2->setAttribute("Boundary", bmark);
        }

        triMesh->addObject(edge1);
        triMesh->addObject(edge2);

        btri->setStatus( JMeshEntity::REMOVE );
        boundary->setStatus(JMeshEntity::REMOVE);

        triMesh->pruneFaces();
        triMesh->pruneEdges();
    */

    return 0;
}




