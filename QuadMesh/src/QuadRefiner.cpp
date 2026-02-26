#include "MeshRefine.hpp"
using namespace Jaal;

int JQuadRefiner::refine4( const JFacePtr &oldface)
{
    newNodes.clear();
    newEdges.clear();
    newFaces.clear();

    if( oldface == nullptr ) return 1;
    if( !oldface->isActive() ) return 2;
    if( oldface->getTypeID() != JFace::QUADRILATERAL ) return 3;

    string name = "Steiner";

    JEdgeSequence quadedges(4);
    JNodeSequence edgenodes(4);

    for( int i = 0; i < 4; i++) {
        quadedges[i] = oldface->getEdgeAt(i);
        if( !quadedges[i]->hasAttribute(name) ) {
            const JNodePtr &v0 =  oldface->getNodeAt(i);
            const JNodePtr &v1 =  oldface->getNodeAt(i+1);
            edgenodes[i] = JNodeGeometry::getMidNode(v0,v1);
            quadedges[i]->setAttribute(name, edgenodes[i]);
            newNodes.push_back( edgenodes[i] );
            newEdges.push_back( JSimplex::getEdgeOf(v0, edgenodes[i], 1));
            newEdges.push_back( JSimplex::getEdgeOf(v1, edgenodes[i], 1));
        }
        quadedges[i]->getAttribute( name, edgenodes[i]);
    }

    const JNodePtr &v0 = oldface->getNodeAt(0);
    const JNodePtr &v1 = oldface->getNodeAt(1);
    const JNodePtr &v2 = oldface->getNodeAt(2);
    const JNodePtr &v3 = oldface->getNodeAt(3);

    JNodePtr ev0 = edgenodes[0];
    JNodePtr ev1 = edgenodes[1];
    JNodePtr ev2 = edgenodes[2];
    JNodePtr ev3 = edgenodes[3];

    JNodePtr vcenter;
    if( oldface->hasAttribute(name) )
        oldface->getAttribute(name, vcenter);
    else {
        vcenter = JNode::newObject();
        Point3D pc;
        oldface->getAvgXYZ( pc );
        vcenter->setXYZCoords( pc );
        newNodes.push_back(vcenter);
    }

    newFaces.resize(4);
    newFaces[0] = JQuadrilateral::newObject( vcenter, ev3, v0, ev0);
    newFaces[1] = JQuadrilateral::newObject( vcenter, ev0, v1, ev1);
    newFaces[2] = JQuadrilateral::newObject( vcenter, ev1, v2, ev2);
    newFaces[3] = JQuadrilateral::newObject( vcenter, ev2, v3, ev3);

    newEdges.push_back( JSimplex::getEdgeOf(vcenter,ev0,1));
    newEdges.push_back( JSimplex::getEdgeOf(vcenter,ev1,1));
    newEdges.push_back( JSimplex::getEdgeOf(vcenter,ev2,1));
    newEdges.push_back( JSimplex::getEdgeOf(vcenter,ev3,1));

    /*
        JEdgePtr edge;
        int err, bmark;
        err = edge0->getAttribute("BoundaryTag", bmark);
        if( !err) {
            edge = newFaces[0]->getEdgeAt(2);
            edge->setAttribute("BoundaryTag", bmark);
            edge = newFaces[1]->getEdgeAt(1);
            edge->setAttribute("Boundary", bmark);
        }

        err = edge2->getAttribute("Boundary", bmark);
        if( !err) {
            edge = newFaces[1]->getEdgeAt(2);
            edge->setAttribute("Boundary", bmark);
            edge = newFaces[2]->getEdgeAt(1);
            edge->setAttribute("Boundary", bmark);
        }

        err = edge3->getAttribute("Boundary", bmark);
        if( !err) {
            edge = newFaces[2]->getEdgeAt(2);
            edge->setAttribute("Boundary", bmark);
            edge = newFaces[3]->getEdgeAt(1);
            edge->setAttribute("Boundary", bmark);
        }

        err = edge4->getAttribute("Boundary", bmark);
        if( !err) {
            edge = newFaces[3]->getEdgeAt(2);
            edge->setAttribute("Boundary", bmark);
            edge = newFaces[0]->getEdgeAt(1);
            edge->setAttribute("Boundary", bmark);
        }
    */

    // Old face is removed ...
    oldface->setStatus( JMeshEntity::REMOVE);

    for( int i = 0; i < 4; i++)  {
        if( quadedges[i]->getNumRelations(2) == 0)
            quadedges[i]->setStatus(JMeshEntity::REMOVE);
    }

    if( outmesh ) {
        outmesh->addObjects( newNodes);
        outmesh->addObjects( newEdges);
        outmesh->addObjects( newFaces);
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JQuadRefiner::refine5(const JFacePtr &oldface, JNodeSequence &vnew, JFaceSequence &fnew)
{
    int err = refine5(oldface);

    vnew = newNodes;
    fnew = newFaces;

    return err;
}
////////////////////////////////////////////////////////////////////////////////

int JQuadRefiner::refine5(const JFacePtr &oldface)
{
    newNodes.clear();
    newFaces.clear();

    if( oldface == nullptr ) return 1;
    if( !oldface->isActive() ) return 2;
    if( oldface->getTypeID() != JFace::QUADRILATERAL ) return 3;

    if (JQuadGeometry::isConvex(oldface))
        refine15_convex(oldface);
    else
        refine15_concave(oldface);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JQuadRefiner :: refine( const JFacePtr &oldface)
{
    newNodes.clear();
    newFaces.clear();

    if( oldface == nullptr ) return 1;
    if( !oldface->isActive() ) return 2;
    if( oldface->getTypeID() != JFace::QUADRILATERAL ) return 3;

    const JNodePtr &v0 = oldface->getNodeAt(0);
    const JNodePtr &v1 = oldface->getNodeAt(1);
    const JNodePtr &v2 = oldface->getNodeAt(2);
    const JNodePtr &v3 = oldface->getNodeAt(3);

    int  ori, nsize;
    JNodeSequence nodes, qnodes;
    JNodeSequence boundnodes[4];

    for( int j = 0; j < 4; j++) {
        const JEdgePtr &edge = oldface->getEdgeAt(j);
        if( edge->hasAttribute("Steiner") ) {
            edge->getAttribute("Steiner", nodes);
            nsize  =  nodes.size();
            boundnodes[j].resize(nsize+2);
            boundnodes[j][0] = edge->getNodeAt(0);
            for( int i = 0; i < nsize; i++)
                boundnodes[j][i+1] = nodes[i];
            boundnodes[j][nsize+1] = edge->getNodeAt(1);
        } else {
            boundnodes[j].resize(2);
            boundnodes[j][0] = edge->getNodeAt(0);
            boundnodes[j][1] = edge->getNodeAt(1);
        }
        ori = oldface->getOrientation( edge );
        if( ori < 0) reverse(boundnodes[j].begin(), boundnodes[j].end() );
    }

    if( boundnodes[0].size() != boundnodes[2].size() ) {
        cout << "Warning: Opposite edge has different# of nodes" << endl;
        return 1;
    }

    // #of nodes on left side and right side must be equal...
    if( boundnodes[1].size() != boundnodes[3].size() ) {
        cout << "Warning: Opposite edge has different# of nodes" << endl;
        return 1;
    }

    int nx = boundnodes[0].size();
    if( nx < 2) return 2;

    int ny = boundnodes[1].size();
    if( ny < 2) return 2;


    qnodes.resize(nx*ny);

    qnodes[0] = v0;
    qnodes[nx-1] = v1;
    qnodes[nx*ny-1] = v2;
    qnodes[nx*ny-nx] = v3;

    for( int i = 1; i < nx-1; i++) {
        qnodes[i] = boundnodes[0][i];
        int offset = (ny-1)*nx + i;
        qnodes[offset] = boundnodes[2][i];
    }

    for( int j = 1; j < ny-1; j++) {
        int offset = j*nx + nx -1;
        qnodes[offset] = boundnodes[1][j];
        qnodes[j*nx] = boundnodes[3][j];
    }

    // Get all the internal nodes using TFI....
    for( int j = 1; j < ny-1; j++) {
        for( int i = 1; i < nx-1; i++) {
            const JNodePtr &vtx = JNode::newObject();
            qnodes[j*nx + i] = vtx;
            set_tfi_coords(i, j, nx, ny, qnodes);
            newNodes.push_back( vtx );
        }
    }

    if( outmesh ) outmesh->addObjects( newNodes);

    // Once all the nodes have been generated, create Quad faces ....
    int numfaces = (nx-1)*(ny-1);
    newFaces.resize(numfaces);
    int index = 0;
    for( int j = 0; j < ny-1; j++) {
        for( int i = 0; i < nx-1; i++) {
            int offset = j*nx + i;
            JNodePtr v0 = qnodes[offset];
            JNodePtr v1 = qnodes[offset+1];
            offset = (j+1)*nx + i;
            JNodePtr v3 = qnodes[offset];
            JNodePtr v2 = qnodes[offset+1];
            JFacePtr f = JQuadrilateral::newObject( v0, v1, v2, v3);
            newFaces[index++] =  f;
        }
    }

    if( outmesh ) outmesh->addObjects( newFaces);

    int xcells = nx-1;
    int ycells = ny-1;

    JEdgePtr edge;
    int err, bmark;
    err =oldface->getEdgeAt(0)->getAttribute("Boundary", bmark);
    if( !err ) {
        for( int i = 0; i < xcells; i++) {
            JFacePtr face = newFaces[i];
            edge = face->getEdgeAt(0);
            edge->setAttribute("Boundary", bmark);
        }
    }
    err = oldface->getEdgeAt(2)->getAttribute("Boundary", bmark);
    if( !err ) {
        for( int i = 0; i < xcells; i++) {
            JFacePtr face = newFaces[xcells*ycells- i-1];
            edge = face->getEdgeAt(2);
            edge->setAttribute("Boundary", bmark);
        }
    }

    err = oldface->getEdgeAt(1)->getAttribute("Boundary", bmark);
    if( !err ) {
        for( int j = 0; j < ycells; j++) {
            JFacePtr face = newFaces[j*xcells+ xcells-1];
            edge = face->getEdgeAt(1);
            edge->setAttribute("Boundary", bmark);
        }
    }
    err = oldface->getEdgeAt(3)->getAttribute("Boundary", bmark);
    if( !err ) {
        for( int j = 0; j < ycells; j++) {
            JFacePtr face = newFaces[j*xcells];
            edge = face->getEdgeAt(3);
            edge->setAttribute("Boundary", bmark);
        }
    }

    oldface->setStatus( JMeshEntity::REMOVE);

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int JQuadRefiner::refine( const JFacePtr &face, int *dim)
{
    if( !face->isActive() || dim == nullptr) return 1;

    int nx = dim[0];
    int ny = dim[1];
    if( nx < 2 || ny < 2 ) return 1;

    const JEdgePtr &edge0 = face->getEdgeAt(0);
    if( !edge0->hasAttribute("Steiner") )
        JEdgeGeometry::generateLinearNodes(edge0, nx, mesh);

    const JEdgePtr &edge1 = face->getEdgeAt(1);
    if( !edge1->hasAttribute("Steiner") )
        JEdgeGeometry::generateLinearNodes(edge1, ny, mesh);

    const JEdgePtr &edge2 = face->getEdgeAt(2);
    if( !edge2->hasAttribute("Steiner") )
        JEdgeGeometry::generateLinearNodes(edge2, nx, mesh);

    const JEdgePtr &edge3 = face->getEdgeAt(3);
    if( !edge3->hasAttribute("Steiner") )
        JEdgeGeometry::generateLinearNodes(edge3, ny, mesh);
    refine(face);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JQuadRefiner ::refine15_concave(const JFacePtr &face)
{
    if( !face->isActive() ) return 1;

    newNodes.resize(4);
    newFaces.resize(5);

    JNodeSequence rotatedNodes(4), tnodes(3);

    Point3D xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    JNodeSequence nodes = face->getNodes();
    JQuadrilateral::quad_tessalate(nodes, rotatedNodes);
    xyz = JNodeGeometry::getMidPoint(rotatedNodes[0], rotatedNodes[2], 1.0 / 3.0);
    newNodes[0] = JNode::newObject();
    newNodes[0]->setXYZCoords(xyz);

    xyz = JNodeGeometry::getMidPoint(rotatedNodes[0], rotatedNodes[2], 2.0 / 3.0);
    newNodes[2] = JNode::newObject();
    newNodes[2]->setXYZCoords(xyz);

    JTriangle triface;
    tnodes[0] = rotatedNodes[0];
    tnodes[1] = rotatedNodes[1];
    tnodes[2] = rotatedNodes[2];
    triface.setNodes(tnodes);
    triface.getAvgXYZ( xyz );
    newNodes[1] = JNode::newObject();
    newNodes[1]->setXYZCoords(xyz);

    tnodes[0] = rotatedNodes[0];
    tnodes[1] = rotatedNodes[2];
    tnodes[2] = rotatedNodes[3];
    triface.setNodes(tnodes);
    triface.getAvgXYZ( xyz );
    newNodes[3] = JNode::newObject();
    newNodes[3]->setXYZCoords(xyz);

    nodes[0] = rotatedNodes[0];
    nodes[1] = rotatedNodes[1];
    nodes[2] = newNodes[1];
    nodes[3] = newNodes[0];
    newFaces[0] = JQuadrilateral::newObject();
    newFaces[0]->setNodes(nodes);

    nodes[0] = rotatedNodes[1];
    nodes[1] = rotatedNodes[2];
    nodes[2] = newNodes[2];
    nodes[3] = newNodes[1];
    newFaces[1] = JQuadrilateral::newObject();
    newFaces[1]->setNodes(nodes);

    nodes[0] = rotatedNodes[2];
    nodes[1] = rotatedNodes[3];
    nodes[2] = newNodes[3];
    nodes[3] = newNodes[2];
    newFaces[2] = JQuadrilateral::newObject();
    newFaces[2]->setNodes(nodes);

    nodes[0] = rotatedNodes[3];
    nodes[1] = rotatedNodes[0];
    nodes[2] = newNodes[0];
    nodes[3] = newNodes[3];
    newFaces[3] = JQuadrilateral::newObject();
    newFaces[3]->setNodes(nodes);

    nodes[0] = newNodes[0];
    nodes[1] = newNodes[1];
    nodes[2] = newNodes[2];
    nodes[3] = newNodes[3];
    newFaces[4] = JQuadrilateral::newObject();
    newFaces[4]->setNodes(nodes);

    // Old face is removed ...
    face->setStatus( JMeshEntity::REMOVE);

    if( outmesh ) {
        outmesh->addObjects( newNodes );
        outmesh->addObjects( newFaces );
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JQuadRefiner ::refine15_convex(const JFacePtr &face)
{
    newNodes.clear();
    newEdges.clear();
    newFaces.clear();

    newNodes.resize(4);
    newEdges.resize(8);
    newFaces.resize(5);

    if( !face->isActive() ) return 1;

    Point3D xyz;
    vector<double> weight;
    vector<double> xc(4), yc(4), zc(4);
    for (int i = 0; i < 4; i++) {
        JNodePtr v = face->getNodeAt(i);
        xyz = v->getXYZCoords();
        xc[i] = xyz[0];
        yc[i] = xyz[1];
        zc[i] = xyz[2];
    }

    double dl = 2.0 / 3.0;

    JQuadGeometry::bilinear_weights(-1.0 + dl, -1.0 + dl, weight);
    xyz[0] = JSimplex::linear_interpolation(xc, weight);
    xyz[1] = JSimplex::linear_interpolation(yc, weight);
    xyz[2] = JSimplex::linear_interpolation(zc, weight);
    newNodes[0] = JNode::newObject();
    newNodes[0]->setXYZCoords(xyz);

    JQuadGeometry::bilinear_weights(-1.0 + 2 * dl, -1.0 + dl, weight);
    xyz[0] = JSimplex::linear_interpolation(xc, weight);
    xyz[1] = JSimplex::linear_interpolation(yc, weight);
    xyz[2] = JSimplex::linear_interpolation(zc, weight);
    newNodes[1] = JNode::newObject();
    newNodes[1]->setXYZCoords(xyz);

    JQuadGeometry::bilinear_weights(-1.0 + 2 * dl, -1.0 + 2 * dl, weight);
    xyz[0] = JSimplex::linear_interpolation(xc, weight);
    xyz[1] = JSimplex::linear_interpolation(yc, weight);
    xyz[2] = JSimplex::linear_interpolation(zc, weight);
    newNodes[2] = JNode::newObject();
    newNodes[2]->setXYZCoords(xyz);

    JQuadGeometry::bilinear_weights(-1.0 + dl, -1.0 + 2.0 * dl, weight);
    xyz[0] = JSimplex::linear_interpolation(xc, weight);
    xyz[1] = JSimplex::linear_interpolation(yc, weight);
    xyz[2] = JSimplex::linear_interpolation(zc, weight);
    newNodes[3] = JNode::newObject();
    newNodes[3]->setXYZCoords(xyz);

    JNodeSequence nodes(4), oldNodes;
    oldNodes = face->getNodes();

    nodes[0] = oldNodes[0];
    nodes[1] = oldNodes[1];
    nodes[2] = newNodes[1];
    nodes[3] = newNodes[0];
    newFaces[0] = JQuadrilateral::newObject(nodes);

    nodes[0] = oldNodes[1];
    nodes[1] = oldNodes[2];
    nodes[2] = newNodes[2];
    nodes[3] = newNodes[1];
    newFaces[1] = JQuadrilateral::newObject(nodes);

    nodes[0] = oldNodes[2];
    nodes[1] = oldNodes[3];
    nodes[2] = newNodes[3];
    nodes[3] = newNodes[2];
    newFaces[2] = JQuadrilateral::newObject(nodes);

    nodes[0] = oldNodes[3];
    nodes[1] = oldNodes[0];
    nodes[2] = newNodes[0];
    nodes[3] = newNodes[3];
    newFaces[3] = JQuadrilateral::newObject(nodes);

    nodes[0] = newNodes[0];
    nodes[1] = newNodes[1];
    nodes[2] = newNodes[2];
    nodes[3] = newNodes[3];
    newFaces[4] = JQuadrilateral::newObject(nodes);

    newEdges.push_back( JSimplex::getEdgeOf(newNodes[0], oldNodes[0], 1));
    newEdges.push_back( JSimplex::getEdgeOf(newNodes[1], oldNodes[1], 1));
    newEdges.push_back( JSimplex::getEdgeOf(newNodes[2], oldNodes[2], 1));
    newEdges.push_back( JSimplex::getEdgeOf(newNodes[3], oldNodes[3], 1));

    newEdges.push_back( JSimplex::getEdgeOf(newNodes[0], newNodes[1], 1));
    newEdges.push_back( JSimplex::getEdgeOf(newNodes[1], newNodes[2], 1));
    newEdges.push_back( JSimplex::getEdgeOf(newNodes[2], newNodes[3], 1));
    newEdges.push_back( JSimplex::getEdgeOf(newNodes[3], newNodes[0], 1));

    if( outmesh ) {
        outmesh->addObjects( newNodes );
        outmesh->addObjects( newEdges );
        outmesh->addObjects( newFaces );
    }

    // Old face is removed ...
    face->setStatus( JMeshEntity::REMOVE);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JQuadRefiner :: refineAll( int scheme )
{
    if( mesh == nullptr) return 1;

    for( int iter = 0; iter < numIters; iter++)  {
        size_t numFaces = mesh->getSize(2);
        if( scheme == QUAD14 ) {
            for( size_t i = 0; i < numFaces; i++) {
                const JFacePtr &face = mesh->getFaceAt(i);
                refine4(face);
            }
            mesh->deleteEdgeAttribute("Steiner");
        }

        if( scheme == QUAD15 ) {
            for( size_t i = 0; i < numFaces; i++) {
                const JFacePtr &face = mesh->getFaceAt(i);
                refine5(face);
            }
        }
    }

    mesh->pruneAll();
    mesh->getTopology()->searchBoundary();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JQuadRefiner :: refineAll( const JFaceSequence &fselected, int scheme )
{
    if( mesh == nullptr) return 1;

    for( int iter = 0; iter < numIters; iter++)  {
        size_t numFaces = mesh->getSize(2);
        if( scheme == QUAD14 ) {
            for( const JFacePtr &face : fselected)
                refine4(face);
            mesh->deleteEdgeAttribute("Steiner");
        }
        if( scheme == QUAD15 ) {
            for( const JFacePtr &face : fselected)
                refine5(face);
        }
    }

    mesh->pruneAll();
    mesh->getTopology()->searchBoundary();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////


int JQuadRefiner :: refineAll( int *dim)
{
    size_t numFaces = mesh->getSize(2);

    mesh->deleteEdgeAttribute("Steiner");

    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        refine( face, dim);
    }
    mesh->deleteEdgeAttribute("Steiner");
    mesh->pruneAll();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JQuadRefiner :: refine2( const JFacePtr &oldface)
{
    if( !oldface->isActive() ) return 1;

    const JNodePtr &v0 = oldface->getNodeAt(0);
    const JNodePtr &v1 = oldface->getNodeAt(1);
    const JNodePtr &v2 = oldface->getNodeAt(2);
    const JNodePtr &v3 = oldface->getNodeAt(3);

    const JEdgePtr &edge0 = oldface->getEdgeAt(0);
    const JEdgePtr &edge1 = oldface->getEdgeAt(1);
    const JEdgePtr &edge2 = oldface->getEdgeAt(2);
    const JEdgePtr &edge3 = oldface->getEdgeAt(3);

    JNodePtr v4, v5;
    JNodePtr v6 = JFaceGeometry::getCentroid(v0, v1, v2, v3);

    int attrib[4];

    attrib[0] = edge0->hasAttribute("Steiner");
    attrib[1] = edge1->hasAttribute("Steiner");
    attrib[2] = edge2->hasAttribute("Steiner");
    attrib[3] = edge3->hasAttribute("Steiner");

    if( attrib[0] && attrib[1] ) {
        edge0->getAttribute( "Steiner", v4 );
        edge1->getAttribute( "Steiner", v5 );
        newNodes.push_back(v6);
        newFaces.resize(3);
        newFaces[0] = JQuadrilateral::newObject(v0, v4, v6, v3);
        newFaces[1] = JQuadrilateral::newObject(v1, v5, v6, v4);
        newFaces[2] = JQuadrilateral::newObject(v2, v3, v6, v5);
    }

    if( attrib[1] && attrib[2] ) {
        edge1->getAttribute( "Steiner", v4 );
        edge2->getAttribute( "Steiner", v5 );
        newNodes.push_back(v6);
        newFaces.resize(3);
        newFaces[0] = JQuadrilateral::newObject(v0, v1, v4, v6);
        newFaces[1] = JQuadrilateral::newObject(v2, v5, v6, v4);
        newFaces[2] = JQuadrilateral::newObject(v3, v0, v6, v5);
    }

    if( attrib[2] && attrib[3] ) {
        edge2->getAttribute( "Steiner", v4 );
        edge3->getAttribute( "Steiner", v5 );
        newNodes.push_back(v6);
        newFaces.resize(3);
        newFaces[0] = JQuadrilateral::newObject(v0, v1, v6, v5);
        newFaces[1] = JQuadrilateral::newObject(v1, v2, v4, v6);
        newFaces[2] = JQuadrilateral::newObject(v3, v5, v6, v4);
    }


    if( attrib[0] && attrib[3] ) {
        edge0->getAttribute( "Steiner", v4 );
        edge3->getAttribute( "Steiner", v5 );
        newNodes.push_back(v6);
        newFaces.resize(3);
        newFaces[0] = JQuadrilateral::newObject(v0, v4, v6, v5);
        newFaces[1] = JQuadrilateral::newObject(v1, v2, v6, v4);
        newFaces[2] = JQuadrilateral::newObject(v2, v3, v5, v6);
    }

    if( outmesh ) {
        outmesh->addObjects( newNodes );
        outmesh->addObjects( newFaces );
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JQuadRefiner :: refine3( const JFacePtr &face )
{
    if( !face->isActive() ) return 1;

    JNodeSequence edgenodes;
    const JEdgePtr &edge0 = face->getEdgeAt(0);
    const JEdgePtr &edge1 = face->getEdgeAt(1);
    const JEdgePtr &edge2 = face->getEdgeAt(2);
    const JEdgePtr &edge3 = face->getEdgeAt(3);

    int ncount = 0;
    if( edge0->hasAttribute("Steiner" ) ) {
        edge0->getAttribute( "Steiner", edgenodes );
        if( edgenodes.size() != 2 ) return 1;
        ncount++;
    }

    if( edge1->hasAttribute("Steiner" ) ) {
        edge1->getAttribute( "Steiner", edgenodes );
        if( edgenodes.size() != 2 ) return 1;
        ncount++;
    }

    if( edge2->hasAttribute("Steiner" ) ) {
        edge2->getAttribute( "Steiner", edgenodes );
        if( edgenodes.size() != 2 ) return 1;
        ncount++;
    }

    if( edge3->hasAttribute("Steiner" ) ) {
        edge3->getAttribute( "Steiner", edgenodes );
        if( edgenodes.size() != 2 ) return 1;
        ncount++;
    }

    // If all the edges have been refined. We get simple structured mesh..
    if( ncount == 4 ) {
    }

    // If only two oppsite edges have been refined. then also remeshing is simple...
    if( ncount == 2 ) {
    }

    // Use the general case....
    const JNodePtr &v0 = face->getNodeAt(0);
    const JNodePtr &v1 = face->getNodeAt(1);
    const JNodePtr &v2 = face->getNodeAt(2);
    const JNodePtr &v3 = face->getNodeAt(3);

    JNodePtr v4, v5;
    JNodePtr v6 = JNode::newObject();
    JNodePtr v7 = JNode::newObject();

    JFacePtr f0, f1, f2, f3;

    Point3D  xyz;
    Point2D paramCoords;
    if( edge0 ->hasAttribute("Steiner") ) {
        edge0->getAttribute( "Steiner", edgenodes );

        int ori = face->getOrientation(edge0);
        if( ori == 1) {
            v4 =  edgenodes[0];
            v5 =  edgenodes[1];
        } else {
            v4 =  edgenodes[1];
            v5 =  edgenodes[0];
        }

        paramCoords[0] = 0.5;
        paramCoords[1] = 0.0;
        JQuadGeometry::getBilinearCoords( face,  paramCoords, xyz);
        v6->setXYZCoords(xyz);

        paramCoords[0] = -0.5;
        paramCoords[1] =  0.0;
        JQuadGeometry::getBilinearCoords( face, paramCoords, xyz);
        v7->setXYZCoords(xyz);

        newNodes.push_back(v6);
        newNodes.push_back(v7);

        f0 = JQuadrilateral::newObject(v0, v4, v7, v3);
        f1 = JQuadrilateral::newObject(v1, v2, v6, v5);
        f2 = JQuadrilateral::newObject(v2, v3, v7, v6);
        f3 = JQuadrilateral::newObject(v4, v5, v6, v7);

        refine3( f0 );
        refine3( f1 );
        refine3( f2 );

        newFaces.push_back(f3);

        if( outmesh) {
            outmesh->addObjects( newNodes );
            outmesh->addObjects( newFaces );
        }

        return 0;
    }

    if( edge1 ->hasAttribute("Steiner") ) {
        edge1->getAttribute( "Steiner", edgenodes );
        int ori = face->getOrientation(edge1);
        if( ori == 1) {
            v4 =  edgenodes[0];
            v5 =  edgenodes[1];
        } else {
            v4 =  edgenodes[1];
            v5 =  edgenodes[0];
        }

        paramCoords[0] = 0.0;
        paramCoords[1] = 0.5;
        JQuadGeometry::getBilinearCoords( face,  paramCoords, xyz);
        v6->setXYZCoords(xyz);

        paramCoords[0] =  0.0;
        paramCoords[1] = -0.5;
        JQuadGeometry::getBilinearCoords( face, paramCoords, xyz);
        v7->setXYZCoords(xyz);

        newNodes.push_back(v6);
        newNodes.push_back(v7);
        f0 = JQuadrilateral::newObject(v0, v1, v4, v7);
        f1 = JQuadrilateral::newObject(v2, v3, v6, v5);
        f2 = JQuadrilateral::newObject(v0, v7, v6, v3);
        f3 = JQuadrilateral::newObject(v4, v5, v6, v7);

        refine3( f0);
        refine3( f1);
        refine3( f2);
        newFaces.push_back( f3 );

        return 0;
    }

    if( edge2 ->hasAttribute("Steiner") ) {
        edge2->getAttribute( "Steiner", edgenodes );
        int ori = face->getOrientation(edge2);
        if( ori == 1) {
            v4 =  edgenodes[0];
            v5 =  edgenodes[1];
        } else {
            v4 =  edgenodes[1];
            v5 =  edgenodes[0];
        }

        paramCoords[0] = -0.5;
        paramCoords[1] =  0.0;
        JQuadGeometry::getBilinearCoords( face,  paramCoords, xyz);
        v6->setXYZCoords(xyz);

        paramCoords[0] = +0.5;
        paramCoords[1] =  0.0;
        JQuadGeometry::getBilinearCoords( face, paramCoords, xyz);
        v7->setXYZCoords(xyz);

        newNodes.push_back(v6);
        newNodes.push_back(v7);

        f0 =  JQuadrilateral::newObject(v0, v1, v7, v6);
        f1 =  JQuadrilateral::newObject(v2, v4, v7, v1);
        f2 =  JQuadrilateral::newObject(v3, v0, v6, v5);
        f3 =  JQuadrilateral::newObject(v4, v5, v6, v7);

        refine3( f0);
        refine3( f1);
        refine3( f2);

        newFaces.push_back(f3);

        return 0;
    }

    if( edge3 ->hasAttribute("Steiner") ) {
        edge3->getAttribute( "Steiner", edgenodes );
        int ori = face->getOrientation(edge3);
        if( ori == 1) {
            v4 =  edgenodes[0];
            v5 =  edgenodes[1];
        } else {
            v4 =  edgenodes[1];
            v5 =  edgenodes[0];
        }

        paramCoords[0] = 0.0;
        paramCoords[1] = 0.5;
        JQuadGeometry::getBilinearCoords( face,  paramCoords, xyz);
        v6->setXYZCoords(xyz);

        paramCoords[0] =  0.0;
        paramCoords[1] = -0.5;
        JQuadGeometry::getBilinearCoords( face, paramCoords, xyz);
        v7->setXYZCoords(xyz);

        newNodes.push_back(v6);
        newNodes.push_back(v7);

        f0 =  JQuadrilateral::newObject(v0, v1, v7, v5);
        f1 =  JQuadrilateral::newObject(v1, v2, v6, v7);
        f2 =  JQuadrilateral::newObject(v2, v3, v4, v6);
        f3 =  JQuadrilateral::newObject(v4, v5, v7, v6);

        refine3( f0);
        refine3( f1);
        refine3( f2);

        newFaces.push_back(f3);

        return 0;
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int JQuadRefiner :: refineAt( const JNodePtr &vertex )
{
    if( mesh == nullptr ) return 1;
    if( mesh->getAdjTable(0,2) == 0)
        mesh->buildRelations(0,2);

    mesh->deleteEdgeAttribute("Steiner");

    JFaceSequence faceneighs;
    JNode::getRelations( vertex, faceneighs );
    int nSize = faceneighs.size();

    JEdgeSequence edges;
    JNodePtr vmid;
    for( int i = 0; i < nSize; i++) {
        int  pos = faceneighs[i]->getPosOf(vertex);
        const JEdgePtr &edge0 = faceneighs[i]->getEdgeAt(i+pos);
        const JEdgePtr &edge1 = faceneighs[i]->getEdgeAt(nSize-1-pos);
        assert( edge0 );
        assert( edge1 );
        if( !edge0->hasAttribute("Steiner")  ) {
            vmid = JNodeGeometry::getMidNode( edge0->getNodeAt(0), edge0->getNodeAt(1) ) ;
            edge0->setAttribute( "Steiner", vmid);
            outmesh->addObject( vmid );
        }

        if( !edge1->hasAttribute("Steiner")  ) {
            vmid = JNodeGeometry::getMidNode( edge1->getNodeAt(0), edge1->getNodeAt(1) ) ;
            edge1->setAttribute( "Steiner", vmid);
            outmesh->addObject( vmid );
        }
    }

    for( int i = 0; i < nSize; i++)
        refine2(faceneighs[i]);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JQuadRefiner::insert_boundary_pillows()
{
    if( mesh == nullptr ) return 1;

    mesh->getTopology()->searchBoundary();
    mesh->buildRelations(0,2);

    map<JNodePtr, JNodePtr> dupnode;
    JNodeSequence neighs;

    JEdgeSequence boundEdges;
    mesh->getTopology()->getBoundary(boundEdges);
    if( boundEdges.empty() ) {
        cout << "Warning: There are no boundary edges detected: inserting pillows failed " << endl;
        return 1;
    }

    JNodeSequence boundNodes;
    mesh->getTopology()->getBoundary(boundNodes);

    if( boundNodes.empty() ) {
        cout << "Warning: There are no boundary nodes detected: inserting pillows failed " << endl;
        return 1;
    }

    /*  Probabky wrong here. Delete it...
        int grp;
        JNodeSequence fixedNodes;
        size_t numnodes = mesh->getSize(0);
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( !vtx->hasAttribute("Constraint") ) {
                vtx->setAttribute("Constraint", grp);
                fixedNodes.push_back(vtx);
            }
        }
    */

    JNodeSequence newnodes;

    newnodes.reserve( boundNodes.size() );
    newfaces.reserve( boundNodes.size() );

    JNodeSequence concaveCorners;
    for( size_t i = 0; i < boundNodes.size(); i++) {
        JNodePtr v = JNode::newObject();
        Point3D xyz = boundNodes[i]->getXYZCoords();
        v->setXYZCoords( xyz );
        dupnode[ boundNodes[i] ] = v;
        newnodes.push_back(v);
        double angle  = JNodeGeometry::getSpanAngleAt(boundNodes[i], ANGLE_IN_DEGREES);
        if( angle > 181.0) concaveCorners.push_back(v);
    }

    outmesh->addObjects(newnodes);
    outmesh->enumerate(0);

    JFaceSequence faceneighs;
    JFaceSet oldSet;
    for( size_t i = 0; i < boundNodes.size(); i++) {
        JNode::getRelations( boundNodes[i], faceneighs );
        for( const JFacePtr &face:faceneighs) oldSet.insert( face );
    }

    map<JNodePtr, JNodePtr>::iterator it;

    JFaceSequence newfaces;

    JNodeSequence facenodes;
    for( const JFacePtr &oldface: oldSet) {
        JFacePtr newface = oldface->getClone();
        facenodes = oldface->getNodes();
        for( size_t i = 0; i < facenodes.size(); i++) {
            const JNodePtr &oldvertex = oldface->getNodeAt(i);
            it = dupnode.find(oldvertex);
            if( it != dupnode.end() ) facenodes[i] = it->second;
        }
        newface->setNodes( facenodes);
        newfaces.push_back( newface);
    }

    for( size_t i = 0; i < boundEdges.size(); i++) {
        JEdge::getRelations( boundEdges[i], faceneighs );
        assert( faceneighs.size() == 1);
        const JNodePtr &v0 = boundEdges[i]->getNodeAt(0);
        const JNodePtr &v1 = boundEdges[i]->getNodeAt(1);
        const JNodePtr &v2 = dupnode[v1];
        const JNodePtr &v3 = dupnode[v0];
        assert(v2);
        assert(v3);
        JQuadrilateralPtr q = JQuadrilateral::newObject(v0,v1,v2,v3);
        newfaces.push_back(q);
    }
    outmesh->addObjects(newfaces);

    JEdgeSequence faceedges;
    int nedges;
    for( const JFacePtr oldface: oldSet) {
        nedges = oldface->getSize(0);
        faceedges.resize(nedges);
        for( int i = 0; i < nedges; i++)
            faceedges[i] = oldface->getEdgeAt(i);
        oldface->setStatus( JMeshEntity::REMOVE);
        for( const JEdgePtr &edge : faceedges) {
            if( edge->getNumRelations(2) == 0)
                edge->setStatus(JMeshEntity::REMOVE);
        }
    }

    mesh->pruneAll();
    mesh->enumerate(1);
    mesh->enumerate(2);

    /*
        for( const JNodePtr &v : concaveCorners) {
             v->setAttribute("Constraint", grp);
        }

        JLaplaceMeshSmoother mopt;
        mopt.setMesh(mesh);
        mopt.setNumIterations(100);
        mopt.smoothAll();

        for( const JNodePtr &v : concaveCorners)
             v->deleteAttribute("Constraint");

        JQuadMeshOptimizer mopt;
        mopt.setMesh(mesh);
        mopt.setNumIterations(1);
        mopt.smooth( newnodes );
    */

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

