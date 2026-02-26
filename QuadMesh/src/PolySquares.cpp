#include  "PolySquares.hpp"

///////////////////////////////////////////////////////////////////////////
JPolySquares :: JPolySquares()
{
    valid_topology = 0;
    deformDirection = -1;
    limDeformer.reset( new JLocallyInjectiveMap);
}

///////////////////////////////////////////////////////////////////////////

void JPolySquares :: setMesh( const JMeshPtr &m)
{
    inmesh = m;
    if( m == nullptr )   {
        cout << "Warning: Empty object passed " << endl;
        return;
    }

    inmesh->getTopology()->searchBoundary();

}
/////////////////////////////////////////////////////////////////////////////////

void JPolySquares :: clear()
{
    singularGraph.reset();
    voxmesh.reset();
}

/////////////////////////////////////////////////////////////////////////////////
int  JPolySquares :: reorient( const JEdgePtr &edge)
{
    if( !edge->isActive() )   return 1;

    const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
    const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();

    double dx = p1[1] - p0[1];
    double dy = p0[0] - p1[0];

    double angle = 180.0*atan2(dy,dx)/M_PI;

    if( angle < 0.0) angle = 360.0 + angle;

    int side = -1;
    if( angle >= 0   && angle <= 45.0  ) side = 1;
    if( angle >= 315 && angle <= 360.0 ) side = 1;
    if( angle >= 45  && angle <= 135.0 ) side = 2;
    if( angle >= 135 && angle <= 225.0 ) side = 3;
    if( angle >= 225 && angle <= 315.0 ) side = 0;

    edge->setAttribute("SquareSide", side);

    return 0;
}

///////////////////////////////////////////////////////////////////////////

void JPolySquares :: boundarySegmentation()
{
    if( inmesh == nullptr) return;

    inmesh->deleteNodeAttribute("SingularNode");
    inmesh->deleteEdgeAttribute("SquareSide");

    Vec3F vec;

    edgeNormal.resize(6);

    vec[0] = 0.0;
    vec[1] =-1.0;
    vec[2] = 0.0;
    edgeNormal[0] = vec;

    vec[0] =  1.0;
    vec[1] =  0.0;
    vec[2] =  0.0;
    edgeNormal[1] = vec;

    vec[0] =  0.0;
    vec[1] =  1.0;
    vec[2] =  0.0;
    edgeNormal[2] = vec;

    vec[0] =  -1.0;
    vec[1] =   0.0;
    vec[2] =   0.0;
    edgeNormal[3] = vec;

    int dim = inmesh->getTopology()->getDimension();
    inmesh->getTopology()->searchBoundary();

    JEdgeSequence boundedges;
    if( dim == 2)
        inmesh->getTopology()->getBoundary(boundedges);
    else
        boundedges = inmesh->getEdges();

    size_t numedges = boundedges.size();
    assert( numedges );

    for( size_t i = 0; i < numedges; i++)
        reorient( boundedges[i] );
}

//////////////////////////////////////////////////////////////////////////////

void JPolySquares :: getNodesInBetween( const JNodeSequence &nodes,
                                        const JNodePtr &vstart, JNodePtr &vend,
                                        JNodeSequence &outnodes)
{
    outnodes.clear();
    int pos1 = -1;
    int pos2 = -1;
    int nsize = nodes.size();
    for( int i = 0; i <  nsize; i++) {
        if( nodes[i] == vstart) {
            pos1 = i;
            for ( int j = i+1; j < nsize+2; j++) {
                if( nodes[j%nsize] == vend ) {
                    pos2 = j%nsize;
                    break;
                }
            }
            break;
        }
    }
    assert( pos1 != -1 && pos2 != -1);
    for( int i = pos1; i <= pos2; i++) {
        const JNodePtr &vtx = nodes[i%nsize];
        outnodes.push_back( vtx );
    }
}
//////////////////////////////////////////////////////////////////////////////

void JPolySquares :: setLinePoints( const JNodeSequence &nodes)
{
    size_t vid;

    const JNodePtr &v0 = nodes.front();
    const JNodePtr &v1 = nodes.back();

    int nsize  = nodes.size();
    double  dt = 1.0/(double)(nsize-1);

    vid = v0->getID();
    const Point3D &p0 = polymesh->getNodeAt(vid)->getXYZCoords();

    vid = v1->getID();
    const Point3D &p1 = polymesh->getNodeAt(vid)->getXYZCoords();

    Point3D xyz;
    for( int i = 1; i < nsize-1; i++) {
        double t = i*dt;
        xyz[0] = (1-t)*p0[0] + t*p1[0];
        xyz[1] = (1-t)*p0[1] + t*p1[1];
        xyz[2] = (1-t)*p0[2] + t*p1[2];
        vid    = nodes[i]->getID();
        const JNodePtr &vtx = polymesh->getNodeAt(vid);
        vtx->setXYZCoords(xyz);
    }
}

//////////////////////////////////////////////////////////////////////////////
void JPolySquares :: integerSnap()
{
    int dim[2];
    double org[2];
    double len[2];

    size_t numSingular = singularGraph->getSize(0);

    if( numSingular == 0 ) return;

    vector<double> xCorner, yCorner;
    for( size_t i = 0; i < numSingular; i++) {
        const JNodePtr &vtx = singularGraph->getNodeAt(i);
        const Point3D &xyz = vtx->getXYZCoords();
        xCorner.push_back( xyz[0] );
        yCorner.push_back( xyz[1] );
        vtx->setAttribute("TargetPos", xyz);
    }

    // Processs X -Direction ...
    boost::sort( xCorner );
    xCorner.erase( std::unique( xCorner.begin(), xCorner.end()), xCorner.end() );

    int nSize = xCorner.size();
    vector<double> xRelativeDist(nSize);
    xRelativeDist[0] = 0.0;
    for( int i = 1; i < nSize; i++)
        xRelativeDist[i] = xCorner[i] - xCorner[i-1];

    vector<int> xInteger(nSize);
    xInteger[0] = floor( xCorner[0] );

    for( int i = 1; i < nSize; i++)
        xInteger[i] = ceil(xInteger[i-1] + xRelativeDist[i]);

    xmap.clear();
    for( int i = 0; i < nSize; i++) {
        xmap[xCorner[i]] = xInteger[i];
    }

    dim[0] = xInteger[nSize-1] - xInteger[0] + 1;
    len[0] = xInteger[nSize-1] - xInteger[0];
    org[0] = xInteger[0];

    // Process y-Direction ....
    boost::sort( yCorner );
    yCorner.erase( std::unique( yCorner.begin(), yCorner.end()), yCorner.end() );

    nSize = yCorner.size();
    vector<double> yRelativeDist(nSize);
    yRelativeDist[0] = 0.0;
    for( int i = 1; i < nSize; i++)
        yRelativeDist[i] = yCorner[i] - yCorner[i-1];

    vector<int> yInteger(nSize);
    yInteger[0] = floor( yCorner[0] );

    for( int i = 1; i < nSize; i++)
        yInteger[i] = ceil(yInteger[i-1] + yRelativeDist[i]);

    ymap.clear();
    for( int i = 0; i < nSize; i++) {
        ymap[yCorner[i]] = yInteger[i];
    }

    dim[1] = yInteger[nSize-1] - yInteger[0] + 1;
    len[1] = yInteger[nSize-1] - yInteger[0];
    org[1] = yInteger[0];

    for( size_t i = 0; i < numSingular; i++) {
        const JNodePtr &vtx = singularGraph->getNodeAt(i);
        Point3D xyz = vtx->getXYZCoords();
        xyz[0]  = xmap[ xyz[0]];
        xyz[1]  = ymap[ xyz[1]];
        vtx->setXYZCoords( xyz );
    }
    voxmesh = AllQuadMeshGenerator::getStructuredMesh( dim, len, org);
}

//////////////////////////////////////////////////////////////////////////////
int JPolySquares :: setOriginalPosition( const JNodePtr &vtx)
{
    double eps = 1.0E-06;
    Point3D xyz = vtx->getXYZCoords();

    JFacePtr face = pointLoc.searchFace(xyz,1);
    if( face == nullptr) {
        cout << "Error: the following point in not found in any triangle" << endl;
        cout << xyz[0] << " " << xyz[1] << "  " << xyz[2] << endl;
        return 1;
    }

    assert( JTriGeometry::isInside(face, xyz,1));

    Point3D uvw;
    JTriGeometry :: getBaryCoordinates( face, xyz, uvw);
    assert( uvw[0] > -eps && uvw[0] < 1.0 + eps);
    assert( uvw[1] > -eps && uvw[1] < 1.0 + eps);
    assert( uvw[2] > -eps && uvw[2] < 1.0 + eps);

    int err;
    Point3D p0, p1, p2;
    err = face->getNodeAt(0)->getAttribute("TargetPos", p0);
    assert(!err);
    err = face->getNodeAt(1)->getAttribute("TargetPos", p1);
    assert(!err);
    err = face->getNodeAt(2)->getAttribute("TargetPos", p2);
    assert(!err);

    JTriGeometry::getXYZCoordinates(&p0[0], &p1[0], &p2[0], &uvw[0], &xyz[0]);
    xyz[2] = 0.0;

    vtx->setXYZCoords(xyz);
    return 0;
}
//////////////////////////////////////////////////////////////////////////////
int JPolySquares :: setOriginalPosition()
{
    map<int,double>   imap, jmap;
    map<double,int>::iterator it;

    for( it = xmap.begin(); it != xmap.end(); ++it) {
        double key = it->first;
        int    val = it->second;
        imap[val] = key;
    }

    for( it = ymap.begin(); it != ymap.end(); ++it) {
        double key = it->first;
        int    val = it->second;
        jmap[val] = key;
    }

    map<JNodePtr, JNodePtr> src2dst;
    map<JNodePtr, JNodePtr> dst2src;

    size_t numnodes    = voxmesh->getSize(0);
    size_t numSingular = singularGraph->getSize(0);
    for( size_t i = 0; i < numSingular; i++) {
        const JNodePtr &ivtx = singularGraph->getNodeAt(i);
        const Point3D  &ixyz = ivtx->getXYZCoords();
        for( size_t j = 0; j < numnodes; j++) {
            const JNodePtr &jvtx = voxmesh->getNodeAt(j);
            Point3D  jxyz = jvtx->getXYZCoords();
            if( (ixyz[0] == jxyz[0])  && (ixyz[1] == jxyz[1]) ) {
                src2dst[ivtx] = jvtx;
                dst2src[jvtx] = ivtx;
            }
        }
    }

    //  Now duplicate the source singularity graph on the deformed shape.
    pointLoc.setMesh( singularGraph );

    for( size_t i = 0; i < numnodes; i++) {
        int err = setOriginalPosition( voxmesh->getNodeAt(i) );
        if( err ) return 1;
    }
    return 0;
}
//////////////////////////////////////////////////////////////////////////////

void JPolySquares :: removeOutsideVoxels()
{
    JPointLocation  pointloc;
    pointloc.setMesh( singularGraph);

    size_t numfaces = voxmesh->getSize(2);
    Point3D xyz;
    for( size_t i = 0; i <  numfaces; i++) {
        const JFacePtr &f = voxmesh->getFaceAt(i);
        f->getAvgXYZ(xyz);
        if( pointloc.searchFace(xyz,0) == nullptr) {
            f->setStatus(JMeshEntity::REMOVE);
        }
    }
    voxmesh->pruneFaces();

    size_t numedges = voxmesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = voxmesh->getEdgeAt(i);
        if( edge->getNumRelations(2) == 0)
            edge->setStatus(JMeshEntity::REMOVE);
    }
    voxmesh->pruneEdges();

    voxmesh->buildRelations(0,2);
    size_t numnodes = voxmesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = voxmesh->getNodeAt(i);
        if( vtx->getNumRelations(2) == 0)
            vtx->setStatus(JMeshEntity::REMOVE);
    }
    voxmesh->pruneNodes();
    voxmesh->enumerate(0);
}
//////////////////////////////////////////////////////////////////////////////
void JPolySquares :: smoothSingularGraph(JEdgeSequence &edges)
{
    /*
        vector<int>  cornerID;
        // Start the edgeloop from one of the singular node ....
        double xmin = std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();

        int pos = 0;
        for( const JNodePtr &vtx : singularNodes) {
             const Point3D &xy = vtx->getXYZCoords();
             if( xy[1] < ymin ) {
                 if( xy[0] < xmin ) {
                     xmin = xy[0];
                     ymin = xy[1];
                     pos++;
                 }
             }
        }

        // Calculate the distance between the singular nodes ...
        JNodePtr node_start_from = singularNodes[pos];
        JEdge::getChain(edges, node_start_from);

        numCorners = singularNodes.size();
        for( int i = 0; i < numCorners; i++) {
            const JNodePtr &v0 = singularNodes[i];
            const JNodePtr &v1 = singularNodes[(i+1)%numCorners];
            double len = JNodeGeometry::getLength(v0,v1);
            v0->setAttribute("ChainLength", len);
        }

        // Now place all the singular nodes on integer locations.
        integerSnap();

        JNodeSequence nodes;
        JEdge::getChainNodes( edges, nodes);

        // Collect all the nodes between two singular nodes and align them along the
        // axis. We use uniform distribution of points....
        JNodeSequence segnodes;
        for( int i = 0; i < numCorners-1; i++) {
            getNodesInBetween( nodes, singularNodes[i], singularNodes[i+1], segnodes);
            genPolyEdge( segnodes);
        }

        // The last segment is left out in the previous step. Complete it here...
        int numnodes = nodes.size();
        pos = -1;
        for( int i = 0; i < numnodes; i++) {
            if( nodes[i] == singularNodes[numCorners-1] )  {
                pos = i;
                break;
            }
        }
        segnodes.clear();
        for( int i = pos; i < numnodes; i++)
            segnodes.push_back( nodes[i] );
        segnodes.push_back( nodes[0] );
        genPolyEdge( segnodes);
    */
}
///////////////////////////////////////////////////////////////////////////////////////////

void JPolySquares :: cornerSingularGraph(JEdgeSequence &edges)
{
    int err;
    assert( !edges.empty() );

    assert( JEdgeTopology::isChain(edges)  );
    assert( JEdgeTopology::isClosed(edges) );

    JNodeSequence nodes;
    JEdgeTopology::getChainNodes(edges, nodes);
    int val = 0;
    for( const JNodePtr &vtx: nodes) {
        double angle = JNodeGeometry::getSpanAngleAt(vtx, ANGLE_IN_DEGREES);
        if( angle < 180 && angle > 80) {
            vtx->setAttribute("SingularNode", val);
        }

        if( angle > 180 && angle >260) {
            vtx->setAttribute("SingularNode", val);
        }
    }

}

//////////////////////////////////////////////////////////////////////////////
void JPolySquares :: getSingularGraph(JEdgeSequence &edges)
{
    cornerSingularGraph(edges);
//   smoothSingularGraph(edges);

    JNodeSequence nodes;
    JEdgeTopology::getChainNodes(edges, nodes);

    JNodeSequence snodes;

    for( const JNodePtr &vtx: nodes) {
        if( vtx->hasAttribute("SingularNode") )
            snodes.push_back(vtx);
    }

    int numSingular = snodes.size();
    if( numSingular == 0)  {
        cout << "Warning: there are no singular nodes on the boundary " << endl;
        return;
    }

    singularGraph->addObjects(snodes);
    for( int i = 0; i < numSingular; i++) {
        const JNodePtr &v0 = snodes[i];
        const JNodePtr &v1 = snodes[(i+1)%numSingular];
        JEdgePtr edge = JEdge::newObject(v0,v1);
        singularGraph->addObject(edge);
    }
}

//////////////////////////////////////////////////////////////////////////////

const JMeshPtr & JPolySquares :: getSingularGraph()
{
    assert( inmesh );

    if( singularGraph) return singularGraph;

    JDelaunayMesh2D  delmesher;
    if( inmesh->getSize(2) == 0) {
        delmesher.setMesh(inmesh);
        inmesh = delmesher.getSimpleMesh();
    }

    inmesh->buildRelations(0,2);

    // Find out how many loops ( or holes ) in the geometry ....
    vector<JEdgeSequence> edgeloops;
    inmesh->getTopology()->getBoundary( edgeloops );

    singularGraph = JMesh::newObject();

    // Build all the information for a given loop...
    int numloops = edgeloops.size();
    assert( numloops > 0);
    for( int i = 0; i < numloops; i++)
        getSingularGraph( edgeloops[i] );

    if( singularGraph->getSize(0) == 0) {
        cout << "Warning: there are no singular nodes in the mesh " << endl;
        return singularGraph;
    }

    if( singularGraph->getSize(1) == 0) {
        cout << "Warning: there are no singular edges in the mesh " << endl;
        return singularGraph;
    }

    inmesh->deleteSurfaceMesh();

    delmesher.setMesh(singularGraph);
    JMeshPtr m = delmesher.getSimpleMesh();
    if( m  ) singularGraph = m;

    return singularGraph;
}

//////////////////////////////////////////////////////////////////////////////

bool JPolySquares ::  isIntegerMapped( const JMeshPtr &mesh) const
{
    // Determine if the singularities of a polycube are integer mapped.
    double xmin =  std::numeric_limits<double>::max();
    double ymin =  std::numeric_limits<double>::max();
    double zmin =  std::numeric_limits<double>::max();

    int integer_mapped = 1;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &anode = mesh->getNodeAt(i);
        if( anode->hasAttribute("PartitionCorner") ) {
            Point3D &xyz = anode->getXYZCoords();
            if( fabs((int)xyz[0] -xyz[0]) > 1.0E-10 ) integer_mapped = 0;
            if( fabs((int)xyz[1] -xyz[1]) > 1.0E-10 ) integer_mapped = 0;
            if( fabs((int)xyz[2] -xyz[2]) > 1.0E-10 ) integer_mapped = 0;
            xmin = min( xmin, xyz[0] );
            ymin = min( ymin, xyz[1] );
            zmin = min( zmin, xyz[2] );
        }
    }

    return integer_mapped;
}
//////////////////////////////////////////////////////////////////////////////

const JMeshPtr& JPolySquares :: getQuadMesh()
{
    getSingularGraph();
    integerSnap();
    removeOutsideVoxels();
    setOriginalPosition();
    return voxmesh;
}

///////////////////////////////////////////////////////////////////////////

#ifdef CSV
int JPolySquares :: targetPositions( const JEdgeSequence &edges)
{

    if( edges.empty() ) return 1;

    JNodeSequence nodes;
    JEdge::getChainNodes(edges, nodes);

    size_t numnodes = nodes.size();

    double dt = 1.0/(double)( numnodes-1);

    const JNodePtr &v0 = edges.front()->getNodeAt(0);
    const Point3D  &p0 = v0->getXYZCoords();

    const JNodePtr &v1 = edges.back()->getNodeAt(1);
    const Point3D  &p1 = v1->getXYZCoords();

    Point3D p3d;
    for( size_t i = 1; i < numnodes-1; i++) {
        double t = i*dt;
        p3d[0] = (1-t)*p0[0] + t*p1[0];
        p3d[1] = (1-t)*p0[1] + t*p1[1];
        p3d[2] = (1-t)*p0[2] + t*p1[2];
        nodes[i]->setXYZCoords(p3d);
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////
int  JPolySquares :: deformSource( int dir)
{
    /*
      if( dir == 1)
           tetmesh1 = AllTetMeshGenerator::getQualityMesh(insurfmesh);
      else
           tetmesh1 = AllTetMeshGenerator::getQualityMesh(polysurfmesh2);

      JMeshPartitioner mp;
      mp.setMesh(tetmesh1);

      for( const JNodePtr &vtx: singularNodes) {
           int id = vtx->getID();
           const JNodePtr &vsrc = insurfmesh->getNodeAt(id);
           const JNodePtr &vdst = polysurfmesh->getNodeAt(id);
           vsrc->setXYZCoords( vdst->getXYZCoords();
           if( eform() ) return;
      }

      JEdgeSequence edges;

      for( int i = 0; i < numInterfaces; i++) {
           mp.getInterace(edges);
           Simplex::getNodes( edges, nodes);
      for( const JNodePtr &vtx: nodes) {
           int id = vtx->getID();
           const JNodePtr &vsrc = insurfmesh->getNodeAt(id);
           const JNodePtr &vdst = polysurfmesh->getNodeAt(id);
           vsrc->setXYZCoords( vdst->getXYZCoords();
           if( eform() ) return;
      }
      }

      JFaceSequence faces;
      for( int i = 0; i < numPatches; i++) {
           mp.getRegion(faces);
           Simplex::getNodes( faces, nodes);
      for( const JNodePtr &vtx: nodes) {
           int id = vtx->getID();
           const JNodePtr &vsrc = insurfmesh->getNodeAt(id);
           const JNodePtr &vdst = polysurfmesh->getNodeAt(id);
           vsrc->setXYZCoords( vdst->getXYZCoords();
           if( eform() ) return;
      }
      }
    */
}

///////////////////////////////////////////////////////////////////////////
#endif

