#include  "PolyCubes.hpp"

///////////////////////////////////////////////////////////////////////////
JPolyCubes :: JPolyCubes()
{
    numInterfaces  = 0;
    numPatches     = 0;
    valid_topology = 0;
    deformDirection = -1;
}
///////////////////////////////////////////////////////////////////////////

void JPolyCubes :: setModelMesh( const JMeshPtr &m)
{
    if( m == nullptr) return;

    insurfmesh.reset();
    tetmesh.reset();

    limDeformer.reset( new JLocallyInjectiveMap);

    int topDim = m->getTopology()->getDimension();
    if( topDim == 2) insurfmesh = m;

    if( topDim == 3) {
        tetmesh = m;
        insurfmesh = tetmesh->getTopology()->getSurfaceMesh();
        limDeformer->setMesh(tetmesh);
    }
    valid_topology = 0;

    if( deformDirection == 1) {
        sourcemesh = insurfmesh;
        targetmesh = polysurfmesh;
    } else {
        sourcemesh = polysurfmesh;
        targetmesh = insurfmesh;
    }
}

///////////////////////////////////////////////////////////////////////////

int  JPolyCubes :: reorient( const JFacePtr &face)
{
    if( !face->isActive() ) return 1;

    Vec3F normal;
    int err = face->getAttribute("Normal", normal);
    if( err ) {
        cout << "Warning: Normal on the face is unknown " << endl;
        return 2;
    }

    int side = 6;
    face->setAttribute("CubeSide", side);

    double xn = fabs(normal[0]);
    double yn = fabs(normal[1]);
    double zn = fabs(normal[2]);

    if( xn > 0.8 && xn > yn && xn > zn) {
        if( normal[0] >  0.0 ) {
            side = JHexahedron::RIGHT_SIDE;
            face->setAttribute( "CubeSide", side);
        }
        if( normal[0] <  0.0 ) {
            side = JHexahedron::LEFT_SIDE;
            face->setAttribute( "CubeSide", side);
        }
    }

    if( yn > 0.8 && yn > xn && yn > zn ) {
        if( normal[1] >  0.0 ) {
            side = JHexahedron::TOP_SIDE;
            face->setAttribute( "CubeSide", side);
        }
        if( normal[1] <  0.0 ) {
            side = JHexahedron::BOTTOM_SIDE;
            face->setAttribute( "CubeSide", side);
        }
    }

    if( zn > 0.8 && zn > xn && zn > yn ) {
        if( normal[2] >  0.0 ) {
            side = JHexahedron::FRONT_SIDE;
            face->setAttribute( "CubeSide", side);
        }
        if( normal[2] <  0.0 ) {
            side = JHexahedron::BACK_SIDE;
            face->setAttribute( "CubeSide", side);
        }
    }

    face->setAttribute("Partition", side);

    return 0;
}

///////////////////////////////////////////////////////////////////////////

void JPolyCubes :: seedClusters( const JMeshPtr &amesh)
{
    if( amesh == nullptr) return;

    amesh->deleteNodeAttribute("PartitionCorner");
    amesh->deleteEdgeAttribute("Interface");
    amesh->deleteEdgeAttribute("CubeEdge");
    amesh->deleteFaceAttribute("CubeSide");
    amesh->deleteFaceAttribute("Partition");

    Vec3F vec;

    cubeNormal.resize(6);

    vec[0] = 1.0;
    vec[1] = 0.0;
    vec[2] = 0.0;
    cubeNormal[0] = vec;

    vec[0] = -1.0;
    vec[1] =  0.0;
    vec[2] =  0.0;
    cubeNormal[1] = vec;

    vec[0] =  0.0;
    vec[1] =  1.0;
    vec[2] =  0.0;
    cubeNormal[2] = vec;

    vec[0] =   0.0;
    vec[1] =  -1.0;
    vec[2] =   0.0;
    cubeNormal[3] = vec;

    vec[0] =   0.0;
    vec[1] =   0.0;
    vec[2] =   1.0;
    cubeNormal[4] = vec;

    vec[0] =    0.0;
    vec[1] =    0.0;
    vec[2] =   -1.0;
    cubeNormal[5] = vec;

    size_t numfaces = amesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++)
        reorient( amesh->getFaceAt(i) );
}

int JPolyCubes:: buildPatches( const JMeshPtr &amesh)
{
    size_t numfaces = amesh->getSize(2);
    if( numfaces == 0) {
        cout << "Warning: there are no boundary faces in the mesh " << endl;
        return 1;
    }

    int side = -1, side1, side2;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = amesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("CubeSide", side);
            if( side < 0 || side > 5 ) {
                cout << "Warning: face is not assigned cube side " << endl;
                return 1;
            }
        }
    }

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = amesh->getFaceAt(i);
        if( face->isActive() ) {
            int err =  face->getAttribute("CubeSide",  side);
            if( !err ) face->setAttribute("Partition", side);
        }
    }

    JMeshPartitioner mp;
    mp.setMesh(amesh);
    mp.searchComponents();

    mp.searchInterfaces();
    mp.searchCorners();
    mp.getCorners( singularNodes );
    if( singularNodes.empty() ) {
        cout << "Warning: Singular nodes not identified: Mesh will not move " << endl;
        return 4;
    }
    mp.searchRegions();
    numInterfaces = mp.getNumInterfaces();
    numPatches    = mp.getNumPartitions();

    size_t numedges  = amesh->getSize(1);
    if( numedges == 0) {
        cout << "Warning: There are no interface edges on the boundary " << endl;
        return 2;
    }

    JFaceSequence boundfaces;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = amesh->getEdgeAt(i);
        if( edge->isActive() ) {
            if( edge->hasAttribute("Interface") ) {
                JEdge::getRelations(edge, boundfaces);
                int numneighs = boundfaces.size();
                if( numneighs != 2) {
                    cout << "Fatal Error: Every boundary edge must have two boundary faces" << endl;
                    return 1;
                }
                boundfaces[0]->getAttribute("CubeSide", side1);
                boundfaces[1]->getAttribute("CubeSide", side2);
                int side = JHexahedron::getCommonEdgeID(side1, side2);
                if( side < 0 || side > 11) {
                    cout << "Fatal Error: Incorrect edge side detected" << endl;
                    return 1;
                }
                edge->setAttribute("CubeEdge", side);
            }
        }
    }

    /*
        int numSegments = mp.getNumInterfaces();
        JEdgeSequence newsegment;
        JNodeSequence segmentnodes;

            int gid;
            JNodePtr vtx, start_node;
            for( int i = 0; i < numSegments; i++) {
                mp.getInterface(i, newsegment);
                int nSize = newsegment.size();

                start_node = nullptr;
                for( int j = 0; j < nSize; j++) {
                    JNodePtr v0 =   newsegment[j]->getNodeAt(0);
                    if( v0->hasAttribute("PartitionCorner") ) {
                        start_node = v0;
                        break;
                    }
                    JNodePtr v1 =   newsegment[j]->getNodeAt(1);
                    if( v1->hasAttribute("PartitionCorner") ) {
                        start_node = v1;
                        break;
                    }
                }

                if( start_node == nullptr)  {
                    cout << "Warning: Starting singular node not found " << endl;
                    return 1;
                }

                JEdge::getChain( newsegment, start_node);
                const JEdgePtr &edge1  = newsegment.front();
                vtx = edge1->getNodeAt(0);
                if( !vtx->hasAttribute("PartitionCorner") )  {
                    cout << "Error: Edge segment does not begin at a corner " << endl;
                    return 2;
                }

                const JEdgePtr &edge2  = newsegment.back();
                vtx = edge2->getNodeAt(1);
                if( !vtx->hasAttribute("PartitionCorner") ) {
                    cout << "Error: Edge segment does not end at a corner " << endl;
                    return 2;
                }

                if( JEdge::isClosed(newsegment) ) {
                    cout << "Error: A segment must be an open curve " << endl;
                    return 2;
                }

                interfaces.push_back(newsegment);
                targetPositions(newsegment);
                for( size_t j = 0; j < newsegment.size(); j++) {
                    const JNodePtr &v0 = newsegment[j]->getNodeAt(0);
                    const JNodePtr &v1 = newsegment[j]->getNodeAt(1);
                    gid = v0->getID();
                    v0->setAttribute("Constraint", gid);
                    gid = v1->getID();
                    v1->setAttribute("Constraint", gid);
                }
            }
    */

    return 0;
}

bool JPolyCubes ::  isIntegerMapped( const JMeshPtr &mesh) const
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

    if( integer_mapped ) {
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &anode = mesh->getNodeAt(i);
            Point3D xyz = anode->getXYZCoords();
            xyz[0] -= xmin;
            xyz[1] -= ymin;
            xyz[2] -= zmin;
            anode->setXYZCoords(xyz);
        }
    }

    return integer_mapped;
}

///////////////////////////////////////////////////////////////////////////////

int JPolyCubes :: setPolycubes(const JMeshPtr &pmesh)
{
    if( pmesh == nullptr) return 1;

    if( insurfmesh) {
        if( !insurfmesh->getTopology()->isSameAs(pmesh) ) {
            cout << "Warning: Polycube mesh has different topology than input mesh" << endl;
            return 2;
        }
    }

    pmesh->getGeometry()->setFacesNormal();

    seedClusters(pmesh);

    int err = buildPatches(pmesh);
//    if( err) return 3;

    polysurfmesh = pmesh;
    if( insurfmesh == nullptr ) return 0;

    size_t numfaces = pmesh->getSize(2);

    // Copy the sides of polycubes to the model ...
    int side, pid;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &aface = pmesh->getFaceAt(i);
        const JFacePtr &bface = insurfmesh->getFaceAt(i);
        err = aface->getAttribute("CubeSide", side);
        assert( !err);
        bface->setAttribute("CubeSide", side);

        err = aface->getAttribute("Partition", pid);
        assert(!err);
        bface->setAttribute("Partition", pid);
    }

    boost::any anyval;
    size_t numnodes = pmesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &anode = pmesh->getNodeAt(i);
        if( anode->hasAttribute("PartitionCorner") ) {
            const JNodePtr &bnode = insurfmesh->getNodeAt(i);
            bnode->setAttribute("PartitionCorner", anyval);
//            const Point3D &xyz = bnode->getXYZCoords();
        }
    }

    cout << "Is Polycube Integer Mapped " << isIntegerMapped(pmesh) << endl;

    if( deformDirection == 1) {
        sourcemesh = insurfmesh;
        targetmesh = polysurfmesh;
    } else {
        sourcemesh = polysurfmesh;
        targetmesh = insurfmesh;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////

int JPolyCubes::buildCubicalTopology()
{
    /*
        if( surfmesh == nullptr) return 1;

        if( valid_topology == 1 ) return 0;

        size_t numfaces = surfmesh->getSize(2);
        if( numfaces == 0) {
            cout << "Warning: there are no boundary faces in the mesh " << endl;
            return 1;
        }

        int side, side1, side2;
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = surfmesh->getFaceAt(i);
            if( face->isActive() ) {
                int err = face->getAttribute("CubeSide", side);
                if( side < 0 || side > 5 ) {
                    cout << "Warning: face is not assigned cube side " << endl;
                    return 1;
                }
            }
        }

        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = surfmesh->getFaceAt(i);
            if( face->isActive() ) {
                int err =  face->getAttribute("CubeSide",  side);
                if( !err ) face->setAttribute("Partition", side);
            }
        }

        JMeshPartitioner mp;
        mp.setMesh(surfmesh);

        mp.searchInterfaces();
        mp.searchCorners();
        mp.getCorners( singularNodes );
        if( singularNodes.empty() ) {
            cout << "Warning: Singular nodes not identified: Mesh will not move " << endl;
            return 4;
        }
        mp.searchRegions();

        size_t numedges  = surfmesh->getSize(1);
        if( numedges == 0) {
            cout << "Warning: There are no interface edges on the boundary " << endl;
            return 2;
        }

        JFaceSequence boundfaces;
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = surfmesh->getEdgeAt(i);
            if( edge->isActive() ) {
                if( edge->hasAttribute("Interface") ) {
                    JEdge::getRelations(edge, boundfaces);
                    int numneighs = boundfaces.size();
                    if( numneighs != 2) {
                        cout << "Fatal Error: Every boundary edge must have two boundary faces" << endl;
                        return 1;
                    }
                    boundfaces[0]->getAttribute("CubeSide", side1);
                    boundfaces[1]->getAttribute("CubeSide", side2);
                    int side = Hexahedron::getCommonEdgeID(side1, side2);
                    if( side < 0 || side > 11) {
                        cout << "Fatal Error: Incorrect edge side detected" << endl;
                        return 1;
                    }
                    edge->setAttribute("CubeEdge", side);
                }
            }
        }

        // We will reset the values after calculating the target positions...
        vector<double> vCoords;
        vector<size_t> l2g;
        surfmesh->getGeometry()->getCoordsArray( vCoords, l2g);

        int numsegments = mp.getNumInterfaces();
        JEdgeSequence newsegment;
        JNodeSequence segmentnodes;

        int gid;
        JNodePtr vtx, start_node;
        for( int i = 0; i < numsegments; i++) {
            mp.getInterface(i, newsegment);
            int nSize = newsegment.size();

            start_node = nullptr;
            for( int j = 0; j < nSize; j++) {
                 JNodePtr v0 =   newsegment[j]->getNodeAt(0);
                 if( v0->hasAttribute("PartitionCorner") ) {
                     start_node = v0;
                     break;
                 }
                 JNodePtr v1 =   newsegment[j]->getNodeAt(1);
                 if( v1->hasAttribute("PartitionCorner") ) {
                     start_node = v1;
                     break;
                 }
             }

             if( start_node == nullptr)  {
                 cout << "Warning: Starting singular node not found " << endl;
                 return 1;
             }

            JEdge::getChain( newsegment, start_node);
            const JEdgePtr &edge1  = newsegment.front();
            vtx = edge1->getNodeAt(0);
            if( !vtx->hasAttribute("PartitionCorner") )  {
                cout << "Error: Edge sement does not begin at a corner " << endl;
                return 2;
            }

            const JEdgePtr &edge2  = newsegment.back();
            vtx = edge2->getNodeAt(1);
            if( !vtx->hasAttribute("PartitionCorner") ) {
                cout << "Error: Edge sement does not end at a corner " << endl;
                return 2;
            }

            if( JEdge::isClosed(newsegment) ) {
                cout << "Error: A segment must be an open curve " << endl;
                return 2;
            }

            interfaces.push_back(newsegment);
            targetPositions(newsegment);
            for( size_t j = 0; j < newsegment.size(); j++) {
                 const JNodePtr &v0 = newsegment[j]->getNodeAt(0);
                 const JNodePtr &v1 = newsegment[j]->getNodeAt(1);
                 gid = v0->getID();
                 v0->setAttribute("Constraint", gid);
                 gid = v1->getID();
                 v1->setAttribute("Constraint", gid);
             }
        }

        JLaplaceMeshSmoother  laplace;
        laplace.setMesh(surfmesh);
        laplace.setNumIterations(100);
        laplace.smoothAll();

        size_t numnodes = surfmesh->getSize(0);
        for( size_t i = 0; i < numnodes; i++) {
             const JNodePtr &vtx = surfmesh->getNodeAt(i);
             const Point3D  &xyz = vtx->getXYZCoords();
             vtx->setAttribute("TargetPos", xyz);
        }

        surfmesh->getGeometry()->setCoordsArray( vCoords, l2g);
        surfmesh->delete_node_attribute("Constraint");
        valid_topology = 1;
    */
    return 0;
}


void JPolyCubes :: initialSegmentation()
{
    if( insurfmesh == nullptr ) return;

    seedClusters(insurfmesh);
    buildCubicalTopology();
}

///////////////////////////////////////////////////////////////////////////
void JPolyCubes :: regionGrow()
{
    if( insurfmesh == nullptr) return;

    insurfmesh->buildRelations(1,2);

    deque<JFacePtr> faceQ;
    JFaceSequence fneighs;

    // First search all the faces, which have been assigned CubeSide
    // ond of their neighbors is undecided ...

    int ival, jval, side;
    size_t numfaces = insurfmesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = insurfmesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("CubeSide", ival);
            if( ival >= 0 && ival < 6) {
                JFace::getRelations12(face, fneighs);
                for( JFacePtr fneigh: fneighs) {
                    fneigh->getAttribute("CubeSide", jval);
                    if( jval == 6)  {
                        faceQ.push_back(face);
                        break;
                    }
                }
            }
        }
    }

    while( !faceQ.empty() ) {
        JFacePtr currface = faceQ.front();
        faceQ.pop_front();
        currface->getAttribute("CubeSide",  side);
        assert( side >= 0 && side < 6) ;
        JFace::getRelations12(currface, fneighs);
        for( JFacePtr fneigh: fneighs) {
            fneigh->getAttribute("CubeSide", jval);
            if( jval == 6)  {
                fneigh->setAttribute("CubeSide", side);
                faceQ.push_back(fneigh);
            }
        }
    }

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = insurfmesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("CubeSide",  side);
            if( side < 0 || side > 5) {
                cout << "Error: One of the face side is still undecided " << endl;
                return;
            }
        }
    }

    buildCubicalTopology();
}

///////////////////////////////////////////////////////////////////////////

void JPolyCubes :: integerMap( const vector<double> &dvalues, map<double,int> &imap) const
{
    int n = dvalues.size();
    vector<double>  dv(n);
    vector<int>     di(n);

    for(int i = 0; i < n; i++)
        dv[i] = dvalues[i];

    sort( dv.begin(), dv.end() );

    vector<double> reldist(n-1);
    for( int i = 0; i < n-1; i++)
        reldist[i] = max(1.0E-10, dv[i] - dv[i-1]);

    di[0] = floor( dv[0] );
    for( int i = 1; i < n; i++)
        di[i] = di[i-1] + ceil(reldist[i]);

    cout << "Integer Snap " << endl;
    for( int i = 0; i < n; i++) {
        imap[dv[i]] = di[i];
        cout << dv[i] << "  " << di[i] << endl;
    }
}

///////////////////////////////////////////////////////////////////////////

int JPolyCubes :: targetPositions( const JEdgeSequence &edges)
{

    if( edges.empty() ) return 1;

    JNodeSequence nodes;
    JEdgeTopology::getChainNodes(edges, nodes);

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

int JPolyCubes :: alignAlongXYZPlanes()
{
    if( insurfmesh == nullptr ) return 1;

    if( !valid_topology)  return 1;

    JEdgeSequence edges;
    JFaceSequence faces;

    if( tetmesh == nullptr) {
        AllTetMeshGenerator tetmesher;
        tetmesh = tetmesher.getConstrainedMesh(sourcemesh);
    }

    for( const JNodePtr &vtx : singularNodes) {
        int id = vtx->getID();
        vtx->setAttribute("Constraint", id);
    }

    JLocallyInjectiveMap limDeform;
    limDeform.setMesh(tetmesh);
    limDeform.solve();
    return 0;
}

///////////////////////////////////////////////////////////////////////////
int JPolyCubes :: flatten( const JFaceSequence &patch)
{
    /*
        double xsum = 0.0;
        double ysum = 0.0;
        double zsum = 0.0;
    */
    return 0;
}
///////////////////////////////////////////////////////////////////////////

JMeshPtr JPolyCubes :: getPolycubeMesh1()
{
    /*
        if( insurfmesh == nullptr ) return nullptr;

        vector<double> vCoords;
        vector<size_t> l2g;
        insurfmesh->getGeometry()->getCoordsArray( vCoords, l2g);

        JMeshPartitioner mp;
        mp.setMesh(insurfmesh);

        JMeshPtr target1 = insurfmesh->deepCopy();

        insurfmesh->getGeometry()->setCoordsArray( vCoords, l2g);
        return target1;
    */
    return nullptr;
}

///////////////////////////////////////////////////////////////////////////

JMeshPtr JPolyCubes :: integerSnap()
{
    if( polysurfmesh == nullptr) {
        cout << "Warning: A tight polycubes must exist at this stage " << endl;
        return nullptr;
    }

    JNodeSequence nodes;
    size_t numnodes = polysurfmesh->getSize(0);

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = polysurfmesh->getNodeAt(i);
        if( vtx->hasAttribute("PartitionCorner") )
            nodes.push_back(vtx);
    }

    if( nodes.empty() ) {
        cout << "Warning: Singular points not found in the polycubes " << endl;
        return nullptr;
    }

    std::map<double, JNodeSequence> xmap, ymap, zmap;

    for (size_t i = 0; i < nodes.size(); i++) {
        const Point3D &xyz = nodes[i]->getXYZCoords();
        xmap[xyz[0]].push_back(nodes[i]);
        ymap[xyz[1]].push_back(nodes[i]);
        zmap[xyz[2]].push_back(nodes[i]);
        nodes[i]->setAttribute("TargetPos", xyz);
    }

    int nsize;

    vector<double> xval, yval, zval;
    std::map<double,int> ival, jval, kval;

    nsize = xmap.size();
    assert( nsize);
    for( auto keyVal : xmap)  xval.push_back(keyVal.first);
    integerMap(xval, ival);

    nsize = ymap.size();
    assert( nsize);
    for( auto keyVal : ymap)  yval.push_back(keyVal.first);
    integerMap(yval, jval);

    nsize = zmap.size();
    assert( nsize);
    for( auto keyVal : zmap)  zval.push_back(keyVal.first);
    integerMap(zval, kval);

    Point3D xyz;
    for( auto keyVal : xmap)  {
        double x =  keyVal.first;
        const JNodeSequence &xnodes = keyVal.second;
        for( const JNodePtr &vtx : xnodes) {
            vtx->getAttribute("TargetPos", xyz);
            xyz[0] = ival[x];
            vtx->setAttribute("TargetPos", xyz);
        }
    }

    for( auto keyVal : ymap)  {
        double y =  keyVal.first;
        const JNodeSequence &ynodes = keyVal.second;
        for( const JNodePtr &vtx : ynodes) {
            vtx->getAttribute("TargetPos", xyz);
            xyz[1] = jval[y];
            vtx->setAttribute("TargetPos", xyz);
        }
    }

    for( auto keyVal : zmap)  {
        double z =  keyVal.first;
        const JNodeSequence &znodes = keyVal.second;
        for( const JNodePtr &vtx : znodes) {
            vtx->getAttribute("TargetPos", xyz);
            xyz[2] = kval[z];
            vtx->setAttribute("TargetPos", xyz);
        }
    }

    JLocallyInjectiveMap limDeform;
    limDeform.setMesh(polysurfmesh);
    limDeform.solve();
    return polysurfmesh;
}

///////////////////////////////////////////////////////////////////////////
int  JPolyCubes :: deformSource( int dir)
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
    return 0;
}

///////////////////////////////////////////////////////////////////////////
void JPolyCubes :: getIntegerPoints( const JCellPtr &tet, vector<Point3D> &iPoints)
{
    iPoints.clear();

    double xmin =  0.99*std::numeric_limits<double>::max();
    double ymin =  0.99*std::numeric_limits<double>::max();
    double zmin =  0.99*std::numeric_limits<double>::max();
    double xmax = -xmin;
    double ymax = -ymin;
    double zmax = -zmin;

    for( int it = 0; it < 4; it++) {
        const Point3D &xyz = tet->getNodeAt(it)->getXYZCoords();
        xmin = min( xmin, xyz[0] );
        ymin = min( ymin, xyz[1] );
        zmin = min( zmin, xyz[2] );

        xmax = max( xmax, xyz[0] );
        ymax = max( ymax, xyz[1] );
        zmax = max( zmax, xyz[2] );
    }

    int i0 = floor(xmin)-1;
    int j0 = floor(ymin)-1;
    int k0 = floor(zmin)-1;

    int i1 = ceil(xmax)+1;
    int j1 = ceil(ymax)+1;
    int k1 = ceil(zmax)+1;

    Point3D xyz;
    for( int k = k0; k < k1; k++) {
        for( int j = j0; j < j1; j++) {
            for( int i = i0; i < i1; i++) {
                xyz[0] = i;
                xyz[1] = j;
                xyz[2] = k;
                if( TetGeometry::isInside(tet, xyz) )
                    iPoints.push_back(xyz);
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////

JMeshPtr JPolyCubes :: getPolyHexMesh()
{
    if( polysurfmesh == nullptr) return nullptr;

    double xmin =  0.99*std::numeric_limits<double>::max();
    double ymin =  0.99*std::numeric_limits<double>::max();
    double zmin =  0.99*std::numeric_limits<double>::max();
    double xmax = -xmin;
    double ymax = -ymin;
    double zmax = -zmin;

    int integer_mapped = 1;
    size_t numnodes = polysurfmesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &anode = polysurfmesh->getNodeAt(i);
        if( anode->hasAttribute("PartitionCorner") ) {
            Point3D &xyz = anode->getXYZCoords();
            if( fabs((int)xyz[0] -xyz[0]) > 1.0E-10 ) integer_mapped = 0;
            if( fabs((int)xyz[1] -xyz[1]) > 1.0E-10 ) integer_mapped = 0;
            if( fabs((int)xyz[2] -xyz[2]) > 1.0E-10 ) integer_mapped = 0;
            xmin = min( xmin, xyz[0] );
            ymin = min( ymin, xyz[1] );
            zmin = min( zmin, xyz[2] );

            xmax = max( xmax, xyz[0] );
            ymax = max( ymax, xyz[1] );
            zmax = max( zmax, xyz[2] );
        }
    }
    assert( integer_mapped);

    double length[3];
    length[0] = xmax-xmin;
    length[1] = ymax-ymin;
    length[2] = zmax-zmin;

    double origin[3];
    origin[0] = xmin;
    origin[1] = ymin;
    origin[2] = zmin;

    int nodeDim[3];
    nodeDim[0] = xmax - xmin + 1;
    nodeDim[1] = ymax - ymin + 1;
    nodeDim[2] = zmax - zmin + 1;

    JMeshPtr hmesh  = AllHexMeshGenerator::getStructuredMesh(nodeDim, length, origin);
    if( tetmesh == nullptr) {
        AllTetMeshGenerator tetmesher;
        tetmesh = tetmesher.getConstrainedMesh(polysurfmesh);
    }
    JMeshIO::saveAs(tetmesh, "tet.xml");

    size_t nCount = tetmesh->getGeometry()->getNumOfInvertedElements();
    assert( nCount == 0);

    size_t numNodes = tetmesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = tetmesh->getNodeAt(i);
        vtx->setVisitBit(0);
    }

    size_t numCells = tetmesh->getSize(3);

    size_t numHexNodes = nodeDim[0]*nodeDim[1]*nodeDim[2];
    cout << "#hex nodes " << numHexNodes << endl;

    vector<Point3D> iPoints;
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = tetmesh->getCellAt(i);
        getIntegerPoints( cell, iPoints);
        for( size_t j = 0; j < iPoints.size(); j++) {
            int x = iPoints[j][0];
            int y = iPoints[j][1];
            int z = iPoints[j][2];
            size_t id = z*nodeDim[0]*nodeDim[1] + y*nodeDim[0] + x;
            if( id < numHexNodes) {
                const JNodePtr &vtx = hmesh->getNodeAt(id);
                vtx->setVisitBit(1);
            }
        }
    }

    numCells = hmesh->getSize(3);
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = hmesh->getCellAt(i);
        bool visited = 0;
        for( int j = 0; j < 8; j++) {
            const JNodePtr &vtx = cell->getNodeAt(j);
            if( vtx->getVisitBit() == 1) {
                visited = 1;
                break;
            }
        }
        cell->setVisitBit(visited);
    }

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = hmesh->getCellAt(i);
        if( cell->getVisitBit() == 0)
            cell->setStatus( JMeshEntity::REMOVE);
    }

    JNodeSet nodeSet;
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = hmesh->getCellAt(i);
        if( cell->isActive() ) {
            for( int j = 0; j < 8; j++)
                nodeSet.insert( cell->getNodeAt(j) );
        }
    }

    JNodeSequence activenodes;
    map<JNodePtr, JNodePtr> mapnodes;
    size_t index = 0;
    for( const JNodePtr &oldnode : nodeSet) {
        JNodePtr newnode = JNode::newObject();
        Point3D  xyz = oldnode->getXYZCoords();
        newnode->setXYZCoords(xyz);
        newnode->setID( index++);
        mapnodes[oldnode] = newnode;
        activenodes.push_back( newnode);
    }

    JCellSequence activeCells;
    JNodeSequence connect(8);

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = hmesh->getCellAt(i);
        if( cell->isActive() ) {
            for( int j = 0; j < 8; j++)
                connect[j] = mapnodes[cell->getNodeAt(j)];
            JHexahedronPtr h = JHexahedron::newObject();
            h->setNodes( connect);
            activeCells.push_back(h);
        }
    }

    JMeshPtr amesh = JMesh::newObject();
    amesh->addObjects( activenodes);
    amesh->addObjects( activeCells);
    return amesh;

}

///////////////////////////////////////////////////////////////////////////
void JPolyCubes :: saveForMatlab()
{
    if( sourcemesh == nullptr) {
        cout << "Warning: source surface mesh not provided " << endl;
        return;
    }

    if( targetmesh == nullptr) {
        cout << "Warning: target surface mesh not provided " << endl;
        return;
    }
    assert( sourcemesh->getSize(0) == targetmesh->getSize(0) );

    if( tetmesh == nullptr) {
        AllTetMeshGenerator tetmesher;
        tetmesh = tetmesher.getConstrainedMesh( sourcemesh );
        JMeshIO::saveAs(tetmesh, "tri.xml");
    }

    ofstream ofile;
    ofile.open("X_source.dat", ios::out);

    size_t numNodes = tetmesh->getSize(0);
    size_t numBound = sourcemesh->getSize(0);
    for( size_t i = 0; i < numBound; i++) {
        const JNodePtr &vtx = sourcemesh->getNodeAt(i);
        assert( vtx->getID() == i );
        const Point3D  &xyz = vtx->getXYZCoords();
        ofile << xyz[0] << "  " << xyz[1] << "  " << xyz[2] << endl;
    }
    for(size_t i = numBound; i < numNodes; i++) {
        const JNodePtr &vtx = tetmesh->getNodeAt(i);
        assert( vtx->getID() == i );
        const Point3D  &xyz = vtx->getXYZCoords();
        ofile << xyz[0] << "  " << xyz[1] << "  " << xyz[2] << endl;
    }
    ofile.close();

    ofile.open("X_target.dat", ios::out);
    for( size_t i = 0; i < numBound; i++) {
        const JNodePtr &vtx = targetmesh->getNodeAt(i);
        assert( vtx->getID() == i );
        const Point3D  &xyz = vtx->getXYZCoords();
        ofile << xyz[0] << "  " << xyz[1] << "  " << xyz[2] << endl;
    }

    for(size_t i = numBound; i < numNodes; i++) {
        const JNodePtr &vtx = tetmesh->getNodeAt(i);
        assert( vtx->getID() == i );
        const Point3D  &xyz = vtx->getXYZCoords();
        ofile << xyz[0] << "  " << xyz[1] << "  " << xyz[2] << endl;
    }
    ofile.close();

    ofile << "# name: tri" << endl;
    ofile << "# type: matrix "  << endl;
    ofile << "# rows: " << tetmesh->getSize(3) << endl;
    ofile << "# columns: 4" << endl;

    ofile.open("tri.dat", ios::out);
    size_t numCells = tetmesh->getSize(3);
    for(size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell= tetmesh->getCellAt(i);
        ofile << cell->getNodeAt(0)->getID() + 1 << "  "
              << cell->getNodeAt(1)->getID() + 1 << "  "
              << cell->getNodeAt(2)->getID() + 1 << "  "
              << cell->getNodeAt(3)->getID() + 1 << endl;
    }
    ofile.close();

    ofile.open("constraint_matrix.dat", ios::out);
    for( size_t i = 0; i < numBound; i++) {
        const JNodePtr &vtx = sourcemesh->getNodeAt(i);
        assert( vtx->getID() == i);
//        int colid = 3*vtx->getID() + 1;
        ofile << 3*i+1 << "  " <<  3*i+1 << " 1 " << endl;
        ofile << 3*i+2 << "  " <<  3*i+2 << " 1 " << endl;
        ofile << 3*i+3 << "  " <<  3*i+3 << " 1 " << endl;
    }
    ofile.close();

    ofile.open("constraint_rhs.dat", ios::out);
    for( size_t i = 0; i < numBound; i++) {
        const JNodePtr &vtx = sourcemesh->getNodeAt(i);
        const Point3D  &xyz = vtx->getXYZCoords();
        ofile << xyz[0] << endl;
        ofile << xyz[1] << endl;
        ofile << xyz[2] << endl;
    }
    ofile.close();
}
///////////////////////////////////////////////////////////////////////////
