#include "MSTQuadMesher.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesher :: setMesh( const JMeshPtr &m)
{
    mesh.reset();
    if( m->getSize(2) == 0) {
        JDelaunayMesh2D delmesher;
        delmesher.setMesh(m);
        mesh = delmesher.getSimpleMesh();
    } else
        mesh = m;

    if( mesh == nullptr) return;

    mesh->deleteNodeAttribute("Interface");
    mesh->deleteEdgeAttribute("Interface");

    singularMesh.reset();
    quadmesh.reset();
    mesh->buildRelations(0,2);

    // Set Target edges length based on boundary edges ...
    mesh->getTopology()->getBoundary(domainBoundaryNodes);

    planarMesh = 1;
    if( mesh->getTopology()->isClosed() ) planarMesh = 0;

    double elen = mesh->getGeometry()->getMeanEdgeLength();
    mesh->setAttribute("TargetEdgeLength", elen);
}

////////////////////////////////////////////////////////////////////////////////////

JMeshPtr JMSTQuadMesher :: getTemplate( const vector<int> &segments)
{
    //************************************************************************
    // This module will generate a quad mesh in a canonical shape. Total
    // number of segments must be even..
    //************************************************************************
    int nsize = segments.size();
    if( nsize < 2 || nsize > 6) {
        cout << "Info: Current polygon size " << nsize << endl;
        cout << "Warning: Quad template valid for 2,3,4,5, and 6 sides only" << endl;
        return nullptr;
    }

    int nsum = 0;
    for( int i = 0; i < nsize; i++) {
        nsum += segments[i];
        if( segments[i] < 1 ) return nullptr;
    }


    if( nsum < 4 ) {
        cout << "Warning: Too few number of segments" << endl;
        return nullptr;
    }

    if( nsum%2 != 0 ) {
        cout << "Warning: number of segments must be even" << endl;
        for( int i : segments) cout << i << endl;
        return nullptr;
    }

    //**************************************************************************
    // For the time being: we use System command to execute the program ..
    //**************************************************************************
    ofstream ofile("patch.dat", ios::out);
    if( ofile.fail() ) return nullptr;

    ofile << nsize << endl;
    for( int i = 0; i < nsize; i++)
        ofile << segments[i] << " ";
    ofile.close();

    //**************************************************************************
    // The algorithm is guaranteed to produce quad elements, when the number of
    // segments are even, and no new boundary vertices are introduced ...
    // Whether we accept this mesh or not is decided later ...
    //**************************************************************************

    system("mstquadpatch patch.dat patch.off");

    // Start the patch from lower left corner and make sure that the boundary
    // nodes are numbered first starting from 0 in the anti-clockwise direction.

    JMeshPtr patch = JMeshIO::readFile("patch.off");
    if( patch == nullptr) return nullptr;

    patch->getTopology()->searchBoundary();

    double start_angle = 0.0;

    switch( nsize) {
    case 3:
        start_angle = 210;
        break;
    case 4:
        start_angle = 225;
        break;
    case 5:
        start_angle = 234;
        break;
    case 6:
        start_angle = 240;
        break;
    }

    JNodeSequence boundnodes;
    patch->getTopology()->getBoundary(boundnodes);

    JNodePtr start_node;
    for( const JNodePtr &vtx : boundnodes) {
        const Point2D &xy = vtx->getXYCoords();
        double angle  = atan2(xy[1], xy[0] );
        if( angle < 0.0) angle += 2.0*M_PI;
        angle *= 180.0/M_PI;
        if( fabs(angle-start_angle) < 1.0E-03) {
            start_node = vtx;
            break;
        }
    }
    assert( start_node );

    JEdgeSequence boundedges;
    patch->getTopology()->getBoundary(boundedges);
    JEdgeTopology::getChain( boundedges, start_node);
    JEdgeTopology::getChainNodes( boundedges, boundnodes);

    size_t index = 0;
    for( const JNodePtr &vtx : boundnodes)
        vtx->setID(index++);

    size_t numnodes = patch->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = patch->getNodeAt(i);
        if( !vtx->isBoundary() ) vtx->setID( index++);
    }

    JEntityCompare  *jcomp = new JEntityIDCompare;
    patch->sort(jcomp);
    delete jcomp;

    return patch;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMSTQuadMesher :: getPatch( const JNodeSequence &srcnodes, const vector<int> &segments)
{
    // Get the quad mesh. It will always give a quad mesh.
    JMeshPtr meshtemplate = getTemplate(segments);
    if( meshtemplate == nullptr) return nullptr;

    // Now map the patch to the physical patch ...
    UV2XYZMap(meshtemplate, srcnodes, segments);

    return meshtemplate;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMSTQuadMesher :: getLoftPatch( const JNodeSequence &side1,
        const JNodeSequence &side2,
        const JNodeSequence &side3,
        const JNodeSequence &side4)
{
    assert( side1.front() == side4.front() );
    assert( side1.back()  == side2.front() );
    assert( side3.back()  == side2.back()  );
    assert( side4.back()  == side3.front() );

    vector<int> segments(4);
    segments[0] = side1.size()-1;
    segments[1] = side2.size()-1;
    segments[2] = side3.size()-1;
    segments[3] = side4.size()-1;

    int nsum = segments[0] + segments[1] + segments[2] + segments[3];

    JNodeSequence srcnodes(nsum);

    int index = 0;

    int np = side1.size();
    for ( int i = 0; i < np-1; i++)
        srcnodes[index++] = side1[i];

    np = side2.size();
    for ( int i = 0; i < np-1; i++)
        srcnodes[index++] = side2[i];

    np = side3.size();
    for ( int i = 0; i < np-1; i++)
        srcnodes[index++] = side3[np-1-i];

    np = side4.size();
    for ( int i = 0; i < np-1; i++)
        srcnodes[index++] = side4[np-1-i];

    JMeshPtr loftmesh = getPatch(srcnodes, segments);
    return loftmesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMSTQuadMesher :: getAdaptivePatch( const JNodeSequence &srcNodes,
        const vector<int> &segments,
        double factor)
{
    int  numBound = srcNodes.size();
    int  numSides = segments.size();
    if ( numBound < 3 || numSides < 3 || numSides > 6 || numBound%2 ) return nullptr;

    int index = 0;
    // get the outer patch ...
    JMeshPtr meshtemplate1 = getTemplate(segments);
    JMeshAffineTransform affine;
    affine.setMesh( meshtemplate1 );
    affine.scale(1.5, 1.5, 1.0);
    JNodeSequence outerNodes(numBound);
    for( int i = 0; i < numBound; i++)
        outerNodes[i] = meshtemplate1->getNodeAt(i);

    vector<JNodeSequence> outerSideNodes(numSides);
    int nsize = outerNodes.size();
    index = 0;
    for( int i = 0; i < numSides; i++) {
        for( int j = 0; j <= segments[i]; j++) {
            outerSideNodes[i].push_back(outerNodes[index%numBound]);
            index++;
        }
        index--;
    }
    meshtemplate1->deleteFaces();

    // Create new segments and keep the parity same ...
    int numNewBound = 0;
    vector<int>  newSegments(numSides);
    for( int i = 0; i < numSides; i++) {
        int nedges = factor*segments[i];
        int nsides = max(1, nedges);
        if( nsides%2 != segments[i]%2) nsides++;
        newSegments[i] = nsides;
        numNewBound += nsides;
    }

    // get the innermost patch ...
    JMeshPtr meshtemplate2 = getTemplate(newSegments);
    vector<JNodeSequence > innerSideNodes(numSides);
    index = 0;
    for( int i = 0; i < numSides; i++) {
        for( int j = 0; j <= newSegments[i]; j++) {
            const JNodePtr &vtx = meshtemplate2->getNodeAt(index%numNewBound);
            innerSideNodes[i].push_back(vtx);
            index++;
        }
        index--;
    }
    JFaceSequence innerFaces = meshtemplate2->getFaces();
    int id = 1;
    for( const JFacePtr &face : innerFaces)
        face->setAttribute("Partition", id);

//  Genrate bridge : The edges which connected outer to inner nodes.
//  int nbridge = max(1.0, fabs(1.0-factor)*srcNodes.size()/8.0);
    int nbridge = 1;

    vector<JNodeSequence> bridgeNodes(numSides);
    for( int i = 0; i < numSides; i++) {
        JEdgeGeometry::generateLinearNodes( outerSideNodes[i][0], innerSideNodes[i][0], nbridge+1,
                                            bridgeNodes[i] );
        for( size_t j = 1; j < bridgeNodes[i].size()-1; j++)
            meshtemplate2->addObject( bridgeNodes[i][j] );
    }

    meshtemplate2->addObjects( outerNodes );
    JNodeSequence innerNodes;
    for( int i = 0; i < numSides; i++) {
        JMeshPtr padmesh = getLoftPatch( outerSideNodes[i], bridgeNodes[(i+1)%numSides],
                                         innerSideNodes[i], bridgeNodes[i]);
        innerFaces = padmesh->getFaces();
        padmesh->getTopology()->getInternal( innerNodes);
        meshtemplate2->addObjects( innerNodes );
        meshtemplate2->addObjects( innerFaces );
        id = i + 2;
        for( const JFacePtr &face:innerFaces)
            face->setAttribute("Partition", id);

    }
    meshtemplate2->enumerate(0);

    meshtemplate2->getTopology()->collectEdges();
    meshtemplate2->getTopology()->searchBoundary();

    JLloydMeshOptimizer mopt;
    mopt.setMesh(meshtemplate2);
    mopt.setNumIterations(500);
    mopt.smoothAll();

    index = 0;
    for( const JNodePtr &vtx : outerNodes) {
        vtx->setID(index++);
        assert( vtx->isBoundary() );
    }

    innerFaces = meshtemplate2->getFaces();
    meshtemplate2->getTopology()->getSubmeshInternal(innerFaces, innerNodes);
    for( const JNodePtr &vtx : innerNodes) {
        vtx->setID(index++);
        assert( !vtx->isBoundary() );
    }

    JEntityCompare  *jcomp = new JEntityIDCompare;
    meshtemplate2->sort(jcomp);
    delete jcomp;

    assert(  outerNodes.size() == srcNodes.size() );
    UV2XYZMap(meshtemplate2, srcNodes, segments);

    return meshtemplate2;
}

///////////////////////////////////////////////////////////////////////////////

void  JMSTQuadMesher :: UV2XYZMap(const JMeshPtr &meshtemplate,
                                      const JNodeSequence &srcnodes,
                                      const vector<int> &segments)
{
    JNodeSequence dstnodes = meshtemplate->getNodes();
    int numbounds = srcnodes.size();

    int nsegs = segments.size();
    int numCagePoints = max(3, nsegs);

    // Do the affine mapping between parametric to the physical space ...
    JMeanValueCoordinates mvc;

    vector<Point2D>  cageCoords(numCagePoints);
    vector<int>      cageIndex( numCagePoints);

    int pos = 0;
    if( nsegs < 3) {
        for( int i = 0; i < numCagePoints; i++) {
            cageIndex[i] = pos;
            pos += numbounds/3;
        }
    } else {
        for( int i = 0; i < numCagePoints; i++) {
            cageIndex[i] = pos;
            pos += segments[i];
        }
    }

    for( int i = 0; i < numCagePoints; i++) {
        pos = cageIndex[i];
        const Point3D &xyz = dstnodes[pos]->getXYZCoords();
        cageCoords[i][0] = xyz[0];
        cageCoords[i][1] = xyz[1];
    }
    mvc.setBaryCoords( meshtemplate, cageCoords);

    patchCorners.resize( numCagePoints );
    for( int i = 0; i < numCagePoints; i++) {
        pos = cageIndex[i];
        const Point3D &xyz = srcnodes[pos]->getXYZCoords();
        patchCorners[i]  = srcnodes[pos];
        cageCoords[i][0] = xyz[0];
        cageCoords[i][1] = xyz[1];
    }
    mvc.setXYCoords( meshtemplate, cageCoords);

    int grp = 0;
    for( int i = 0; i < numbounds; i++) {
        Point3D p0 = dstnodes[i]->getXYZCoords();
        dstnodes[i]->setAttribute("Constraint", grp);
        assert( dstnodes[i]->isBoundary() );
        Point3D p1 = srcnodes[i]->getXYZCoords();
        dstnodes[i]->setXYZCoords(p1);
    }

    if( optimize ) {
        JLloydMeshOptimizer mopt;
        mopt.setMesh(meshtemplate);
        mopt.setNumIterations(100);
        mopt.smoothAll();
    }

    stitchBoundary( meshtemplate, srcnodes);
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMSTQuadMesher :: getTemplate( const JMeshPtr &m, int nside)
{
    mesh = m;
    if( mesh == nullptr ) return nullptr;

    // Collect the boundary edges of the model ( for which we want a quad mesh)
    // and collect all the nodes in anti-clockwise direction. For a quad mesh
    // to exist, the number of boundary edges must be even ...

    JEdgeSequence edges;
    mesh->getTopology()->getBoundary( edges);

    size_t numedges = edges.size();
    if( numedges%2 ) return nullptr;
    JEdgeTopology::getChain(edges);

    JNodeSequence srcnodes;
    JEdgeTopology::getChainNodes( edges, srcnodes);

    // Now decide how the edges will be split into "nsides"

    vector<int> segments(nside);
    for( int i = 0; i < nside; i++)
        segments[i] =  numedges/nside;
    segments[nside-1] += numedges%nside;

    // Get the quad mesh. It will always give a quad mesh.
    JMeshPtr meshtemplate = getPatch( srcnodes, segments);
    if( meshtemplate == nullptr) return nullptr;

    /*
        // Now we deform the mesh to match with the model boundary...
        size_t numnodes = srcnodes.size();
        int grp = 0;
        for( size_t i = 0; i < numnodes; i++) {
            const Point3D &xyz = srcnodes[i]->getXYZCoords();
            dstnodes[i]->setXYZCoords(xyz);
            dstnodes[i]->setAttribute("Constraint", grp);
        }
        JLaplaceMeshSmoother mlap;
        mlap.setMesh( meshtemplate);
        mlap.setNumIterations(100);
        mlap.smoothAll();
        mlap.untangle();

        JMeshIO::saveAs( meshtemplate, "tri.off");

        size_t nCount = meshtemplate->getGeometry()->getNumInvertedElements();
        if( nCount) {
            JMeshPtr tmesh = AllTriMeshGenerator::getFromQuadMesh(meshtemplate,2);
            JMeshNonlinearOptimization mopt;
            mopt.setMesh(tmesh);
    //      mopt.untangle();
            nCount = tmesh->getGeometry()->getNumInvertedElements();
            if( nCount == 0) mopt.improveShapes();
        }
    */
    return  meshtemplate;

}

////////////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesher :: stitchBoundary( const JMeshPtr &meshtemplate,
        const JNodeSequence &chainnodes)
{
    if( meshtemplate == nullptr) return 1;

    JEdgeSequence boundedges; // These will be removed and will be replaced by mesh edges.
    JNodeSequence boundnodes; // These will be removed and will be replaced by mesh nodes.
    meshtemplate->getTopology()->getBoundary(boundnodes);
    meshtemplate->getTopology()->getBoundary(boundedges);

    assert( boundnodes.size() == chainnodes.size() );
    assert( boundedges.size() == chainnodes.size() );

    JNodeSequence nodes = meshtemplate->getNodes();
    size_t numbound  = chainnodes.size();

    // replace the node of the template ...
    for( size_t i = 0; i < numbound; i++) {
        assert( nodes[i]->isBoundary() ) ;
        nodes[i] = chainnodes[i];
        nodes[i]->setID(i);
    }

    meshtemplate->addObjects(chainnodes);

    // Create the new faces of the template ...
    size_t numfaces = meshtemplate->getSize(2);
    JNodeSequence newConnect(4);
    int id;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = meshtemplate->getFaceAt(i);
        for( int j = 0; j < 4; j++)
            newConnect[j] = nodes[ face->getNodeAt(j)->getID() ];
        JFacePtr newface = JQuadrilateral::newObject( newConnect);
        int err = face->getAttribute("Partition", id);
        if( !err) newface->setAttribute("Partition", id);
        meshtemplate->addObject( newface );
    }

    // We can delete the faces of the mesh templates now..
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = meshtemplate->getFaceAt(i);
        face->setStatus( JMeshEntity::REMOVE);
    }

    for( const JEdgePtr &edge : boundedges)
        edge->setStatus( JMeshEntity::REMOVE);

    for( const JNodePtr &node : boundnodes)
        node->setStatus( JMeshEntity::REMOVE);

    meshtemplate->getTopology()->collectEdges();

    meshtemplate->pruneAll();
    meshtemplate->enumerate(0);
    meshtemplate->enumerate(1);
    meshtemplate->enumerate(2);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesher :: nsidedQuads( const JFacePtr &face)
{
    int nside = face->getSize(0);

    JNodeSequence chainnodes, edgenodes;
    vector<int>  segments(nside);

    for( int i = 0; i < nside; i++) {
        const JNodePtr &node = face->getNodeAt(i);
        const JEdgePtr &edge = face->getEdgeAt(i);
        int err = edge->getAttribute("Steiner", edgenodes);
        assert(!err);
        chainnodes.push_back(node);
        int nnodes = edgenodes.size();
        segments[i] = nnodes+1;
        if( face->getOrientation(edge) > 0) {
            for( int j = 0; j < nnodes; j++)
                chainnodes.push_back( edgenodes[j] );
        } else {
            for( int j = 0; j < nnodes; j++)
                chainnodes.push_back( edgenodes[nnodes-j-1] );
        }
    }
    JMeshPtr meshtemplate = getPatch( chainnodes, segments);
    if( meshtemplate == nullptr) return 1;

    JNodeSequence nodes = meshtemplate->getNodes();
    size_t numbound  = chainnodes.size();
    for( size_t i = 0; i < numbound; i++)
        nodes[i] = chainnodes[i];

    size_t numnodes = meshtemplate->getSize(0);
    for( size_t i = numbound; i < numnodes; i++) {
        const JNodePtr &vtx = meshtemplate->getNodeAt(i);
        quadmesh->addObject( vtx);
    }
    size_t numfaces = meshtemplate->getSize(2);
    JNodeSequence newConnect(4);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = meshtemplate->getFaceAt(i);
        for( int j = 0; j < 4; j++)
            newConnect[j] = nodes[ face->getNodeAt(j)->getID() ];
        JFacePtr newface = JQuadrilateral::newObject( newConnect);
        quadmesh->addObject( newface );
    }

    face->setStatus(JMeshEntity::REMOVE);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesher :: getConvexPatches()
{
    // Decompose the geometric model convex elements. Original mesh
    // will be deleted ....
    JPolyPartitioner polypart;
    polypart.setMesh( singularMesh);
    singularMesh = polypart.getPartitions();

    // All the polygons must have 3,4,5,6 nodes. If not split the face ....
    int numfaces = singularMesh->getSize(2);
    for( int i = 0; i < numfaces; i++) {
        const JFacePtr &face = singularMesh->getFaceAt(i);
        if( face->getSize(0) > 6 ) {
            splitPolygon(face);
        }
    }
    singularMesh->pruneFaces();
}

/////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesher :: genEdgeNodes()
{
    size_t numedges = mesh->getSize(1);

    double sum = 0.0;
    for( size_t i = 0;  i < numedges; i++) {
        const JEdgePtr &edge =  mesh->getEdgeAt(i);
        sum +=  JNodeGeometry::getLength(edge->getNodeAt(0), edge->getNodeAt(1) );
    }
    double maxlen = sum/numedges;

    // Divided each edge into number of "even" segments....
    JNodeSequence newnodes;

    numedges = singularMesh->getSize(1);
    int numSegments;
    for( size_t i = 0;  i < numedges; i++) {
        const JEdgePtr &edge =  singularMesh->getEdgeAt(i);
        if( !edge->hasAttribute("Steiner") ) {
            double len = JNodeGeometry::getLength(edge->getNodeAt(0), edge->getNodeAt(1) );
            numSegments = len/maxlen;
            numSegments  = max(2, numSegments);
            if( numSegments%2 ) numSegments++;
            JEdgeGeometry::generateLinearNodes( edge, numSegments+1, newnodes);
            edge->setAttribute("Steiner", newnodes);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
int JMSTQuadMesher :: make_even_segments(JEdgeSequence &edgeloop, const JEdgePtr &edge)
{
    int numSegments = edgeloop.size();

    const JNodePtr &vstart = edge->getNodeAt(0);
    const JNodePtr &vend   = edge->getNodeAt(1);

    int pos0 = -1;
    for( size_t i = 0; i < numSegments; i++) {
        if( edgeloop[i]->getNodeAt(0) == vstart) {
            pos0 = i;
            break;
        }
    }

    int pos1 = -1;
    for( int j = pos0+1; j < numSegments+1; j++) {
        if( edgeloop[j%numSegments]->getNodeAt(0) == vend) {
            pos1 = j;
            break;
        }
    }

    JNodeSequence steiner;
    steiner.clear();
    for( int j = pos0+1; j < pos1; j++)  {
        const JNodePtr &vtx = edgeloop[j%numSegments]->getNodeAt(0);
        steiner.push_back(vtx);
    }

    JNodePtr vnew;
    if( steiner.empty() ) {
        vnew = JNodeGeometry::getMidNode( vstart, vend);
        steiner.push_back(vnew);
    }

    if( steiner.size()%2 ) {
        edge->setAttribute("Steiner", steiner);
        return 0;
    }

    vnew = JNodeGeometry::getMidNode( vstart, steiner[0]);
    steiner.insert( steiner.begin(), vnew);

    edge->setAttribute("Steiner", steiner);

    return 0;
}
/////////////////////////////////////////////////////////////////////////////

QDefectivePatchPtr JMSTQuadMesher :: getDefectivePatch( const JNodePtr &vtx)
{
    /*
        if( geodesic_path ) {
            mesh->buildRelations(0,0);
            JMeshGeodesics djk;
            djk.setMesh(mesh, nullptr);
            JNodePtr vnear;
            JNodeSequence nodes = getInteriorSingularNodes();
            double mindist = numeric_limits<double>::max();
            for( const JNodePtr &v : nodes) {
                if( v != vtx ) {
                double dist = JNodeGeometry::getLength(vtx,v);
                if( dist < mindist) {
                    mindist = dist;
                    vnear   = v;
                }
                }
            }
            JNodeSequence nodepath = djk.getPath(vtx,vnear);
            return getDefectivePatch(nodepath);
        }
    */

    QDefectivePatchPtr  patch = boost::shared_ptr<QDefectivePatch>(new QDefectivePatch);
    patch->setMesh( mesh );
    patch->setSeed( vtx );
    patch->setMinSingularNodes(minSingularNodes);
    patch->build();
    return patch;
}

/////////////////////////////////////////////////////////////////////////////

QDefectivePatchPtr JMSTQuadMesher :: getDefectivePatch( const JFacePtr &face)
{
    /*
        if( geodesic_path ) {
            mesh->buildRelations(0,0);
            DijkstraShortestPath  djk;
            djk.setMesh(mesh);
            JNodePtr vnear;
            JNodeSequence nodes = getInteriorSingularNodes();
            double mindist = numeric_limits<double>::max();
            for( const JNodePtr &v : nodes) {
                if( v != vtx ) {
                double dist = JNodeGeometry::getLength(vtx,v);
                if( dist < mindist) {
                    mindist = dist;
                    vnear   = v;
                }
                }
            }
            JNodeSequence nodepath = djk.getPath(vtx,vnear);
            return getDefectivePatch(nodepath);
        }
    */
    QDefectivePatchPtr  patch = boost::shared_ptr<QDefectivePatch>(new QDefectivePatch);

    patch->setMesh( mesh );
    patch->setSeed( face );
    patch->build();
    return patch;
}

/////////////////////////////////////////////////////////////////////////////////


QDefectivePatchPtr JMSTQuadMesher :: getDefectivePatch( const JNodeSequence &nodepath)
{
    QDefectivePatchPtr  patch = boost::shared_ptr<QDefectivePatch>(new QDefectivePatch);

    patch->setMesh( mesh );
    patch->setInitialPath( nodepath );
    patch->setMinSingularNodes(minSingularNodes);

    patch->build();
    return patch;
}
/////////////////////////////////////////////////////////////////////////////
void JMSTQuadMesher :: reorderSingularities()
{
    if( singularNodesQ.empty() ) return;

    if(singularity_removal_policy) {
        boost::sort(singularNodesQ, []( const JNodePtr &v1, const JNodePtr &v2)
        {
            double dist1 = 0.0;
            v1->getAttribute("Distance", dist1);
            double dist2 = 0.0;
            v2->getAttribute("Distance", dist2);
            return dist1 < dist2;
        });

        if( singularity_removal_policy == INTERIOR_REMOVAL)
            boost::reverse(singularNodesQ);
    }
}
/////////////////////////////////////////////////////////////////////////////

QDefectivePatchPtr JMSTQuadMesher :: getAnyDefectivePatch()
{
    if(select_global_patch) return getGlobalPatch();

    if( singularNodesQ.empty() ) {
        JNodeSequence nodes = getSingularNodes();
        for( const JNodePtr &v : nodes) singularNodesQ.push_back(v);

        JNodeSequence boundnodes;
        mesh->getTopology()->getBoundary(boundnodes);
        JGraphGeodesics geodesic;
        geodesic.setMesh(mesh);
        geodesic.setDistanceField(boundnodes);
        reorderSingularities();
    }

    if( singularNodesQ.empty() ) return nullptr;

    JNodePtr apex;
    while( !singularNodesQ.empty() ) {
        apex = singularNodesQ.front();
        singularNodesQ.pop_front();
        if( apex->isActive() ) break;
    }

    QDefectivePatchPtr pch = getDefectivePatch(apex);
    return pch;
}
/////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesher :: remesh( const QDefectivePatchPtr  &patch)
{
    // A patch must be a topological disk ...
    if( patch == nullptr)  return 1;
    if( !patch->isValid()) return 1;

    // Collect all the faces, edges and nodes that we will remove from the mesh...
    JFaceSequence oldfaces = patch->getFaces();
    JNodeSequence oldInternalNodes;
    JEdgeSequence oldInternalEdges;
    mesh->getTopology()->getSubmeshInternal(oldfaces, oldInternalEdges, oldInternalNodes);

    // Get the boudary edges of the patch and ensure that number of segments are even.
    // Evenness is necessary condition for quadrangulation....
    JEdgeSequence boundedges = patch->getBoundEdges();
    JNodeSequence  boundnodes;
    JEdgeTopology::getChainNodes(boundedges, boundnodes);

    if( boundnodes.size()%2 ) {
        cout << "Warning: number of segments must be even for quad meshing " << endl;
        return 1;
    }
    if( boundnodes.size() < 3) {
        cout << "Warning: Minimum three segments must be provided for quad meshing " << endl;
        return 2;
    }

    // Create the new MST patch now ...
    JEdgeGeometry::smooth(boundedges,5);
    vector<int> segments  = patch->getSegments();
    JMeshPtr meshtemplate = getPatch(boundnodes, segments);

    if( meshtemplate == nullptr) {
        patch->clear();
        return 2;
    }

    // Collect new internal edges and  nodes from the patch..
    // Remove entities of the patch...(all faces, internal edges, and internal nodes);
    for( const JFacePtr &f : oldfaces)
        f->setStatus( JMeshEntity::REMOVE);

    for( const JEdgePtr &e : oldInternalEdges)
        e->setStatus( JMeshEntity::REMOVE);

    for( const JNodePtr &v : oldInternalNodes)
        v->setStatus( JMeshEntity::REMOVE);
    mesh->pruneAll();

    // Now after removal of patch, there is a hole in the mesh. We will fill it
    // with the new patch. Note that the next stage is valid only when there is
    // a hole in the mesh..

    JFaceSequence newInternalFaces = meshtemplate->getFaces();
    JNodeSequence newInternalNodes;
    JEdgeSequence newInternalEdges;
    meshtemplate->getTopology()->getSubmeshInternal(newInternalFaces, newInternalEdges, newInternalNodes);

    if( meshtemplate->getSize(2) < 1 || newInternalEdges.empty() ) {
        cout << "Fatal error: The new mesh template is invalid" << endl;
        cout << "Writing template in debug.off" << endl;
        cout << meshtemplate->getSize(0) << " " << meshtemplate->getSize(1) << " "
             << meshtemplate->getSize(2) << endl;
        meshtemplate->enumerate(0);
        JMeshIO::saveAs(meshtemplate, "debug.off");
        exit(0);
    }

    mesh->addObjects(newInternalNodes);
    mesh->addObjects(newInternalEdges);
    mesh->addObjects(newInternalFaces);

    // Smooth the entire mesh ...
    JLloydMeshOptimizer  mopt;
    mopt.setMesh(mesh);
    mopt.setNumIterations(10);
    mopt.smoothAll();

    // Update the position of irregular nodes from the boundary (using Euclidean distance)
    meshtemplate->enumerate(0);
    meshtemplate->buildRelations(0,2);
    size_t numnodes = meshtemplate->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = meshtemplate->getNodeAt(i);
        if( vtx->getNumRelations(2) != 4) updateDistanceMap(vtx);
    }

    //Store the template inside the patch object ...
    patch->setTemplateMesh(meshtemplate);

    // Enumerate objects ...
    mesh->enumerate(2);
    mesh->enumerate(1);
    mesh->enumerate(0);

    // Ready to make new patch....
    patch->clear();

    return 0;
}

/////////////////////////////////////////////////////////////////////////////
void  JMSTQuadMesher :: straighten( const JEdgePtr &edge)
{
    JNodeSequence nodes;
    int err = edge->getAttribute("Steiner", nodes);
    if( err ) return;

    const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
    const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();

    int numnodes = nodes.size();
    double dt = 1.0/(double)(numnodes+1);

    Point3D xyz;
    for( int i = 0; i < numnodes; i++) {
        double t = (i+1)*dt;
        xyz[0] = (1-t)*p0[0] + t*p1[0];
        xyz[1] = (1-t)*p0[1] + t*p1[1];
        xyz[2] = (1-t)*p0[2] + t*p1[2];
        nodes[i]->setXYZCoords(xyz);
    }
}

/////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesher :: getBoundarySingularGraph( JEdgeSequence &edgeloop)
{
    if( !JEdgeTopology::isChain( edgeloop ))
        JEdgeTopology::getChain(edgeloop);

    for( const JEdgePtr &edge: edgeloop) {
        const JNodePtr &vtx = edge->getNodeAt(0);
        if( vtx->hasAttribute("SingularNode") ) {
            JEdgeTopology::getChain(edgeloop, vtx);
            cout << "Search for defective node " << endl;
            break;
        }
    }

    JNodeSequence nodes;
    JEdgeTopology::getChainNodes(edgeloop, nodes);

    int numCorners = 0;
    int nodeGroup  = 0;
    for( const JNodePtr &vtx : nodes) {
        if( vtx->hasAttribute("CornerAngle")) {
            vtx->setAttribute("SingularNode", nodeGroup);
            numCorners++;
        }
    }

    /*
        JMeshContour mc;
        mc.setSegments(edgeloop);

        JNodeSequence corners;
        if( numCorners == 0 ) {
            corners.resize(4);
            corners[0] = mc.getNearestNode( 0.0);
            corners[1] = mc.getNearestNode( 0.25);
            corners[2] = mc.getNearestNode( 0.50);
            corners[3] = mc.getNearestNode( 0.75);
            for( const JNodePtr &vtx : corners)
                vtx->setAttribute("SingularNode", nodeGroup);
        }

        if( numCorners == 1 ) {
            corners.resize(4);
                JEdge::getChain(edgeloop, corners[0]);
            corners[1] = mc.getNearestNode( 0.25);
            corners[2] = mc.getNearestNode( 0.50);
            corners[3] = mc.getNearestNode( 0.75);
            for( const JNodePtr &vtx : corners)
                vtx->setAttribute("SingularNode", nodeGroup);
        }

        if( numCorners == 2 ) {
            corners.resize(4);
            JEdge::getChain(edgeloop, corners[0]);
            corners[2] = corners[1];
            corners[1] = mc.getNearestNode( 0.25);
            corners[3] = mc.getNearestNode( 0.75);
            for( const JNodePtr &vtx : corners)
                vtx->setAttribute("SingularNode", nodeGroup);
        }
    */

    JNodeSequence snodes;

    for( size_t i = 0; i < nodes.size(); i++) {
        const JNodePtr &vtx  = nodes[i];
        if( nodes[i]->hasAttribute("SingularNode") ) {
            snodes.push_back(vtx);
        }
    }
    singularMesh->addObjects( snodes);

    JNodeSequence innodes;
    int nSize = snodes.size();
    for( int i = 0; i < nSize; i++) {
        const JNodePtr &v0 = snodes[i];
        const JNodePtr &v1 = snodes[(i+1)%nSize];
        JEdgePtr edge = JSimplex::getEdgeOf(v0, v1, 1);
        singularMesh->addObject(edge);
        make_even_segments(edgeloop, edge);
        straighten(edge);
    }
}

/////////////////////////////////////////////////////////////////////////////
JMeshPtr JMSTQuadMesher :: getBoundarySingularGraph()
{
    if( mesh == nullptr) return nullptr;
    mesh->getTopology()->searchBoundary();

    mesh->deleteNodeAttribute("TargetPos");
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isBoundary() ) {
            const Point3D &xyz = vtx->getXYZCoords();
            vtx->setAttribute("TargetPos", xyz);
        }
    }

    segmentBoundary();

    /*
        rectDomain = 1;
        double angle;
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            int err = vtx->getAttribute("CornerAngle", angle);
            if( !err) {
                if( angle < 100  && fabs(angle-90) > 1.0) {
                    rectDomain = 0;
                    break;
                }
                if( angle > 260  && fabs(angle-270) > 1.0) {
                    rectDomain = 0;
                    break;
                }
            }
        }
        if( rectDomain) return nullptr;
    */

    vector<JEdgeSequence> edgeloops;
    mesh->getTopology()->getBoundary( edgeloops);

    singularMesh = JMesh::newObject();

    // Go to each edgeloop and complete the singular graph...
    int numloops = edgeloops.size();
    for( int i = 0; i < numloops; i++)
        getBoundarySingularGraph( edgeloops[i] );

    // From the singular graph, construct convex partitions of the
    // model...
    getConvexPatches();

    return singularMesh;
}

/////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesher :: segmentBoundary()
{
    // Find out how many holes in the mesh ....
    vector<JEdgeSequence> edgeloops;
    mesh->getTopology()->getBoundary(edgeloops);

    double sum = 0;
    size_t numSegments = 0;
    for( size_t i = 0; i < edgeloops.size(); i++) {
        for( const JEdgePtr &edge : edgeloops[i] ) {
            const JNodePtr &v0 = edge->getNodeAt(0);
            const JNodePtr &v1 = edge->getNodeAt(1);
            sum += JNodeGeometry::getLength(v0,v1);
            numSegments++;
        }
    }
    avg_edgelen = sum/(double)numSegments;

    // Search for the corners on the boundary.Each corner is considered a
    // singular nodes. In noisy model, we should smooth the boundary before
    // this step to reduced the initial singular nodes ...
    JNodeSequence corners = mesh->getGeometry()->getBoundaryCorners(cornerAngle);

    int val = 0;
    for( const JNodePtr &vtx : corners)
        vtx->setAttribute("SingularNode", val);
    return 0;
}

/////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesher :: getVoxelMesh()
{
    size_t numnodes = mesh->getSize(0);

    vector<double> xCorner, yCorner;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->hasAttribute("SingularNode") ) {
            const Point3D &xyz = vtx->getXYZCoords();
            xCorner.push_back( xyz[0] );
            yCorner.push_back( xyz[1] );
        }
    }

    boost::sort( xCorner );
    boost::sort( yCorner );

    std::unique( xCorner.begin(), xCorner.end()  );
    std::unique( yCorner.begin(), yCorner.end()  );

    int nSize = xCorner.size();
    vector<double> xRelativeDist(nSize);
    xRelativeDist[0] = 0.0;
    for( int i = 1; i < nSize; i++)
        xRelativeDist[i] = xCorner[i] - xCorner[i-1];

    vector<int> xInteger(nSize);
    xInteger[0] = floor( xCorner[0] );

    for( int i = 1; i < nSize; i++)
        xInteger[i] = ceil(xInteger[i-1] + xRelativeDist[i]);

    nSize = yCorner.size();
    vector<double> yRelativeDist(nSize);
    yRelativeDist[0] = 0.0;
    for( int i = 1; i < nSize; i++)
        yRelativeDist[i] = yCorner[i] - yCorner[i-1];

    vector<int> yInteger(nSize);
    xInteger[0] = floor( yCorner[0] );

    for( int i = 1; i < nSize; i++)
        yInteger[i] = ceil(yInteger[i-1] + yRelativeDist[i]);

    map<double,int> xmap, ymap;
    for( int i = 0; i < nSize; i++) {
        xmap[xCorner[i]] = xInteger[i];
        ymap[yCorner[i]] = yInteger[i];
    }

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->hasAttribute("SingularNode") ) {
            Point3D xyz = vtx->getXYZCoords();
            xyz[0]  = xmap[ xyz[0]];
            xyz[1]  = ymap[ xyz[1]];
        }
    }

    /*
        int dim[2];
        dim[0] = xInteger[nSize] - xInteger[0] + 1;
        dim[1] = yInteger[nSize] - yInteger[0] + 1;

        JDelaunayMesh2D trimesher;

        JEdgeSequence boundedges;
        mesh->getTopology()->getBoundary(boundedges);
    //  quadmesh = trimesher.addSegments( boundedges);
    */

}

/////////////////////////////////////////////////////////////////////////////

/*
JNodeSequence JMSTQuadMesher :: getInteriorSingularNodes()
{
    JNodeSequence nodes;
    if( mesh == nullptr ) return nodes;
    mesh->buildRelations(0,2);

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() && !vtx->isBoundary() )  {
            if( vtx->getNumRelations(2) != 4) nodes.push_back(vtx);
        }
    }
    return nodes;
}
*/

/////////////////////////////////////////////////////////////////////////////
JNodeSequence JMSTQuadMesher :: getSingularNodes()
{
    JNodeSequence nodes;
    if( mesh == nullptr ) return nodes;
    mesh->buildRelations(0,2);

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() )  {
            if( vtx->getNumRelations(2) != 4) {
                if( planarMesh) {
                    if( !vtx->isBoundary() )
                        nodes.push_back(vtx);
                } else {
                    nodes.push_back(vtx);
                }
            }
        }
    }
    return nodes;
}
/////////////////////////////////////////////////////////////////////////////

JMeshPtr JMSTQuadMesher :: getQuadMesh()
{
    if( mesh == nullptr) return nullptr;

    cout << "Info: Constructing Singular graph " << endl;
    getBoundarySingularGraph();

    cout << "Info: Singular graph stored in file : singular.off" << endl;
    JMeshIO::saveAs( singularMesh, "singular.off");

    cout << "Info: Refine edges   " << endl;
    // Discretize all the edges ....
    genEdgeNodes();

    cout << "Info: applying Quad templates  " << endl;
    quadmesh  = JMesh::newObject();

    JNodeSequence nodes = singularMesh->getNodes();
    quadmesh->addObjects( nodes);

    int numedges = singularMesh->getSize(1);

    for( int i = 0; i < numedges; i++) {
        const JEdgePtr &edge = singularMesh->getEdgeAt(i);
        int err = edge->getAttribute("Steiner", nodes);
        assert(!err);
        quadmesh->addObjects( nodes);
    }

#ifdef DEBUG
    cout << "Info: Convex decomposition stored in file : convexdecomp.off" << endl;
    JMeshIO::saveAs( quadmesh, "convexdecomp.off");
#endif

    // Apply template on each face ....
    int numfaces = singularMesh->getSize(2);
    for( int i = 0; i < numfaces; i++) {
        const JFacePtr &face = singularMesh->getFaceAt(i);
        if( face->isActive() ) {
            nsidedQuads(face);
        }
    }
    cout << "Quadmeshing done ... " << endl;

    return quadmesh;
}

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void JMSTQuadMesher :: pillowPatch( const JMeshPtr &patch, const vector<int> &segments, int minQuads)
{
    JQuadRefiner  refiner;
    refiner.setMesh(patch);

    int sum = 0;
    for( const int val : segments) sum += val;

    int numPillows =   (minQuads-patch->getSize(2))/sum;

    if( numPillows < 1) return;

    for( int i = 0; i < numPillows; i++)
        refiner.insert_boundary_pillows();

    JLloydMeshOptimizer mopt;
    mopt.setMesh(patch);
    mopt.setNumIterations(100);
    mopt.smoothAll();
}
///////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesher :: applyLIM( const JMeshPtr &quadmesh)
{
    JMeshUntangle untangle;
    untangle.setMesh(quadmesh);
    untangle.execute();

    size_t numfaces = quadmesh->getSize(2);
    for( int i = 0; i < numfaces; i++) {
        const JFacePtr &face = quadmesh->getFaceAt(i);
        if( JFaceGeometry::isInverted(face) )  {
            cout << "Warning: Locally injective mapping failed to untangle " << endl;
            cout << "Info: lim mesh stored in lim.off for debugging" << endl;
            JMeshIO::saveAs(quadmesh, "lim.off");
            return 1;
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesher :: updateDistanceMap( const JNodePtr &vi)
{
    const Point3D &pi = vi->getXYZCoords();
    double mindist = numeric_limits<double>::max();
    for( const JNodePtr &vj: domainBoundaryNodes) {
        const Point3D &pj = vj->getXYZCoords();
        double len = JMath::length(pi,pj);
        mindist  = min(mindist,len);
    }
    distMap.insert(make_pair(mindist,vi));
}
/////////////////////////////////////////////////////////////////////////////

QDefectivePatchPtr JMSTQuadMesher :: getGlobalPatch()
{
    if( distMap.empty() ) {
        JNodeSequence nodes = getSingularNodes();
        if( nodes.empty() ) return nullptr;
        mesh->getTopology()->getBoundary(domainBoundaryNodes);
        for( const JNodePtr &vi: nodes) updateDistanceMap( vi ) ;
    }

    if( distMap.empty() ) return nullptr;

    /*
        while( !distMap.empty() ) {
            JNodePtr vtx = singularNodesQ.front();
            singularNodesQ.pop_front();
            if( vtx->isActive() ) return buildCircularPatch(vtx);
        }
    */
    return nullptr;

}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesher :: splitPolygon( const JFacePtr &oldface)
{
    int nnodes = oldface->getSize(0);

    if( nnodes < 7) return 0;

    JFacePtr face1, face2;

    JNodeSequence nodes;
    for( int i = 0; i < 6; i++)
        nodes.push_back(oldface->getNodeAt(i));

    face1= JPolygon::newObject( nodes);

    nodes.clear();
    nodes.push_back( oldface->getNodeAt(0) );
    for( int i = 5; i < nnodes; i++)
        nodes.push_back( oldface->getNodeAt(i));

    switch( nodes.size() )
    {
    case 3:
        face2 = JTriangle::newObject( nodes);
        break;
    case 4:
        face2 = JQuadrilateral::newObject( nodes);
        break;
    default:
        face2 = JPolygon::newObject( nodes);
        break;
    }
    singularMesh->addObject( face1 );
    singularMesh->addObject( face2 );

    oldface->setStatus(JMeshEntity::REMOVE);
    splitPolygon(face2);
    return 0;
}


JMeshPtr JMSTQuadMesher :: getQuadTemplate( const JFacePtr &face, double edgelen)
{
    /*
       if( face->getSize(0) != 4) return nullptr;

       JQuadRefiner quadRefiner;
       JNodeSequence newnodes;
       JFaceSequence newfaces;
       quadRefiner.refine5( face, newnodes, newfaces);

       const JNodePtr &v0 = face->getNodeAt(0);
       const JNodePtr &v1 = face->getNodeAt(1);

       double l1   = JNodeGeometry::getLength(v0,v1);
       int    n1   = l1/edgelen;
       if( n1%2 == 0) n1++;

       double l2   = JNodeGeometry::getLength(v0,newnodes[0]);
       int    n2   = max(1, (int)l1/edgelen);
    */
}

QDefectivePatchPtr JMSTQuadMesher :: buildCircularPatch(const JNodePtr &vtx)
{
    /*
        QDefectivePatchPtr  patch = boost::shared_ptr<QDefectivePatch>(new QDefectivePatch);

            if( !vtx->isActive()) {
                cout << "Warning: Inactive seed " << endl;
                return nullptr;
            }
            patch->setMesh( mesh );

            double mindist = numeric_limits<double>::max();

            double dist;
            Point3D  pCenter = vtx->getXYZCoords();

            size_t numnodes = mesh->getSize(0);
            for( size_t i = 0; i < numnodes; i++) {
                const JNodePtr &v = mesh->getNodeAt(i);
                if( v->isActive() ) {
                    const Point3D &pj = v->getXYZCoords();
                    dist =  JMath::length(pCenter, pj);
                    if( v->isBoundary() )
                        mindist = min( mindist, dist);
                    v->setAttribute("Distance", dist);
                }
            }
            double radius = mindist;

            JFaceSequence patchfaces;

            size_t numfaces  = mesh->getSize(2);
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &f = mesh->getFaceAt(i);
                if( f->isActive() ) {
                    int nn = f->getSize(0);
                    bool inside = 1;
                    for( int j = 0; j < nn; j++) {
                        const JNodePtr &v = f->getNodeAt(j);
                        v->getAttribute("Distance", dist);
                        if( dist > radius || v->isBoundary() ) {
                            inside = 0;
                            break;
                        }
                    }
                    if( inside ) patchfaces.push_back(f);
                }
            }

            if( patchfaces.empty() )  return  nullptr;
            if( patchfaces.size() < 4) return nullptr;

            patch->setCorners(4);
            patch->build(patchfaces);
            return patch;
    */
}

