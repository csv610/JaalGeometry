#include "Mesh.hpp"

#include "QuadCleanUp.hpp"

using namespace Jaal;

//////////////////////////////////////////////////////////////////////////////

void QDefectivePatch :: setMaximumFaces( int m ) {
    maxSearchedFaces = m;
    /*
        if( m ) stopAtBoundary = 0;
        if( m == 0) {
            stopAtBoundary = 1;
            maxSearchedFaces = std::numeric_limits<int>::max();
        }
    */
}

//////////////////////////////////////////////////////////////////////////////

void QDefectivePatch:: getTemplate()
{
    meshTemplate.reset();

    //************************************************************************
    // This module will generate a quad mesh in a canonical shape. Total
    // number of segments must be even..
    //************************************************************************
    int nsize = segSize.size();
    if( nsize < 2 || nsize > 6) {
        cout << "Info: Current polygon size " << nsize << endl;
        cout << "Warning: Quad template valid for 2,3,4,5, and 6 sides only" << endl;
        return;
    }

    int nsum = 0;
    for( int i = 0; i < nsize; i++) {
        nsum += segSize[i];
        if( segSize[i] < 1 ) return;
    }

    if( nsum < 4 ) {
        cout << "Warning: Too few number of segments" << endl;
        return;
    }

    if( nsum%2 != 0 ) {
        cout << "Warning: number of segments must be even" << endl;
        for( int i : segSize) cout << i << endl;
        return;
    }

    //**************************************************************************
    // For the time being: we use System command to execute the program ..
    //**************************************************************************
    ofstream ofile("patch.dat", ios::out);
    if( ofile.fail() ) return;

    ofile << nsize << endl;
    for( int i = 0; i < nsize; i++)
        ofile << segSize[i] << " ";
    ofile.close();

    //**************************************************************************
    // The algorithm is guaranteed to produce quad elements, when the number of
    // segments are even, and no new boundary vertices are introduced ...
    // Whether we accept this mesh or not is decided later ...
    //**************************************************************************

    system("mstquadpatch patch.dat patch.off");

    // Start the patch from lower left corner and make sure that the boundary
    // nodes are numbered first starting from 0 in the anti-clockwise direction.

    meshTemplate = JMeshIO::readFile("patch.off");

    meshTemplate->getTopology()->searchBoundary();

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
    meshTemplate->getTopology()->getBoundary(boundnodes);

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
    meshTemplate->getTopology()->getBoundary(boundedges);
    JEdgeTopology::getChain( boundedges, start_node);
    int err = JEdgeTopology::getChainNodes( boundedges, boundnodes);
    if( err ) {
        meshTemplate.reset();
        return;
    }

    size_t index = 0;
    for( const JNodePtr &vtx : boundnodes)
        vtx->setID(index++);

    size_t numnodes = meshTemplate->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = meshTemplate->getNodeAt(i);
        if( !vtx->isBoundary() ) vtx->setID( index++);
    }

    JEntityCompare  *jcomp = new JEntityIDCompare;
    meshTemplate->sort(jcomp);
    delete jcomp;
}

///////////////////////////////////////////////////////////////////////

bool QDefectivePatch :: isReplacementProfitable()
{
    //
    // A patch is considered "Profitable" is the remesh patch
    // has fewer number of singulariies ( interior + boundary) than
    // the original submesh...
    //

    // Singularities in the submesh ....
    JNodeSet  patchNodes;
    for( const JFacePtr &face: faceSet) {
        patchNodes.insert( face->getNodeAt(0) );
        patchNodes.insert( face->getNodeAt(1) );
        patchNodes.insert( face->getNodeAt(2) );
        patchNodes.insert( face->getNodeAt(3) );
    }

    int singularitiesBefore = 0;
    for( const JNodePtr &vtx : patchNodes) {
        if( vtx->getNumRelations(2) != 4) {
            singularitiesBefore++;
        }
    }

    // Singularities after the patch replacement ....

    // Get the new template ...
    getTemplate();
    if( meshTemplate == nullptr) return 0;

    // Find the boundary of the template ...
    meshTemplate->getTopology()->searchBoundary();

    meshTemplate->buildRelations(0,2);
    size_t nnodes = meshTemplate->getSize(0);

    // Singularities in the interior  ...
    int singularitiesAfter  = 0;
    for( size_t i = 0; i < nnodes; i++) {
        const JNodePtr  &vtx = meshTemplate->getNodeAt(i);
        if( !vtx->isBoundary() && vtx->getNumRelations(2) != 4)
            singularitiesAfter++;
    }

    // Singularities at the boundary ...
    int nbound = boundNodes.size();
    for( int i = 0; i < nbound; i++) {
        const JNodePtr &vsrc = boundNodes[i];
        int nOutside   = getTopologicalAngle(vsrc) + 2;
        const JNodePtr &vdst = meshTemplate->getNodeAt(i);
        assert( vdst->isBoundary() );
        int nInside  =  vdst->getNumRelations(2);
        if( (nOutside + nInside) != 4 )  singularitiesAfter++;
    }

    return singularitiesAfter < singularitiesBefore;
}
///////////////////////////////////////////////////////////////////////
void QDefectivePatch :: clear()
{
    nodePath.clear();
    cornerNodes.clear();

    innerNodes.clear();
    boundNodes.clear();

    faceSet.clear();
    edgeSet.clear();

    boundEdges.clear();
    segSize.clear();
}

////////////////////////////////////////////////////////////////////
bool QDefectivePatch :: isValid() const
{
    if( boundEdges.empty() ) return 0;
    if( boundEdges.size()%2 ) return 0;
    if( segSize.empty() ) return 0;
    if( validPatch ) return 1;
    return 0;
}

////////////////////////////////////////////////////////////////////
double QDefectivePatch :: getIsoperimeticQuotient() const {
    // This definiation is taken from Wikipedia..
    double A = getArea();
    double L = getPerimeter();
    double q = 4*M_PI*A/(L*L);
    return q;
}

////////////////////////////////////////////////////////////////////
JFaceSequence  QDefectivePatch :: getFaces() const
{
    JFaceSequence result;
    size_t nSize = faceSet.size();

    result.clear();
    if( nSize == 0 ) return result;

    result.resize( nSize );

    int index = 0;
    JFaceSet::const_iterator it;
    for( const JFacePtr &f : faceSet)
        result[index++] = f;
    return result;
}
////////////////////////////////////////////////////////////////////

Point3D QDefectivePatch :: getCenter() const
{
    if( apexNode ) return apexNode->getXYZCoords();

    Point3D pc;
    pc[0] = 0.0;
    pc[1] = 0.0;
    pc[2] = 0.0;
    return pc;
}
////////////////////////////////////////////////////////////////////

int QDefectivePatch::initBlob()
{
    nodes.clear();
    innerNodes.clear();
    boundNodes.clear();
    faceSet.clear();
    edgeSet.clear();
    relations02.clear();

    int err;
    if (apexNode) {
        if( apexNode->isBoundary() ) return 1;
        err = expandBlob( apexNode );
        if( err ) return 1;
    } else {
        for( const JNodePtr &vtx : nodePath) {
            err = expandBlob( vtx );
            if( err ) return 1;
        }
    }

    updateBoundary();

    // If the apex has very large degree, then instead of
    // using "face-open" operation, just enforce that there
    // are four corner points on the patch and remesh it.
    // This is another example of the power of "MST" that
    // we do not need usual location operations ....
    // Note that we do not use "topological angle" concept
    // for this patch--we just pick four points in the
    // boundary ad make them corners.

    specifiedCorners = 0;
    if( apexNode) {
        if( apexNode->getNumRelations(2) > 5) specifiedCorners = 4;
    }

    finalizeBoundary();

    return 0;
}

////////////////////////////////////////////////////////////////////

size_t QDefectivePatch:: getNumSingularNodes(int where) const
{
    assert(mesh->getAdjTable(0, 2));

    size_t ncount = 0;

    if (where == 0) {
        for( const JNodePtr &v : innerNodes) {
            if( v->isActive() && (v->getNumRelations(2) != 4)) ncount++;
        }
    } else {
        for( const JNodePtr &v : boundNodes) {
            if( v->isActive() && (v->getNumRelations(2) == 4)) ncount++;
        }
    }
    return ncount;
}

////////////////////////////////////////////////////////////////////
int
QDefectivePatch::updateBoundary()
{
    boundEdges.clear();
    boundNodes.clear();

    JFaceSequence fneighs;
    int nCount = 0;
    for( const JEdgePtr &edge : edgeSet) {
        JEdge::getRelations(edge, fneighs);
        nCount = 0;
        for( const JFacePtr &f: fneighs)
            if( faceSet.find(f) != faceSet.end() ) nCount++;
        if( nCount == 1) boundEdges.push_back(edge);
    }

    if( boundEdges.empty() ) return 1;

    JEdgeTopology::getChain(boundEdges);
    JEdgeTopology::getChainNodes(boundEdges, boundNodes);

    if( boundNodes.empty() ) return 1;

    JNodeSet aset, bset;
    for( const JNodePtr &vtx : boundNodes) {
        bset.insert(vtx);
        if(vtx->isBoundary() ) {
            boundEdges.clear();
            boundNodes.clear();
            return 1;
        }
    }

    for( const JFacePtr &f : faceSet) {
        aset.insert( f->getNodeAt(0) );
        aset.insert( f->getNodeAt(1) );
        aset.insert( f->getNodeAt(2) );
        aset.insert( f->getNodeAt(3) );
    }

    innerNodes.clear();
    boost::set_difference(aset, bset, std::back_inserter(innerNodes));

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
QDefectivePatch::build()
{
    nodes.clear();
    innerNodes.clear();
    boundNodes.clear();
    faceSet.clear();
    edgeSet.clear();
    relations02.clear();

//  if( apexFace ) return buildTrianglesPatch();
    if( apexNode ) return buildSingularityPatch();

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
QDefectivePatch::buildSingularityPatch()
{
    validPatch = 0;

    assert(mesh->getAdjTable(0, 2));

    /////////////////////////////////////////////////////////////////////////////
    // A valid boundary patch has the following propeties:
    //
    // 1) The loop is closed.
    // 2) The region is simply connected. The Euler characteristic must be
    //    one, which mean there are no holes in the region.
    // 3) The loop has 3, 4 or 5 corner nodes.
    // 4) The loop form almost convex boundary.
    //
    // If these conditions are met, we may try to remesh the region with
    // quadrilateral elements. Sometimes it may not be possible to remesh
    // a reason, or the resulting elements have unacceptable quality that
    // we would like to avoid.
    //
    // A resulting mesh MUST decrease the irregular nodes count.
    //////////////////////////////////////////////////////////////////////////
    int err, topo_convex_region;

    err = initBlob();
    if( err ) return 1;

    // If you encounter very high valence node in the beginning itself, we
    // should first get rid of it....
    if( apexNode ) {
        if( apexNode->getNumRelations(2) > 6 && cornerNodes.size() > 1) {
            validPatch = 1;
            return 0;
        }
    }

    specifiedCorners = 0;

    newFaceSet.clear();

    size_t  nSize;
    while (1) {
        err = updateBoundary();
        if (err) return 2;

        nSize = boundNodes.size();
        assert( nSize > 0);

//      size_t before_expansion = faceSet.size();
        // Expand the blob till topological set is found ..
        newFaceSet.clear();

        topo_convex_region = 1;
        for (size_t i = 0; i < nSize; i++) {
            int topo_angle = getTopologicalAngle(boundNodes[i]);
            if (topo_angle < 0) {
                err = expandBlob(boundNodes[i]);
                if( err ) return 1;
                topo_convex_region = 0;
            }
        }

        // A valid patch must contain atleast two irregular nodes ...
        if( topo_convex_region ) {
            if (getNumSingularNodes(0) < (size_t)minSingularNodes) {
                topo_convex_region = 0;
                err = expandBlob();
                if( err ) return 1;
            }
        }

        // Now finalize the patch for remeshing...
        if (topo_convex_region) {
            finalizeBoundary();
            int nsides = cornerNodes.size();
            if( nsides < 1) return 1;
            if( isReplacementProfitable() ) {
                validPatch = 1;
                return 0;
            }
            cout << "Non-Profitable replacement: expand more " << endl;
            err = expandBlob();
            if( err ) return 1;

        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void QDefectivePatch::getBoundNodes(const JNodePtr &src, const JNodePtr &dst, JNodeSequence &seq)
{
    int start_pos = getPosOf(src);
    int end_pos = getPosOf(dst);
    int nsize = boundNodes.size();

    assert(nsize > 1);

    if (end_pos == 0) end_pos = nsize;
    assert(end_pos > start_pos);

    seq.resize(end_pos - start_pos + 1);
    int index = 0;
    for (int i = start_pos; i <= end_pos; i++)
        seq[index++] = boundNodes[i % nsize];
}

///////////////////////////////////////////////////////////////////////////////

bool QDefectivePatch::isSimple()
{
    size_t V = relations02.size();
    size_t F = faceSet.size();
    size_t E = edgeSet.size();
    if (V - E + F == 1) return 1;
    return 0;
}

/////////////////////////////////////////////////////////////////////////

int
QDefectivePatch::finalizeBoundary()
{
    cornerNodes.clear();
    // Sequence the chain and start from one of the corner...
    int err = JEdgeTopology::getChain(boundEdges);
    if (err) return 2;

    if( JEdgeGeometry::getOrientation(boundEdges) < 0)
        JEdgeTopology::reverse(boundEdges);

    err = JEdgeTopology::getChainNodes( boundEdges, boundNodes);
    if( err ) return 3;

    JNodeSequence cornersFrom;
    if( specifiedCorners == 0) {
        for( const JNodePtr &vtx : boundNodes) {
            if( !vtx->isBoundary() && getTopologicalAngle(vtx) > 0)
                cornersFrom.push_back(vtx);
        }
        int ncorners = cornersFrom.size();
        if( ncorners ) {
            int nselect = min(6,ncorners);
            int nskip   = cornersFrom.size()/nselect;
            cornerNodes.resize(nselect);
            for( int i = 0; i < nselect; i++)
                cornerNodes[i]  = cornersFrom[i*nskip];
        }
    } else {
        for( const JNodePtr &vtx : boundNodes)
            if( !vtx->isBoundary()) cornersFrom.push_back(vtx);
        cornerNodes.resize(specifiedCorners);
        int nskip = cornersFrom.size()/specifiedCorners;
        for( int i = 0; i < specifiedCorners; i++)
            cornerNodes[i]  = cornersFrom[i*nskip];
    }

    if( cornerNodes.size() < 1 ) return 3;

    // Start the chain from one of the cornerNodes.
    err = JEdgeTopology::getChain(boundEdges, *cornerNodes.begin());
    if (err) return 4;

    if( cornerNodes.size() > 6) cornerNodes.resize(6);

    JEdgeTopology::getChainNodes( boundEdges, boundNodes);

    // Split the boundary loop into segments.
    // (i.e. End of the segments are the cornerNodes identified earlier )
    set_boundary_segments();
    return 0;
}

////////////////////////////////////////////////////////////////////

void QDefectivePatch::set_boundary_segments()
{
    int nSize = cornerNodes.size();
    // Although this stage will not come in this algorithm...
    if ( nSize == 0) return;

    vector<int> cornerPos( nSize + 1);

    int index = 0;
    for (const JNodePtr &vtx: cornerNodes)
        cornerPos[index++] = getPosOf(vtx);
    cornerPos[nSize] = boundNodes.size();

    boost::sort(cornerPos);

    segSize.resize( nSize );

    for (int i = 0; i < nSize; i++) {
        segSize[i] = cornerPos[(i + 1)] - cornerPos[i];
    }
}

////////////////////////////////////////////////////////////////////

int QDefectivePatch::getTopologicalAngle(const JNodePtr &vertex)
{
    JFaceSequence vfaces;
    JNode::getRelations( vertex, vfaces );

    // How many faces are outside the regions.
    int ncount = 0;
    for( const JFacePtr &face : vfaces)
        if (faceSet.find(face) == faceSet.end()) ncount++;
    return ncount - 2;
}
////////////////////////////////////////////////////////////////////////////////

int QDefectivePatch::expandBlob(const JNodePtr &vertex)
{
    if( !vertex->isActive() )  return 1;
    if( vertex->isBoundary() ) return 1;

    if( faceSet.size() > maxSearchedFaces) return 1;

    JFaceSequence vfaces;
    JNode::getRelations(vertex, vfaces );

    assert( !vfaces.empty() ) ;

    for( const JFacePtr &face : vfaces) {
        assert(face->isActive() );
        if( face->hasBoundaryNode() ) {
            JFaceSet cset;
            boost::set_difference(faceSet, newFaceSet, std::inserter(cset, cset.end()));
            faceSet = cset;
            return 1;
        }
        if( faceSet.find(face) == faceSet.end() ) newFaceSet.insert(face);
        faceSet.insert(face);
        for (int j = 0; j < 4; j++) {
            const JNodePtr &vtx = face->getNodeAt(j);
            relations02[vtx].insert(face);
            const JEdgePtr &edge = face->getEdgeAt(j);
            edgeSet.insert(edge);
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int QDefectivePatch::expandBlob()
{
    cornerNodes.clear();
    boundEdges.clear();

    size_t nSize = boundNodes.size();
    assert( nSize > 0);

    for (size_t i = 0; i < nSize; i++) {
        int err = expandBlob(boundNodes[i]);
        if( err ) return 1;
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////

int
QDefectivePatch::buildTrianglesPatch()
{
    /*
    if( apexFace == nullptr) return 1;

    if( apexFace->isActive() == 0 || apexFace->getSize(0) != 3) return 1;

    faceSet.clear();

    JMeshDualGraph   mgraph;
    mgraph.setMesh(mesh);
    dualGraph = mgraph.getGraph();
    dualGraph->buildRelations(0,0);

    JNodePtr dualnode;
    apexFace->getAttribute("DualNode", dualnode);

    if( nonQuadFilter  == nullptr) nonQuadFilter.reset( new NonQuadFilter);
    djkPath.setMesh(dualGraph);
    djkPath.setFilter(nonQuadFilter);

       JNodeSequence dualnodes = djkPath.getPath(dualnode);

       JFacePtr pface;
       for( const JNodePtr &v : dualnodes) {
            v->getAttribute("PrimalFace", pface);
            cout << "Element Type " << pface->getSize(0) << endl;
            faceSet.insert(pface);
       }

       deque<JFacePtr>  faceQ;
       faceQ.push_back(apexFace);
       faceSet.insert( apexFace);

       JFaceSequence faceneighs;
       int numTriangles = 0;
       bool   patchfound = 0;
       while( !faceQ.empty() ) {
             const JFacePtr &currface = faceQ.front();
             faceQ.pop_front();
             if( currface->getSize(0) == 3) numTriangles++;
             if( numTriangles%2 == 0) {
                 patchfound = 1;
                 break;
             }
             JFace::getRelations12(currface, faceneighs);
             for( const JFacePtr &nextface : faceneighs) {
                  if( faceSet.find(nextface) == faceSet.end() ) {
                      faceSet.insert(nextface);
                      faceQ.push_back(nextface);
                  }
             }
        }

        if( patchfound == 0) return 1;

        for( const JFacePtr &face : faceSet) {
             edgeSet.insert(face->getEdgeAt(0) );
             edgeSet.insert(face->getEdgeAt(1) );
             edgeSet.insert(face->getEdgeAt(2) );
             edgeSet.insert(face->getEdgeAt(3) );
        }

    updateBoundary();
    finalizeBoundary();
    */

    return 0;

}

