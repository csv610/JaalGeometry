#include "DelaunayMesh.hpp"
#include <boost/lexical_cast.hpp>
#include "PointLocation.hpp"

using namespace Jaal;
using namespace std;

//////////////////////////////////////////////////////////////////////////////
void JDelaunayMesh2D :: build( vector<double> &uvCoords, vector<int> &segments,
                               vector<double> &holeCoords)
{
    struct triangulateio in, out; // triangle sw data structure

    // The following stuff is for "Triangle" software ...
    in.numberofpoints = uvCoords.size() / 2;
    in.pointlist = &uvCoords[0];

    in.numberofsegments = segments.size() / 2;
    if( segments.empty() )
        in.segmentlist = nullptr;
    else
        in.segmentlist = &segments[0];

    in.numberofholes = holeCoords.size()/2;
    if( holeCoords.empty() )
        in.holelist = nullptr;
    else
        in.holelist = &holeCoords[0];

    in.numberofpointattributes = 0;
    in.pointattributelist = nullptr;
    in.pointmarkerlist = nullptr;
    in.segmentmarkerlist = nullptr;
    in.triangleattributelist = nullptr;
    in.numberofregions = 0;
    in.regionlist = nullptr;

    // These needs to declared, I wasted full one day to find out this
    // restriction

    out.pointlist = nullptr;
    out.pointmarkerlist = nullptr;
    out.trianglelist = nullptr;
    out.neighborlist = nullptr;
    out.segmentlist = nullptr;
    out.segmentmarkerlist = nullptr;
    out.edgelist = nullptr;
    out.edgemarkerlist = nullptr;

    char *opt = const_cast<char*>( options.c_str() );
    triangulate(opt, &in, &out, (struct triangulateio *) nullptr);

    double uv[2];
    int nstart = in.numberofpoints;
    Point3D xyz;
    for (int i = nstart; i < out.numberofpoints; i++) {
        uv[0] = out.pointlist[2 * i + 0];
        uv[1] = out.pointlist[2 * i + 1];
        JNodePtr vtx = JNode::newObject();
        xyz[0] = uv[0];
        xyz[1] = uv[1];
        xyz[2] = 0.0;
        vtx->setXYZCoords(xyz);
        vtx->setID(i);
        newmesh->addObject(vtx);
    }

    JNodeSequence connect(3);
    JFacePtr face;
    for (int i = 0; i < out.numberoftriangles; i++) {
        int n1 = out.trianglelist[3 * i + 0];
        connect[0] = newmesh->getNodeAt(n1);

        int n2 = out.trianglelist[3 * i + 1];
        connect[1] = newmesh->getNodeAt(n2);

        int n3 = out.trianglelist[3 * i + 2];
        connect[2] = newmesh->getNodeAt(n3);

        face = JTriangle::newObject( connect );
        newmesh->addObject(face);
    }

    if (out.pointlist) free(out.pointlist);
    if (out.trianglelist) free(out.trianglelist);
}

////////////////////////////////////////////////////////////////////////

void JDelaunayMesh2D:: addSegments( JEdgeSequence  &boundedges)
{
    newmesh = JMesh::newObject();

    vector<int> segments;
    vector<double> uvCoords, holeCoords;

    JNodeSequence boundnodes;
    JMeshTopology::getEntitySet( boundedges, boundnodes);
    size_t numNodes = boundnodes.size();

    uvCoords.resize(2 * numNodes);
    size_t index = 0;

    map<JNodePtr,int> vmap;
    for( const JNodePtr &vtx: boundnodes) {
        const Point3D &xyz = vtx->getXYZCoords();
        uvCoords[2*index+0] = xyz[0];
        uvCoords[2*index+1] = xyz[1];
        vmap[vtx] = index++;
        newmesh->addObject(vtx);
    }

    size_t numEdges = boundedges.size();
    segments.resize(2 * numEdges);

    index = 0;
    for( const JEdgePtr &edge: boundedges) {
        newmesh->addObject(edge);
        segments[2*index+0] = vmap[edge->getNodeAt(0)];
        segments[2*index+1] = vmap[edge->getNodeAt(1)];
        index++;
    }
    build(uvCoords, segments, holeCoords);
}

//##############################################################################

void JDelaunayMesh2D :: addSegments( vector<JEdgeSequence>  &boundedges)
{
    if( boundedges.empty() ) return;

    if( boundedges.size() == 1 ) {
        addSegments(boundedges[0] );
        return;
    }

    int numLoops = boundedges.size();

    // Make sure that the loop is chained ...
    for( int i = 0; i < numLoops; i++) {
        int err = JEdgeTopology::getChain( boundedges[i] );
        assert( !err);
    }

    // Find out which is the outermost loop...
    vector<double> area(numLoops);
    JNodeSequence boundnodes;
    vector<Point2D> poly;
    for( int i = 0; i < numLoops; i++) {
        size_t numsegments = boundedges[i].size();
        poly.resize( 2*numsegments);
        for( size_t j = 0; j < numsegments; j++) {
            const Point3D &xyz = boundedges[i][j]->getNodeAt(0)->getXYZCoords();
            poly[j][0] = xyz[0];
            poly[j][1] = xyz[1];
        }
        area[i] = JGeometry::getSignedArea(poly);
    }

//    double maxarea   = 0.0;
    int    outerLoop = -1;
    for( int i = 0; i < numLoops; i++) {
        if( fabs(area[i]) > maxArea ) {
            maxArea = fabs(area[i]);
            outerLoop = i;
        }
    }

//  options = "BCPpzYQ";
    vector<double> holeCoords;
    Point3D xyz;
    for(int i = 0; i < numLoops; i++) {
        if( i != outerLoop) {
            addSegments( boundedges[i] );
            size_t numfaces = mesh->getSize(2);
            for( size_t j = 0; j < numfaces; j++) {
                const JFacePtr &face = mesh->getFaceAt(j);
                if( face->isActive() ) {
                    double facearea = JFaceGeometry::getArea(face);
                    if( fabs(facearea) > 0.0) {
                        face->getAvgXYZ(xyz);
                        holeCoords.push_back(xyz[0]);
                        holeCoords.push_back(xyz[1]);
                        break;
                    }
                }
            }
            mesh->deleteSurfaceMesh();
        }
    }

    size_t numNodes = 0;
    size_t numEdges = 0;
    for( size_t i = 0; i < boundedges.size(); i++) {
        numNodes += boundedges[i].size();
        numEdges += boundedges[i].size();
    }

    vector<double> uvCoords;
    uvCoords.resize(2 * numNodes);

    size_t index = 0;

    newmesh = JMesh::newObject();

    map<JNodePtr,int> vmap;
    for( size_t i = 0; i < boundedges.size(); i++) {
        JMeshTopology::getEntitySet( boundedges[i], boundnodes);
        for(size_t j = 0; j < boundnodes.size(); j++) {
            const JNodePtr &vtx = boundnodes[j];
            const Point3D &xyz = vtx->getXYZCoords();
            uvCoords[2*index+0] = xyz[0];
            uvCoords[2*index+1] = xyz[1];
            vmap[vtx] = index++;
            newmesh->addObject(vtx);
        }
    }

    vector<int> segments;
    segments.resize(2 * numEdges);
    index = 0;
    for( size_t i = 0; i < boundedges.size(); i++) {
        for( const JEdgePtr &edge: boundedges[i]) {
            newmesh->addObject(edge);
            segments[2*index+0] = vmap[edge->getNodeAt(0)];
            segments[2*index+1] = vmap[edge->getNodeAt(1)];
            index++;
        }
    }
//   holeCoords.clear();

    /*
        cout << uvCoords.size()/2 << "  2 0 0 " << endl;
        for( size_t i = 0; i < uvCoords.size()/2; i++)
            cout << i << "  " << uvCoords[2*i] << "  " << uvCoords[2*i+1] << endl;
        cout << segments.size()/2 <<  "   0 " << endl;
        for( size_t i = 0; i < segments.size()/2; i++)
            cout << i << "  " << segments[2*i] << "  " << segments[2*i+1] << endl;
        cout << holeCoords.size()/2 << endl;
        for( size_t i = 0; i < holeCoords.size()/2; i++)
            cout << i << "  " << holeCoords[2*i] << " " << holeCoords[2*i+1] << endl;
    */

    build(uvCoords, segments, holeCoords);

    newmesh->getTopology()->searchBoundary();

}

////////////////////////////////////////////////////////////////////////

void JDelaunayMesh2D :: addPoints( const vector<Point2D> &pointCloud)
{
    newmesh = JMesh::newObject();
    vector<int> segments;
    vector<double> uvCoords, holeCoords;

    int numNodes = pointCloud.size();

    uvCoords.resize(2 * numNodes);
    size_t index = 0;

    Point3D xyz;
    for( Point2D pd: pointCloud) {
        uvCoords[2*index+0] = pd[0];
        uvCoords[2*index+1] = pd[1];
        xyz[0] = pd[0];
        xyz[1] = pd[1];
        xyz[2] = 0.0;
        JNodePtr vtx = JNode::newObject();
        vtx->setID(index);
        vtx->setXYZCoords(xyz);
        mesh->addObject(vtx);
        index++;
    }

    build(uvCoords, segments, holeCoords);
}

//##############################################################################
bool JDelaunayMesh2D :: isDelaunay(const JFacePtr &face)
{
    int nnodes = face->getSize(0);
    assert(nnodes == 3);

    static JNodeSet  vset;
    static JNodeSequence  vneighs;

    // Collect all the neigboring nodes except those defining the face..
    for( int i = 0; i < nnodes; i++) {
        const JNodePtr &vertex = face->getNodeAt(i);
        JNode::getRelations( vertex, vneighs );
        for( size_t j = 0; j < vneighs.size(); j++)
            vset.insert( vneighs[j] );
    }

    double side;
    const Point3D &pa = face->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = face->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = face->getNodeAt(2)->getXYZCoords();

    JNodeSet::const_iterator it;
    for( it = vset.begin(); it != vset.end(); ++it) {
        const Point3D &pd = (*it)->getXYZCoords();
        double dir = orient2d( const_cast<double*>(&pa[0]),
                               const_cast<double*>(&pb[0]),
                               const_cast<double*>(&pc[0]) );
        if( dir != 0.0) {
            if( dir >  0.0) {
                side = incircle( const_cast<double*>(&pa[0]),
                                 const_cast<double*>(&pb[0]),
                                 const_cast<double*>(&pc[0]),
                                 const_cast<double*>(&pd[0]) );
                if( side > 0.0) return 0;
            } else {
                side = incircle( const_cast<double*>(&pa[0]),
                                 const_cast<double*>(&pc[0]),
                                 const_cast<double*>(&pb[0]),
                                 const_cast<double*>(&pd[0]) );
                if( side > 0.0) return 0;
            }
        }
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

bool JDelaunayMesh2D :: isDelaunay(const JEdgePtr &edge)
{
    if( edge == nullptr ) return  0;

    JFaceSequence faceneighs;
    JEdge::getRelations(edge, faceneighs);

    if( faceneighs.size()  < 1 ) {
        cout << "Warning: Edge-faces relations empty" <<endl;
        return 0;
    }

    const JNodePtr &v0 = edge->getNodeAt(0);
    const JNodePtr &v1 = edge->getNodeAt(1);

    if( faceneighs.size() == 2 ) {
        const JNodePtr &v2 = JTriangle::getOppositeNode(faceneighs[0], v0, v1);
        const JNodePtr &v3 = JTriangle::getOppositeNode(faceneighs[1], v0, v1);
        double l0 = JMath::length2(v0->getXYZCoords(),v1->getXYZCoords());
        double l1 = JMath::length2(v2->getXYZCoords(),v3->getXYZCoords());
        if( l0 < l1 ) return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
size_t JDelaunayMesh2D :: countNonDelaunayEdges()
{
    if( mesh == nullptr) return 0;

    size_t ncount = 0;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge  = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            if( isDelaunay( edge ) == 0 ) ncount++;
        }
    }
    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

size_t JDelaunayMesh2D :: countNonDelaunayFaces()
{
    if( mesh == nullptr) return 0;

    mesh->buildRelations(0,2);

    size_t ncount = 0;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face  = mesh->getFaceAt(i);
        if( face->isActive() ) {
            if( isDelaunay( face ) == 0 ) ncount++;
        }
    }
    mesh->clearRelations(0,2);
    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int JDelaunayMesh2D :: isDelaunay()
{
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge  = mesh->getEdgeAt(i);
        if( !isDelaunay( edge ) ) return 0;
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr JDelaunayMesh2D :: getSimpleMesh()
{
    if( mesh == nullptr ) return nullptr;

    options = "BCPpzYQ";
    vector<JEdgeSequence> edges;

    int dim = mesh->getTopology()->getDimension();

    if( dim == 2 ) {
        mesh->getTopology()->getBoundary(edges);
        mesh->deleteSurfaceMesh();
    }

    if( dim == 1 ) {
        JEdgeSequence alledges = mesh->getEdges();
        JEdgeTopology::getLoops(alledges, edges);
    }

    if( edges.empty() ) {
        cout << "Warning: Boundary edges not found: Delaunay meshing not performed " << endl;
        return nullptr;
    }

    addSegments(edges);
//  JMeshIO::saveAs(newmesh, "tmp.off");

    return newmesh;
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr JDelaunayMesh2D :: getQualityMesh()
{
    if( mesh == nullptr ) return nullptr;

    vector<JEdgeSequence> edges;
    mesh->getTopology()->getBoundary( edges);
    mesh->deleteSurfaceMesh();
    ostringstream oss;
    oss << "BCPpzQq" << minAngle;
    if( maxArea > 0.0) oss << "a" << maxArea;
    if( !boundarySplit)  oss << "Y";
    options = oss.str();
    addSegments(edges);
    return newmesh;
}
///////////////////////////////////////////////////////////////////////////////
JMeshPtr JDelaunayMesh2D :: getConvexHull()
{
    cout << "Not yet implemented " << endl;
    return nullptr;
}
///////////////////////////////////////////////////////////////////////////////
JMeshPtr JDelaunayMesh2D :: getMedialAxis()
{
    if( mesh == nullptr ) return nullptr;

    JMeshPtr tmpmesh = mesh->deepCopy();
    mesh =  tmpmesh;
    getSimpleMesh();

    JPointLocation pointlocater;
    pointlocater.setMesh(newmesh);

    size_t numfaces = newmesh->getSize(2);
    if( numfaces < 1) return nullptr;

    Point3D center;
    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;
    JNodeSequence nodes(numfaces);
    JFaceSequence inface(numfaces);

    newmesh->pruneAll();
    newmesh->enumerate(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = newmesh->getFaceAt(i);
        const Point3D &a = face->getNodeAt(0)->getXYZCoords();
        const Point3D &b = face->getNodeAt(1)->getXYZCoords();
        const Point3D &c = face->getNodeAt(2)->getXYZCoords();
        TriCircumCenter2D( &a[0], &b[0], &c[0], &center[0] );
        JNodePtr vtx = JNode::newObject();
        vtx->setID(i);
        vtx->setXYZCoords(center);
        nodes[i] = vtx;
        inface[i]   = pointlocater.searchFace(center,1);
    }

    JMeshPtr  medialmesh = JMesh::newObject();
    medialmesh->addObjects(nodes);

    size_t numedges = newmesh->getSize(1);
    JFaceSequence neighs;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = newmesh->getEdgeAt(i);
        JEdge::getRelations(edge, neighs);
        if( neighs.size() == 2 ) {
            int fid1 = neighs[0]->getID();
            int fid2 = neighs[1]->getID();
            if( inface[fid1] && inface[fid2] ) {
                JEdgePtr dualedge = JEdge::newObject( nodes[fid1], nodes[fid2] );
                medialmesh->addObject(dualedge);
            }
        }
    }
    return medialmesh;
}

///////////////////////////////////////////////////////////////////////////////
int JDelaunayMesh2D :: getRemeshed()
{
    if( mesh == nullptr) {
        cout << "Warning: A null mesh pointer was passed to Delaunay Mesher" << endl;
        return 1;
    }

    JSwapTriEdge swapper;
    swapper.setMesh(mesh);
    swapper.setCreaseAngle(creaseAngle);

    assert( mesh->getTopology()->isManifold() );

    int totalflipped = 0;
    while(1) {
        size_t numflipped = swapper.execute();
        if( numflipped == 0) break;
        totalflipped += numflipped;
        mesh->pruneAll();
    }

    assert( mesh->getTopology()->isManifold() );

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            if( JFaceGeometry::getSignedArea(face) < 0.0)
                face->reverse();
        }
    }
    return totalflipped;
}
///////////////////////////////////////////////////////////////////////////////

int JDelaunayMesh2D :: getIntrinsicMesh()
{
    if( mesh == nullptr) {
        cout << "Warning: A null mesh pointer was passed to Delaunay Mesher" << endl;
        return 1;
    }

    JSwapTriEdge swapper;
    swapper.setMesh(mesh);

    int totalflipped = 0;
    while(1) {
        size_t numflipped = swapper.execute();
        if( numflipped == 0) break;
        totalflipped += numflipped;
        mesh->pruneAll();
    }
    return totalflipped;
}


#ifdef CSV
///////////////////////////////////////////////////////////////////////////////

JMeshPtr DelaunayMesh :: toTriVoronoi(const JMeshPtr &mesh)
{
    JMeshPtr voromesh = JMesh::newObject();

    double r[3], p[3];
    size_t numfaces = mesh->getSize(2);
    Point3D p3d;

    JNodeSequence ( numfaces );

    /*
        for( size_t i = 0; i < numfaces; i++)
        {
            Face *face  = mesh->getFaceAt(i);
            assert( face->getSize(0) == 3 );
            Point3D pa = face->getNodeAt(0)->getXYZCoords();
            Point3D pb = face->getNodeAt(1)->getXYZCoords();
            Point3D pc = face->getNodeAt(2)->getXYZCoords();
            TriCircumCenter2D(&pa[0], &pb[0], &pc[0], r, p);
            Vertex *voronode = Vertex::newObject();
            p3d[0] = r[0];
            p3d[1] = r[1];
            p3d[2] = 0.0;
            voronode->setXYZCoords( p3d );
            voromesh->addNode(voronode);
            face->setID(i);       // Necessary here...
            voronode->setID(i);
            facenode[i] = voronode;
        }

        int relexist0 =  mesh->buildRelations(0,0);
        int relexist2 =  mesh->buildRelations(0,2);

        size_t numnodes = mesh->getSize(0);
        vector<Edge>  edges;
        for( size_t i = 0; i <  numnodes; i++){
             Vertex *apex = mesh->getNodeAt(i);
             vector<Vertex*> vneighs = apex->getRelations0();
             edges.clear();
             int complete_loop = 1;
             for( int j = 0; j < vneighs.size(); j++) {
                  vector<Face*> efaces = Mesh::getRelations112( apex, vneighs[j]);
                  if( efaces.size() != 2 ) {
                      complete_loop = 0;
                      break;
                  }
                  Vertex *v0 = facenode[ efaces[0]->getID() ];
                  Vertex *v1 = facenode[ efaces[1]->getID() ];
                  edges.push_back( Edge(v0,v1) );
             }

             if( complete_loop ) {
                 Mesh::make_chain( edges );
                 Face *voroface = Face::newObject();
                 voroface->setNodes( Mesh::chain_nodes( edges ));
                 voromesh->addFace(voroface);
             }
         }

        if( !relexist0 ) mesh->clearRelations(0,0);
        if( !relexist2 ) mesh->buildRelations(0,2);

         return voromesh;
    */
}
///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JDelaunayMesh2D :: getConvexHull( const vector<Point2D> &pointCloud)
{
    /*
        mesh = JMesh::newObject();

        vector<int> segments;
        vector<double> uvCoords, holeCoords;

        int numNodes = pointCloud.size();

        uvCoords.resize(2 * numNodes);
        size_t index = 0;

        Point3D xyz;
        for( Point2D pd: pointCloud) {
            uvCoords[2*index+0] = pd[0];
            uvCoords[2*index+1] = pd[1];
            xyz[0] = pd[0];
            xyz[1] = pd[1];
            xyz[2] = 0.0;
            JNodePtr vtx = JNode::newObject();
            vtx->setID(index);
            vtx->setXYZCoords(xyz);
            mesh->addObject(vtx);
            index++;
        }

        build(uvCoords, segments, holeCoords);

        JEdgeSequence boundedges;
        mesh->getTopology()->getBoundary(boundedges);

        JNodeSequence boundnodes;
        JMeshTopology::getEntitySet( boundedges, boundnodes);

        map<JNodePtr,JNodePtr> vmap;
        for( JNodePtr oldvtx: boundnodes)
            vmap[oldvtx] = oldvtx->getClone();

        JEdgeSequence hulledges;
        for( JEdgePtr oldedge: boundedges) {
            JNodePtr v0 = oldedge->getNodeAt(0);
            JNodePtr v1 = oldedge->getNodeAt(1);
            JEdgePtr newedge = JEdge::newObject( vmap[v0], vmap[v1] );
            hulledges.push_back( newedge);
        }
        mesh->deleteAll();
    //  delete mesh;
        mesh = nullptr;
        return hulledges;
    */
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////

