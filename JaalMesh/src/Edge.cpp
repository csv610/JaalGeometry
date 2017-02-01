#include "Mesh.hpp"
#include "GeomPredicates.hpp"

using namespace Jaal;

size_t JEdge::NumObjectsCreated = 0;
std::map<string,string> JEdge::attribInfo;

///////////////////////////////////////////////////////////////////////////////
int JEdge :: registerAttribute( const string &name, const string &type)
{
    int  found = 0;

    if( type =="int"    ) found = 1;
    if( type =="char"   ) found = 1;
    if( type =="float"  ) found = 1;
    if( type =="double" ) found = 1;
    if( type =="uchar"  ) found = 1;

    if( !found) {
        cout << "Warning: invalid attribute type " << type << endl;
        return 2;
    }

    attribInfo[name] = type ;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

string JEdge :: getAttributeTypeName( const string &name)
{
    string str;
    if( attribInfo.find( name ) == attribInfo.end() ) return str;
    return  attribInfo[name];
}
///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JEdge::newObjects(size_t n)
{
    JEdgeSequence edges;
    if( n < 1) return edges;
    edges.resize(n);
    for( size_t i = 0; i < n; i++)
        edges[i] = JEdge::newObject();
    return edges;
}
////////////////////////////////////////////////////////////////////////////////
void JEdge :: getRelations( const JEdgePtr &edge, JFaceSequence &seq)
{
    seq.clear();
    if( edge == nullptr) return;
    edge->getRelations_(seq);
}
////////////////////////////////////////////////////////////////////////////////

void JEdge :: getRelations( const JEdgePtr &edge, JCellSequence &seq)
{
    seq.clear();
    if( edge == nullptr) return;
    edge->getRelations_(seq);
}
///////////////////////////////////////////////////////////////////////////////

bool JEdgeTopology :: lexiCompare(const JEdgePtr &edge1, const JEdgePtr &edge2)
{
    int v0 = edge1->getNodeAt(0)->getID();
    int v1 = edge2->getNodeAt(0)->getID();

    if( v0 == v1) {
        int v2 = edge1->getNodeAt(1)->getID();
        int v3 = edge2->getNodeAt(1)->getID();
        return v2 < v3 ? 1 : 0;
    }

    return v0 < v1 ? 1:0;
}


void
JEdgeTopology::reverse( JEdgeSequence &boundedges)
{
    std::reverse( boundedges.begin(), boundedges.end() );
    // Reverse all the connection
    size_t nSize = boundedges.size();
    for( size_t i = 0; i < nSize; i++)
        boundedges[i]->reverse();
}

///////////////////////////////////////////////////////////////////////////////

void
JEdgeGeometry::generateLinearNodes( const JNodePtr &v0, const JNodePtr &v1, int n, JNodeSequence &nodes)
{
    if( !v0->isActive() ) return;
    if( !v1->isActive() ) return;

    assert(n >= 2);

    const Point3D &xyz0 = v0->getXYZCoords();
    const Point3D &xyz1 = v1->getXYZCoords();

    JNodeSequence newnodes;
    JMesh::generate_objects( n-2, newnodes);
    assert( newnodes.size() == size_t(n-2) );

    nodes.resize(n);

    nodes[0] = v0;
    nodes[n - 1] = v1;

    if( n == 2) return;

    double dt = 2.0 / (double) (n - 1);

    Point3D xyzt;
    for (int i = 1; i < n - 1; i++) {
        double t = -1.0 + i*dt;
        xyzt[0] = TFI::linear_interpolation(t, xyz0[0], xyz1[0]);
        xyzt[1] = TFI::linear_interpolation(t, xyz0[1], xyz1[1]);
        xyzt[2] = TFI::linear_interpolation(t, xyz0[2], xyz1[2]);
        nodes[i] = newnodes[i-1];
        nodes[i]->setXYZCoords(xyzt);
    }
}

/////////////////////////////////////////////////////////////////////////////////
void
JEdgeGeometry::generateLinearNodes( const JEdgePtr &edge, int n, JNodeSequence &newnodes)
{
    newnodes.clear();

    if( !edge->isActive() ) return;

    const JNodePtr &v0 = edge->getNodeAt(0);
    const JNodePtr &v1 = edge->getNodeAt(1);

    JNodeSequence nodes;
    generateLinearNodes(v0,v1,n,nodes);

    if( n > 2) {
        newnodes.resize(n-2);
        for( int i = 0; i < n-2; i++)
            newnodes[i] = nodes[i+1];
        edge->setAttribute("Steiner", newnodes);
    }
}



///////////////////////////////////////////////////////////////////////////////
double JEdgeGeometry :: getCreaseAngle( const JEdgePtr &edge)
{
    double angle = 0.0;
    /*
       JFaceSequence faces;
       JEdge::getRelations( edge, faces);
       if( faces.size() != 2 ) {
           cout << "Warning: An edge is not shared nu two faces " << endl;
           return 0;
       }
       Vec3F avec, bvec;
       int err;
       err = faces[0]->getAttribute("Normal", avec);
       if( err ) JFaceGeometry::getNormal(faces[0], avec);

       err = faces[1]->getAttribute("Normal", bvec);
       if( err ) JFaceGeometry::getNormal(faces[1], bvec);

       double angle = math::getVecAngle( avec, bvec, ANGLE_IN_DEGREES);
    */

    return angle;
}

///////////////////////////////////////////////////////////////////////////////
int JEdgeGeometry :: intersectPos2d( const double *p1, const double *p2,
                                     const double *p3, const double *p4,
                                     double *xy, double *uv)
{
    // Before calling this functiom, call robust function to determine
    // the intersection...
    xy[0] = 0.0;
    xy[1] = 0.0;

    double x1 = p1[0];
    double y1 = p1[1];
    double x2 = p2[0];
    double y2 = p2[1];
    double x3 = p3[0];
    double y3 = p3[1];
    double x4 = p4[0];
    double y4 = p4[1];

    long double denom = (x2-x1)*(y4-y3) - (y2-y1)*(x4-x3);
    if (denom == 0.0) return 1;

    long double ua =  ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3))/denom;
    long double ub =  ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3))/denom;

    double eps = 1.0E-06;
    if( ua < -eps)     return -1;
    if( ua > 1.0+eps)  return -1;
    if( ua < 0.0) ua = 0.0;
    if( ua > 1.0) ua = 1.0;

    if( ub < -eps)     return -1;
    if( ub > 1.0+eps)  return -1;
    if( ub < 0.0) ub = 0.0;
    if( ub > 1.0) ub = 1.0;

    if( xy != nullptr ) {
        xy[0] = x1 + ua*(x2-x1);
        xy[1] = y1 + ua*(y2-y1);
    }

    if( uv != nullptr) {
        uv[0] = ua;
        uv[1] = ub;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool JEdgeGeometry :: isWithin( const double *p0, const double *p1, const double *pq)
{
    double xmin1 = std::min( p0[0], p1[0]);
    double xmax1 = std::max( p0[0], p1[0]);

    double ymin1 = std::min( p0[1], p1[1]);
    double ymax1 = std::max( p0[1], p1[1]);

    if( pq[0] < xmin1 || pq[0] > xmax1) return 0;
    if( pq[1] < ymin1 || pq[1] > ymax1) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////
double JEdgeGeometry :: getU( const double *p0, const double *p1, const double *xy)
{
    //
    // If the point lie on the straight line, return the parametric value with
    // respect to the first point.
    //
    double dx, dy;
    dx = xy[0] - p0[0];
    dy = xy[1] - p0[1];
    double num = sqrt(dx*dx + dy*dy);

    dx = p1[0] - p0[0];
    dy = p1[1] - p0[1];
    double dem = sqrt(dx*dx + dy*dy);

    return num/dem;
}
///////////////////////////////////////////////////////////////////////////////
Point3D JEdgeGeometry :: getMidPoint( const JEdgePtr &edge, double alpha)
{
    return JNodeGeometry::getMidPoint( edge->getNodeAt(0), edge->getNodeAt(1), alpha);
}

JNodePtr JEdgeGeometry :: getMidNode( const JEdgePtr &edge, double alpha)
{
    return JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1), alpha);
}

int JEdgeGeometry :: intersectPredicate2d( const double *p0, const double *p1,
        const double *p2, const double *p3)
{
    ///////////////////////////////////////////////////////////////////////////
    // Checks if the 2D line segments cross each other. This is a robust
    // function based on Jonathan Shewchuls geometrric predicates..
    //
    // Return value:
    // -1  :  Segement do not intersect.
    //  0  :  one of the point is on the line segment, it is upto the
    //        application to decide, if it is inside or not.
    //  1  :  It crosses in the conventional sense...
    ///////////////////////////////////////////////////////////////////////////

    double xmin1 = std::min( p0[0], p1[0]);
    double xmax1 = std::max( p0[0], p1[0]);

    double ymin1 = std::min( p0[1], p1[1]);
    double ymax1 = std::max( p0[1], p1[1]);

    double xmin2 = std::min( p2[0], p3[0]);
    double xmax2 = std::max( p2[0], p3[0]);

    double ymin2 = std::min( p2[1], p3[1]);
    double ymax2 = std::max( p2[1], p3[1]);

    // First simplest check.
    if( xmin2 > xmax1 || xmin1 > xmax2 ) return -1;
    if( ymin2 > ymax1 || ymin1 > ymax2 ) return -1;

    // End points are same ...
    if( p0[0] == p2[0] && p0[1] == p2[1] ) return 0;
    if( p0[0] == p3[0] && p0[1] == p3[1] ) return 0;
    if( p1[0] == p2[0] && p1[1] == p2[1] ) return 0;
    if( p1[0] == p3[0] && p1[1] == p3[1] ) return 0;

    int ori1, ori2;
    // If one segment is on one side then there is no intersection.
    ori1 = JGeomPredicates::getPointOrientation( p0, p1, p2);
    if( ori1 == 0) {
        if( p2[0] < xmin1 || p2[0] > xmax1) return -1.0;
        if( p2[1] < ymin1 || p2[1] > ymax1) return -1.0;
        return 0;
    }
    ori2 = JGeomPredicates::getPointOrientation( p0, p1, p3);
    if( ori2 == 0) {
        if( p3[0] < xmin1 || p3[0] > xmax1) return -1.0;
        if( p3[1] < ymin1 || p3[1] > ymax1) return -1.0;
        return 0;
    }
    if( ori1*ori2 > 0) return -1;

    ori1 = JGeomPredicates::getPointOrientation( p2, p3, p0);
    if( ori1 == 0) {
        if( p0[0] < xmin2 || p0[0] > xmax2) return -1.0;
        if( p0[1] < ymin2 || p0[1] > ymax2) return -1.0;
        return 0;
    }
    ori2 = JGeomPredicates::getPointOrientation( p2, p3, p1);
    if( ori2 == 0) {
        if( p1[0] < xmin2 || p1[0] > xmax2) return -1.0;
        if( p1[1] < ymin2 || p1[1] > ymax2) return -1.0;
        return 0;
    }
    if( ori1*ori2 > 0) return -1;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int JEdgeGeometry :: intersectPredicate2d( const JEdgePtr &edge1, const JEdgePtr &edge2)
{
    int err = -1;

    if( edge1 != edge2) {
        const Point3D &p0 = edge1->getNodeAt(0)->getXYZCoords();
        const Point3D &p1 = edge1->getNodeAt(1)->getXYZCoords();

        const Point3D &p2 = edge2->getNodeAt(0)->getXYZCoords();
        const Point3D &p3 = edge2->getNodeAt(1)->getXYZCoords();

        err = intersectPredicate2d( &p0[0], &p1[0], &p2[0], &p3[0] );
    }
    return err;
}
///////////////////////////////////////////////////////////////////////////////

int JEdgeGeometry :: intersectPos2d( const JEdgePtr &edge1, const JEdgePtr &edge2,
                                     double *xy, double *uv)
{
    int  stat = -1;

    if( edge1 != edge2) {
        const Point3D &p0 = edge1->getNodeAt(0)->getXYZCoords();
        const Point3D &p1 = edge1->getNodeAt(1)->getXYZCoords();

        const Point3D &p2 = edge2->getNodeAt(0)->getXYZCoords();
        const Point3D &p3 = edge2->getNodeAt(1)->getXYZCoords();

        stat = intersectPredicate2d( &p0[0], &p1[0], &p2[0], &p3[0] );
        if( stat < 0) return -1;

        stat = intersectPos2d( &p0[0], &p1[0], &p2[0], &p3[0], xy, uv );
    }

    return stat;
}
///////////////////////////////////////////////////////////////////////////////

void JEdgeGeometry :: smooth(JEdgeSequence &edges, int numIterations)
{
    map<JNodePtr, Point3D> newPos;
    map<JNodePtr, int>     vcount;

    Point3D xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    for( const JEdgePtr &edge: edges) {
        const JNodePtr &v0 = edge->getNodeAt(0);
        if( vcount.find(v0) == vcount.end() )
            vcount[v0] = 1;
        else
            vcount[v0]++;

        const JNodePtr &v1 = edge->getNodeAt(1);
        if( vcount.find(v1) == vcount.end() )
            vcount[v1] = 1;
        else
            vcount[v1]++;
        newPos[v0] = xyz;
        newPos[v1] = xyz;
    }

    Point3D pmid;
    for( int i = 0; i < numIterations; i++) {
        for( auto &key : newPos) key.second = xyz;
        for( const JEdgePtr &edge: edges) {
            const JNodePtr v0 = edge->getNodeAt(0);
            const JNodePtr v1 = edge->getNodeAt(1);
            const Point3D &p0 = v0->getXYZCoords();
            const Point3D &p1 = v1->getXYZCoords();
            newPos[v0][0] += p1[0];
            newPos[v0][1] += p1[1];
            newPos[v0][2] += p1[2];

            newPos[v1][0] += p0[0];
            newPos[v1][1] += p0[1];
            newPos[v1][2] += p0[2];
        }
        for( auto &key: newPos) {
            JNodePtr vtx = key.first;
            int nsize = vcount[vtx];
            if( nsize > 1) {
                pmid[0]   = newPos[vtx][0]/( double)vcount[vtx];
                pmid[1]   = newPos[vtx][1]/( double)vcount[vtx];
                pmid[2]   = newPos[vtx][2]/( double)vcount[vtx];
                vtx->setXYZCoords(pmid);
            }
        }
    }

}
/////////////////////////////////////////////////////////////////////////////////

void
JEdgeGeometry::generateLinearNodes( const JEdgePtr &edge, int n, const JMeshPtr &mesh)
{
    JNodeSequence newnodes;
    generateLinearNodes( edge, n, newnodes);
    if( mesh ) mesh->addObjects(newnodes);
}

/////////////////////////////////////////////////////////////////////////////////

int JEdgeTopology :: getLoops( const JEdgeSequence &edges, vector<JEdgeSequence> &loops)
{
    loops.clear();
//    size_t numEdges = edges.size();
    for( const JEdgePtr &edge : edges)
        edge->setVisitBit(0);

    JEdgeSequence newloop;
    JEdgePtr  curredge, nextedge;

    while(1) {
        newloop.clear();
        curredge.reset();
        for( const JEdgePtr &edge : edges)  {
            if( edge->getVisitBit() == 0) {
                curredge = edge;
                break;
            }
        }
        if( curredge == nullptr) break;

        const JNodePtr &startnode = curredge->getNodeAt(0);
        JNodePtr currnode  = curredge->getNodeAt(1);

        while(1) {
            newloop.push_back(curredge);
            curredge->setVisitBit(1);
            if( currnode == startnode) break;
            nextedge.reset();
            for( const JEdgePtr &edge : edges)  {
                if( edge->getVisitBit() == 0) {
                    if( edge->getNodeAt(0) == currnode) {
                        nextedge = edge;
                        currnode = edge->getNodeAt(1);
                        break;
                    }
                    if( edge->getNodeAt(1) == currnode) {
                        edge->reverse();
                        nextedge = edge;
                        currnode = edge->getNodeAt(1);
                        break;
                    }
                }
            }
            if( nextedge == nullptr) break;
            curredge = nextedge;
        }
        if( !newloop.empty() ) {
            loops.push_back(newloop);
        }
    }

    for( const JEdgePtr &edge : edges) {
        if( edge->getVisitBit() == 0) {
            cout << "Error: Some edges not visited: Edge loops may be wrong " << endl;
            return 1;
        }
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////
bool
JEdgeTopology::isChain( const JEdgeSequence &edges)
{
    // Determine if the nodes in the consecutive edges are connected.
    // If the edges are closed, then the first node of the first segment and
    // second node of the last segment are the same. If the edges are open then
    // the first node of the first segment and the second node of the last segment
    // have only one neighbor. Intermediate nodes can have any number of the neighbors,
    // but, if the same is simple, then every internal node will have only two
    // neighbouring nodes.

    size_t numSegments = edges.size();

    map<JNodePtr, int> nodeCount;
    JNodePtr vertex;
    for( size_t i = 0; i < numSegments; i++) {
        vertex = edges[i]->getNodeAt(0);
        nodeCount[vertex] = 0;
        vertex = edges[i]->getNodeAt(1);
        nodeCount[vertex] = 0;
    }

    for( size_t i = 0; i < numSegments; i++) {
        vertex = edges[i]->getNodeAt(0);
        nodeCount[vertex]++;
        vertex = edges[i]->getNodeAt(1);
        nodeCount[vertex]++;
    }

    JNodePtr startVertex;
    map<JNodePtr,int>::const_iterator it;
    for( it = nodeCount.begin(); it != nodeCount.end(); ++it)
        if( it->second == 1 ) startVertex = it->first;

    if( startVertex == nullptr) {
        for( size_t i = 0; i < numSegments; i++) {
            const JNodePtr &v1 = edges[i]->getNodeAt(1);
            const JNodePtr &v0 = edges[(i+1)%numSegments]->getNodeAt(0);
            if( v0 != v1) return 0;
        }
        return 1;
    }

    for( size_t i = 0; i < numSegments; i++) {
        const JNodePtr &v1 = edges[i]->getNodeAt(1);
        const JNodePtr &v0 = edges[(i+1)%numSegments]->getNodeAt(0);
        if( v0 != v1) return 0;
    }

    return 1;

}
/////////////////////////////////////////////////////////////////////////////////

int
JEdgeGeometry::getOrientation( const JEdgeSequence &edges)
{
    if( !JEdgeTopology::isChain(edges) ) return 0;

    int n = edges.size();

    vector<double>  x(n), y(n);
    int index = 0;
    for( const JEdgePtr &e : edges) {
        const Point3D &p = e->getNodeAt(0)->getXYZCoords();
        x[index] = p[0];
        y[index] = p[1];
        index++;
    }
    double area = JGeometry::getSignedArea( &x[0], &y[0], n);
    if( area < 0.0 ) return -1;
    if( area > 0.0 ) return  1;

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////

int
JEdgeTopology::getChain( JEdgeSequence &edges)
{
    // Starting from the first node of the first edge,
    // link the edge segments to make a chain.

    size_t nSize = edges.size();
    if( nSize == 0) return 1;

    JNodePtr curr_vertex = edges[0]->getNodeAt(1);

    for( size_t i = 1; i < nSize; i++) {
        for( size_t  j = i; j < nSize; j++) {
            const JNodePtr &v0 = edges[j]->getNodeAt(0);
            const JNodePtr &v1 = edges[j]->getNodeAt(1);
            if (v0 == curr_vertex) {
                curr_vertex = v1;
                if( j > i) swap( edges[i], edges[j] );
                break;  // now go to the next segment.
            }

            if (v1 == curr_vertex) {
                curr_vertex = v0;
                edges[j]->reverse();
                if( j > i ) swap( edges[i], edges[j] );
                break; // now go to the next segment.
            }
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
JEdgeTopology::getChain(JEdgeSequence &edges, const JEdgePtr &start_edge)
{
    // starting from the "given" edge, make a chain.

    JEdgeSequence::iterator it;
    it = find(edges.begin(), edges.end(), start_edge);
    if( it == edges.end() ) return 1;

    std::rotate( edges.begin(), it, edges.end() );

    assert( edges[0] == start_edge );
    int stat = JEdgeTopology::getChain( edges);
    return stat;
}

///////////////////////////////////////////////////////////////////////////////
int
JEdgeTopology::getChain(JEdgeSequence &edges, const JNodePtr &start_node)
{
    // Starting from the given vertex, build the chain...
    size_t numedges = edges.size();

    JEdgePtr start_edge = nullptr;
    for( size_t i = 0; i < numedges; i++) {
        if( edges[i]->getNodeAt(0) == start_node) {
            start_edge = edges[i];
            break;
        }
    }

    if( start_edge == nullptr) {
        for( size_t i = 0; i < numedges; i++) {
            if( edges[i]->getNodeAt(1) == start_node) {
                start_edge = edges[i];
                start_edge->reverse();
                break;
            }
        }
    }

    if( start_edge == nullptr ) return 1;

    int stat = JEdgeTopology::getChain(edges, start_edge);
    return stat;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
JEdgeTopology::isCloseable(const JEdgeSequence &edges)
{
    // If all the nodes have atleast two edges, then it is probably closeable.
    // In general, this could be more complicated, but not for our cases, at
    // present...
    std::map<JNodePtr, JEdgeSequence> relations00;

    for (size_t i = 0; i < edges.size(); i++) {
        JNodePtr v0 = edges[i]->getNodeAt(0);
        JNodePtr v1 = edges[i]->getNodeAt(1);
        relations00[v0].push_back(edges[i] );
        relations00[v1].push_back(edges[i] );
    }

    std::map<JNodePtr, JEdgeSequence> ::const_iterator it;
    for (it = relations00.begin(); it != relations00.end(); ++it) {
        JNodePtr v = it->first;
        if (relations00[v].size() != 2) return 0;
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

bool
JEdgeTopology::isTopologicalSimple(const JEdgeSequence &edges)
{
    // In a topological simple, edge curve, every node is atmost has two node
    // neigbours...

    size_t numSegments = edges.size();

    map<JNodePtr, int> nodeCount;
    JNodePtr vertex;
    for( size_t i = 0; i < numSegments; i++) {
        vertex = edges[i]->getNodeAt(0);
        nodeCount[vertex] = 0;
        vertex = edges[i]->getNodeAt(1);
        nodeCount[vertex] = 0;
    }

    for( size_t i = 0; i < numSegments; i++) {
        vertex = edges[i]->getNodeAt(0);
        nodeCount[vertex]++;
        vertex = edges[i]->getNodeAt(1);
        nodeCount[vertex]++;
    }

    map<JNodePtr,int>::const_iterator it;
    for( it = nodeCount.begin(); it != nodeCount.end(); ++it)
        if( it->second > 2) return 0;
    return 1;
}

///////////////////////////////////////////////////////////////////////////////
int
JEdgeTopology::isClosed(const JEdgeSequence  &boundedges)
{
    // Check if the chain is topological circle ...
    JNodePtr first_vertex = boundedges[0]->getNodeAt(0);
    JNodePtr curr_vertex  = boundedges[0]->getNodeAt(1);

    for (size_t i = 1; i < boundedges.size(); i++) {
        if (boundedges[i]->getNodeAt(0) != curr_vertex) return 0;
        curr_vertex = boundedges[i]->getNodeAt(1);
    }
    if (curr_vertex != first_vertex) return 0;

    return 1;
}
///////////////////////////////////////////////////////////////////////////////

int
JEdgeTopology::getChainNodes(const JEdgeSequence &edges, JNodeSequence &bndnodes)
{
    bndnodes.clear();

    JNodePtr start_node = edges[0]->getNodeAt(0);
    bndnodes.push_back(start_node);

    size_t nSize = edges.size();

    JNodePtr currnode = edges[0]->getNodeAt(1);

    for( size_t i = 1; i < nSize; i++) {
        JNodePtr v0 = edges[i]->getNodeAt(0);
        JNodePtr v1 = edges[i]->getNodeAt(1);
        if( v1 == currnode)  swap(v0,v1);
        if( v0 != currnode) {
            cout << "Warning: The edge sequence is not a chain" << endl;
            /*
                        for( size_t i = 1; i < nSize; i++) {
                             const JNodePtr &v0 = edges[i]->getNodeAt(0);
                             const JNodePtr &v1 = edges[i]->getNodeAt(1);
                             cout << v0->getID() << "  " << v1->getID() << endl;
                        }
            */
            bndnodes.clear();
            return 1;
        }
        bndnodes.push_back(v0);
        currnode = v1;
    }

    if( currnode != start_node) bndnodes.push_back(currnode);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JEdge :: remove_unattached_lower_entities()
{
    if( this->getStatus() != JMeshEntity::REMOVE) return 1;

    if( nodes[0]->getNumHigherRelations() == 0 )
        nodes[0]->setStatus(JMeshEntity::REMOVE);

    if( nodes[1]->getNumHigherRelations() )
        nodes[1]->setStatus(JMeshEntity::REMOVE);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JEdgeGeometry :: makeUniform( const JEdgeSequence &edgeseq)
{
    double curve_length = getLength(edgeseq);

    JNodeSequence nodes;
    JEdgeTopology::getChainNodes(edgeseq, nodes);

    size_t numedges  = edgeseq.size();
    size_t numnodes  = nodes.size();

    vector<double> param;
    bool closed_curve = JEdgeTopology::isClosed(edgeseq);

    double dt = 1.0/(double)numedges;
    vector<Point3D> newCoords( numnodes );

    double sum  = 0.0;
    if( closed_curve ) {
        param.resize( numnodes + 1);
        for( size_t i = 0; i < numnodes; i++) {
            param[i] = sum/curve_length;
            sum = sum + JEdgeGeometry::getLength( edgeseq[i] );
        }
        param[numnodes] = 1.0;

        for( size_t i = 1; i < numnodes; i++) {
            double t = i*dt;
            for( size_t j = 0; j < numnodes+1; j++) {
                if( t >= param[j] && t < param[j+1] ) {
                    Point3D p0 = nodes[j]->getXYZCoords();
                    Point3D p1 = nodes[(j+1)%numnodes]->getXYZCoords();
                    double  s  = (t - param[j])/(param[j+1]-param[j]);
                    newCoords[i][0] = (1-s)*p0[0] + s*p1[0];
                    newCoords[i][1] = (1-s)*p0[1] + s*p1[1];
                    newCoords[i][2] = (1-s)*p0[2] + s*p1[2];
                    break;
                }
            }
        }
        for( size_t i = 1; i < numnodes; i++)
            nodes[i]->setXYZCoords(newCoords[i]);
        return 0;
    }

    param.resize( numnodes );
    for( size_t i = 0; i < numnodes-1; i++) {
        param[i] = sum/curve_length;
        sum = sum + JEdgeGeometry::getLength( edgeseq[i] );
    }
    param[numnodes-1] = 1.0;

    for( size_t i = 1; i < numnodes-1; i++) {
        double t = i*dt;
        for( size_t j = 0; j < numnodes; j++) {
            if( t >= param[j] && t < param[j+1] ) {
                Point3D p0 = nodes[j]->getXYZCoords();
                Point3D p1 = nodes[(j+1)]->getXYZCoords();
                double  s  = (t - param[j])/(param[j+1]-param[j]);
                newCoords[i][0] = (1-s)*p0[0] + s*p1[0];
                newCoords[i][1] = (1-s)*p0[1] + s*p1[1];
                newCoords[i][2] = (1-s)*p0[2] + s*p1[2];
                break;
            }
        }
    }

    for( size_t i = 1; i < numnodes-1; i++)
        nodes[i]->setXYZCoords(newCoords[i]);

    return 0;
}


////////////////////////////////////////////////////////////////////////////////
