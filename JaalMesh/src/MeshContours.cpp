#include "MeshContour.hpp"

/*
bool JMeshContour :: isLeft( const Point2D &p0, const Point2D &p1, const Point2D &p2)
{
//  Function : test if a point is Left|On|Right of an infinite 2D line.
//  Input:  three points P0, P1, and P2
//  Return: >0 for P2 left of the line through P0 to P1
//          =0 for P2 on the line
//          <0 for P2 right of the line

 return ( (P1[0] - P0[0]) * (P2[1] - P0[1])
           - (P2[0] - P0[0]) * (P1[1] - P0[1]) );
}
*/
//////////////////////////////////////////////////////////////////////////////////
int
JMeshContour::setMesh( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return 1;

    if( mesh->getTopology()->getDimension() == 1 ) {
        edges = mesh->getEdges();
        setSegments(edges);
    }

    if( mesh->getTopology()->getDimension() == 2 ) {
        mesh->getTopology()->searchBoundary();
        mesh->getTopology()->getBoundary(edges);
        setSegments(edges);
    }

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////
int JMeshContour:: setSegments( const JEdgeSequence &e)
{
    edges = e;
    getChain();
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////
int
JMeshContour::getNumComponents() const
{
    for( const JEdgePtr &edge: edges)
        edge->setVisitBit(0);

    int nCount = 0;
    JNodePtr   search_node;
    JEdgePtr   start_edge;

    /*
    while(1) {
    start_edge.reset();
    for( const JEdgePtr &edge: edges) {
        if( edge->getVisitBit() == 0) {
            start_edge = edge;
            break;
        }
    }
    if( start_edge == nullptr) break;
    nCount++;

           search_node = start_edge->getNodeAt(1);
           start_edge->setVisitBit(1);
           while(1) {
           for( const JEdgePtr &edge: edges) {
                if( edge->getVisitBit() == 0) {
                    if( edge->getNodeAt(0) == search_node)
                        search_node = edge->getNodeAt(1);
                        edge->setVisitBit(1);
                        break;
                    }
                    if( edge->getNodeAt(1) == search_node)
                        search_node = edge->getNodeAt(0);
                        edge->setVisitBit(1);
                        break;
                    }
             }
    }
    */
    return nCount;
}

//////////////////////////////////////////////////////////////////////////////////

int JMeshContour::getEdgeNumber( const JNodePtr &vtx)
{
    size_t numSegments =  edges.size();

    for( size_t i = 0; i < numSegments; i++)
        if( edges[i]->getNodeAt(0) == vtx) return i;
    return -1;
}

//////////////////////////////////////////////////////////////////////////////////

double  JMeshContour :: getCurvedDistance( const JNodePtr &v1, const JNodePtr &v2)
{
    int pos0 = getEdgeNumber(v1);
    int pos1 = getEdgeNumber(v2)-1;

    double sum = 0.0;
    for( int i = pos0; i <= pos1; i++) {
        const JNodePtr &v0 = edges[i]->getNodeAt(0);
        const JNodePtr &v1 = edges[i]->getNodeAt(1);
        sum += JNodeGeometry::getLength(v0,v1);
    }
    return sum;
}

//////////////////////////////////////////////////////////////////////////////////

double  JMeshContour :: getDilation( const JNodePtr &v1, const JNodePtr &v2)
{
    double curvedDistance  = getCurvedDistance(v1,v2);
    double shortestDistance = JNodeGeometry::getLength(v1,v2);

    return curvedDistance/shortestDistance;
}

//////////////////////////////////////////////////////////////////////////////////
JNodePtr  JMeshContour :: getNearestNode( double t)
{
    cout << "NOt implemented " << endl;
    return nullptr;
}
//////////////////////////////////////////////////////////////////////////////////
JNodePtr  JMeshContour :: getMidNode( const JNodePtr &v1, const JNodePtr &v2)
{
    int pos0 = getEdgeNumber(v1);
    int pos1 = getEdgeNumber(v2)-1;

    double sum = 0.0;
    for( int i = pos0; i <= pos1; i++) {
        const JNodePtr &v0 = edges[i]->getNodeAt(0);
        const JNodePtr &v1 = edges[i]->getNodeAt(1);
        sum += JNodeGeometry::getLength(v0,v1);
    }

    double half_distance= 0.5*sum;
    sum = 0.0;
    for( int i = pos0; i <= pos1; i++) {
        const JNodePtr &v0 = edges[i]->getNodeAt(0);
        const JNodePtr &v1 = edges[i]->getNodeAt(1);
        sum += JNodeGeometry::getLength(v0,v1);
        if( sum >= half_distance) return v0;
    }

    return nullptr;
}

//////////////////////////////////////////////////////////////////////////////////

JNodeSequence
JMeshContour::getCorners( double theta )
{
    cornerAngle = theta;
    JNodeSequence corners;
    int numSegments = edges.size();

    Vec3D avec , bvec;
    for( int i = 0; i < numSegments; i++) {
        const JNodePtr  &v0 = edges[i]->getNodeAt(0);
        const JNodePtr  &v1 = edges[(i+1)%numSegments]->getNodeAt(0);
        const JNodePtr  &v2 = edges[(i+2)%numSegments]->getNodeAt(0);
        const Point3D   &p0 = v0->getXYZCoords();
        const Point3D   &p1 = v1->getXYZCoords();
        const Point3D   &p2 = v2->getXYZCoords();
        JMath::make_vector( p0, p1, avec);
        JMath::make_vector( p2, p1, bvec);
        double angle = 180.0-JMath::getVecAngle( bvec, avec, ANGLE_IN_DEGREES);
        if( angle >  cornerAngle ) {
            cout << "Corner" << v1->getID() << " Angle " << angle << endl;
            cout << "Neigh " << v0->getID() << "  " << v2->getID() << endl;
            cout << endl;
            corners.push_back(v1);
        }
    }

    return corners;
}
//////////////////////////////////////////////////////////////////////////////////
bool
JMeshContour::isChain() const
{
    size_t numSegments = edges.size();
    for( size_t i = 0; i < numSegments-1; i++) {
        const JNodePtr &v1 = edges[i]->getNodeAt(1);
        const JNodePtr &v0 = edges[(i+1)%numSegments]->getNodeAt(0);
        if( v0 != v1) return 0;
    }
    return 1;
}
/////////////////////////////////////////////////////////////////////////////////
int
JMeshContour::splitAt( const JNodePtr &midnode, JEdgeSequence &seq1, JEdgeSequence &seq2)
{
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////

JEdgeSequence
JMeshContour ::getChain()
{
    // Starting from the first node of the first edge,
    // link the edge segments to make a chain.

    size_t nSize = edges.size();
    if( nSize == 0) return edges;

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
    return edges;
}

/////////////////////////////////////////////////////////////////////////////////

JEdgeSequence
JMeshContour::getChain(const JEdgePtr &start_edge)
{
    // starting from the "given" edge, make a chain.

    JEdgeSequence::iterator it;
    it = find(edges.begin(), edges.end(), start_edge);
    if( it == edges.end() ) return edges;

    std::rotate( edges.begin(), it, edges.end() );

    assert( edges[0] == start_edge );
    JEdgeTopology::getChain( edges);
    return edges;
}

///////////////////////////////////////////////////////////////////////////////

JEdgeSequence
JMeshContour ::getChain(const JNodePtr &start_node)
{
    // Starting from the given vertex, build the chain...
    size_t numSegments = edges.size();

    JEdgePtr start_edge = nullptr;
    for( size_t i = 0; i < numSegments; i++) {
        if( edges[i]->getNodeAt(0) == start_node) {
            start_edge = edges[i];
            break;
        }
    }

    if( start_edge == nullptr) {
        for( size_t i = 0; i < numSegments; i++) {
            if( edges[i]->getNodeAt(1) == start_node) {
                start_edge = edges[i];
                start_edge->reverse();
                break;
            }
        }
    }

    if( start_edge == nullptr ) {
        JEdgeSequence empty;
        return empty;
    }

    JEdgeTopology::getChain(edges, start_edge);
    return edges;
}

///////////////////////////////////////////////////////////////////////////////

int
JMeshContour::getOrientation() const
{
    // First find the rightmost lowest vertex of the polygon ..

    return 0;
}


/////////////////////////////////////////////////////////////////////////////////
int
JMeshContour ::isCloseable()
{
    // If all the nodes have atleast two inContour, then it is probably closeable.
    // In general, this could be more complicated, but not for our cases, at
    // present...
    std::map<JNodePtr, JEdgeSequence> relations00;

    size_t numSegments = edges.size();
    for (size_t i = 0; i < numSegments; i++) {
        const JNodePtr &v0 = edges[i]->getNodeAt(0);
        const JNodePtr &v1 = edges[i]->getNodeAt(1);
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
JMeshContour ::isTopologicalSimple()
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
JMeshContour ::isClosed()
{
    const JNodePtr &vf = edges.front()->getNodeAt(0);
    const JNodePtr &vl = edges.back()->getNodeAt(1);
    if (vf == vl) return 1;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JNodeSequence
JMeshContour::getChainNodes()
{
    if( !isChain() ) getChain();

    size_t numSegments = edges.size();

    JNodeSequence bndnodes( numSegments);
    for( size_t i = 0; i < numSegments; i++)
        bndnodes[i] = edges[i]->getNodeAt(0);

    return bndnodes;
}

///////////////////////////////////////////////////////////////////////////////

double  JMeshContour :: getLength() const
{
    double sum = 0.0;
    for( const JEdgePtr &edge : edges)
        sum += JNodeGeometry::getLength(edge->getNodeAt(0), edge->getNodeAt(1) );
    return sum;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshContour :: smooth(int numIter)
{
    JNodeSequence nodes = getChainNodes();
    size_t numnodes  = nodes.size();

    vector<Point3D>  points(numnodes);
    for( size_t i = 0; i < numnodes; i++)
        points[i] = nodes[i]->getXYZCoords();

    vector<Point3D> newPos(numnodes);
    Point3D xyz;
    for( int iter = 0; iter < numIter; iter++) {
        for( size_t i = 0; i < numnodes; i++) {
            xyz[0] = 0.5*(points[i][0] + points[(i+2)%numnodes][0]);
            xyz[1] = 0.5*(points[i][1] + points[(i+2)%numnodes][1]);
            xyz[2] = 0.5*(points[i][2] + points[(i+2)%numnodes][2]);
            newPos[(i+1)%numnodes][0] = xyz[0];
            newPos[(i+1)%numnodes][1] = xyz[1];
            newPos[(i+1)%numnodes][2] = xyz[2];
        }
        points = newPos;
    }

    for( size_t i = 0; i < numnodes; i++)
        nodes[i]->setXYZCoords( points[i] );
}
///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JMeshContour :: getDiscretized( JEdgeSequence &subedges, int N)
{
    /*
        double len = getLength();

        bool   closed = isClosed();

        JNodeSequence nodes = getChainNodes();

        int nSize = nodes.size();
        vector<double>  t(nSize);

        double sumlen = 0;
        int    index  = 0;
        for( const JEdgePtr &edge: edges) {
            t[index++] = sumlen/len;
            sumlen += JEdgeGeometry::getLength(edge);
        }

        if(closed) {
           t.push_back(1.0);
           nodes.push_back( nodes[0] );
        }

        JNodeSequence newnodes(N);
        Point3D xyz;

        double dt = 1.0/(double)N;
        for( size_t i = 0; i < N; i++) {
            double u = i*dt;
            bool found = 0;
            for( size_t j = 0; j < nSize; j++) {
                double t1 = t[j];
                double t2 = t[j+1];
                if( u >= t1 && u <= t2 ) {
                    const Point3D &p1 = nodes[j]->getXYZCoords();
                    const Point3D &p2 = nodes[j+1]->getXYZCoords();
                    double tn = (u-t1)/(t2-t1);
                    xyz[0] = (1-tn)*p1[0] + tn*p2[0];
                    xyz[1] = (1-tn)*p1[1] + tn*p2[1];
                    xyz[2] = (1-tn)*p1[2] + tn*p2[2];
                    newnodes[i] = JNode::newObject();
                    newnodes[i]->setXYZCoords(xyz);
                    found = 1;
                    break;
                }
            }
            assert( found );
        }
        newnodes.push_back(newnodes[0]);
        JEdgeSequence edges(N);
        for( size_t i = 0; i < N; i++)
            edges[i] = JEdge::newObject( newnodes[i], newnodes[(i+1)%N]);
        return edges;
    */
}

////////////////////////////////////////////////////////////////////////////////


JEdgeSequence JMeshContour :: getDiscretized( int N, int method )
{
    /*
        double len = getLength();
        bool   closed = isClosed();

        JNodeSequence nodes = getChainNodes();
        int nSize = nodes.size();
        vector<double>  t(nSize);

        double sumlen = 0;
        int    index  = 0;
        for( const JEdgePtr &edge: edges) {
            t[index++] = sumlen/len;
            sumlen += JEdgeGeometry::getLength(edge);
        }

        if(closed) {
           t.push_back(1.0);
           nodes.push_back( nodes[0] );
        }

        vector<JEdgeSequence> subedges;
        int numSegments = 0;
        for( const JNodePtr &vtx : bnodes)
             if( vtx->hasAttribute("Constraint") ) nSegments++;

        if( numSegments ) {
        for( const JEdgePtr &edge: edges) {
        }


        JNodeSequence newnodes(N);
        Point3D xyz;

        double dt = 1.0/(double)N;
        for( size_t i = 0; i < N; i++) {
            double u = i*dt;
            bool found = 0;
            for( size_t j = 0; j < nSize; j++) {
                double t1 = t[j];
                double t2 = t[j+1];
                if( u >= t1 && u <= t2 ) {
                    const Point3D &p1 = nodes[j]->getXYZCoords();
                    const Point3D &p2 = nodes[j+1]->getXYZCoords();
                    double tn = (u-t1)/(t2-t1);
                    xyz[0] = (1-tn)*p1[0] + tn*p2[0];
                    xyz[1] = (1-tn)*p1[1] + tn*p2[1];
                    xyz[2] = (1-tn)*p1[2] + tn*p2[2];
                    newnodes[i] = JNode::newObject();
                    newnodes[i]->setXYZCoords(xyz);
                    found = 1;
                    break;
                }
            }
            assert( found );
        }
        newnodes.push_back(newnodes[0]);
        JEdgeSequence edges(N);
        for( size_t i = 0; i < N; i++)
            edges[i] = JEdge::newObject( newnodes[i], newnodes[(i+1)%N]);
        return edges;
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshContour :: reduceDilation()
{
    /*
        JNodeSequence nodes = getCorners(cornerAngle);
        deque<NodePair>  edgeQ;

        int numCorners = nodes.size();
        for( int i = 0; i < nodes.size(); i++)
            edgeQ.push_back(make_pair(nodes[i], nodes[(i+1)%numCorners]));

        simplified.clear();
        while(!edgeQ.empty()) {
            NodePair nodepair = edgeQ.front();
            edgeQ.pop_front();
            JNodePtr v0 = nodepair.first;
            JNodePtr v1 = nodepair.second;
            double   df = getDilation(v0,v1);
            if( df > dilationFactor) {
                JNodePtr vm = getMidNode(v0,v1);
                edgeQ.push_front(make_pair(vm,v1));
                edgeQ.push_front(make_pair(v0,vm));
            } else
                simplified.push_back(v0);
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

/*
JMeshPtr JImageContour :: getContours( const string &filename)
{
  cv::Mat src = cv::imread( filename.c_str(), 1 );

  cv::Mat canny_output;

  int thresh = 100;
  cv::Canny( src, canny_output, thresh, thresh*2, 3 );

  vector<vector<  cv::Point> > contours;
  vector<cv::Vec4i> hierarchy;
  cv::findContours( canny_output, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

  int numContours = contours.size();
  Point3D xyz;
  xyz[0] = 0.0;
  xyz[1] = 0.0;
  xyz[2] = 0.0;
  JMeshPtr mesh = JMesh::newObject();

  JNodeSequence  nodes;
  JEdgeSequence edges;

  for( int i = 0; i < contours.size(); i++) {
       int npoints = contours[i].size();
       nodes.resize(npoints);
       for( int j = 0; j < contours[i].size(); j++) {
            JNodePtr vtx = JNode::newObject();
            xyz[0] = contours[i][j].x;
            xyz[1] = contours[i][j].y;
            vtx->setXYZCoords(xyz);
            nodes[i] = vtx;
       }
       mesh->addObjects( nodes );
       edges.resize( npoints);
       for( size_t j = 0; j < npoints; j++)
            edges[j]  = JEdge::newObject(nodes[j], nodes[(j+1)%npoints]);
       mesh->addObjects( edges );

  }
  return mesh;
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
vector<double>
JMeshContour::getU()
{
   double total_length = getLength();

   vector<double>  u;
}
*/
