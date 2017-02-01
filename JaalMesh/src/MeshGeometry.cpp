#include "MeshGeometry.hpp"

using namespace std;
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

struct SpatialSort {
    bool operator() ( const JNodePtr &va, const JNodePtr &vb)
    {
        const Point3D &pa = va->getXYZCoords();
        const Point3D &pb = vb->getXYZCoords();
        if( pa[0] < pb[0] ) return 1;
        if( pa[1] < pb[1] ) return 1;
        if( pa[2] < pb[2] ) return 1;
        return 0;
    }
};

///////////////////////////////////////////////////////////////////////////////

JLogger*  JMeshGeometry::logger = JLogger::getInstance();

///////////////////////////////////////////////////////////////////////////////

vector<double> JMeshGeometry :: getEuclideanDistance( const JMeshPtr &mesh1, const JMeshPtr &mesh2)
{
    size_t numnodes = mesh1->getSize(0);

    if( numnodes != mesh2->getSize(0) )
        cout << "Warning: Two meshes are not same " << endl;

    numnodes = min(numnodes, mesh2->getSize(0));

    vector<double> dist(numnodes);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v1 = mesh1->getNodeAt(i);
        const JNodePtr &v2 = mesh2->getNodeAt(i);
        dist[i] = JNodeGeometry::getLength(v1,v2);
    }
    return dist;
}

///////////////////////////////////////////////////////////////////////////////
double JMeshGeometry :: getMeanEdgeLength() const
{
    if( mesh == nullptr) return 0.0;

/*
    JEdgeSequence boundedges;
    mesh->getTopology()->getBoundary(boundedges);
    int numedges = boundedges.size();
    assert( numedges );
    size_t index = 0;
    vector<double>  elen(numedges);
    for( const JEdgePtr &e : boundedges)
        elen[index++] = JEdgeGeometry::getLength(e);
    boost::sort( elen );
    return elen[numedges/2];
*/

    int numedges = mesh->getSize(1);
    assert( numedges );
    size_t index = 0;
    vector<double>  elen;
    elen.reserve(numedges);

    for( size_t i = 0; i < numedges; i++)  {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) 
            elen.push_back( JEdgeGeometry::getLength(edge));
    }
       
    boost::sort( elen );
    return elen[elen.size()/2];
}
///////////////////////////////////////////////////////////////////////////////

double JMeshGeometry :: getMaxDistance( const JMeshPtr &mesh1, const JMeshPtr &mesh2)
{
    size_t numnodes = mesh1->getSize(0);
    if( numnodes != mesh2->getSize(0) )
        cout << "Warning: Two meshes are not same " << endl;

    numnodes = min(numnodes, mesh2->getSize(0));
    double maxdist = 0.0;
    #pragma omp parallel for reduction(max:maxdist)
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v1 = mesh1->getNodeAt(i);
        const JNodePtr &v2 = mesh2->getNodeAt(i);
        maxdist = max( maxdist, JNodeGeometry::getLength(v1,v2));
    }
    return maxdist;
}
///////////////////////////////////////////////////////////////////////////////

int JMeshGeometry :: getDimension() const
{
    if( mesh == nullptr ) return 0;

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive()) {
            if( !vtx->isBoundary() ) {
                const Point3D &xyz = vtx->getXYZCoords();
                if( fabs(xyz[2]) > 0.0) return 3;
            }
        }
    }

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive()) {
            if( !vtx->isBoundary() ) {
                const Point3D &xyz = vtx->getXYZCoords();
                if( fabs(xyz[1]) > 0.0) return 2;
            }
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////
bool JMeshGeometry :: hasPlanarEmbedding() const 
{
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive()) {
            Point3D xyz = vtx->getXYZCoords();
            if( fabs(xyz[2]) > 1.0E-10) return 0;
        }
    }
    return 1;
}
///////////////////////////////////////////////////////////////////////////////
void JMeshGeometry :: addNoise( double maxVal)
{
    if( mesh == nullptr ) return;

    double minVal = 0.0;
    size_t numnodes = mesh->getSize(0);

    #pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive()) {
            if( !vtx->isBoundary() ) {
                Point3D xyz = vtx->getXYZCoords();
                xyz[0] += JMath::random_value(minVal, maxVal);
                xyz[1] += JMath::random_value(minVal, maxVal);
                xyz[2] += JMath::random_value(minVal, maxVal);
                vtx->setXYZCoords(xyz);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

int JMeshGeometry :: spatialSort()
{
    if( mesh == nullptr ) return 1;
    sort( mesh->nodes.begin(), mesh->nodes.end(), SpatialSort() );
    size_t numnodes = mesh->getSize(0);
    int index = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) v->setID( index++);
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeometry :: setBoundaryConstraints(bool preserve_boundary )
{
    if( mesh == nullptr ) return;

    JNodeSequence nodes;
    mesh->getTopology()->getBoundary(nodes);

    int gid = 1;
    if( preserve_boundary ) {
        for( const JNodePtr &vtx : nodes) vtx->setAttribute("Constraint", gid);
    } else {
        for( const JNodePtr &vtx : nodes) vtx->deleteAttribute("Constraint");
    }
}


///////////////////////////////////////////////////////////////////////////////

JBoundingBox
JMeshGeometry ::getBoundingBox() const
{
    logger->setInfo("Calculating the mesh bounding box " );

    JBoundingBox box;

    size_t numnodes = mesh->getSize(0);

    if( numnodes < 2) return box;

    double xmin, xmax, ymin, ymax, zmin, zmax;
    Point3D xyz;
    xyz = mesh->getNodeAt(0)->getXYZCoords();

    xmin = xyz[0];
    xmax = xyz[0];
    ymin = xyz[1];
    ymax = xyz[1];
    zmin = xyz[2];
    zmax = xyz[2];

//#pragma omp for private(xyz) reduction(min:xmin,ymin,zmin) reduction(max:xmax,ymax,zmax)
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) {
            xyz = v->getXYZCoords();
            xmin = min(xmin, xyz[0]);
            xmax = max(xmax, xyz[0]);
            ymin = min(ymin, xyz[1]);
            ymax = max(ymax, xyz[1]);
            zmin = min(zmin, xyz[2]);
            zmax = max(zmax, xyz[2]);
        }
    }

    xyz[0] = xmin;
    xyz[1] = ymin;
    xyz[2] = zmin;
    box.setLower(xyz);

    xyz[0] = xmax;
    xyz[1] = ymax;
    xyz[2] = zmax;
    box.setUpper(xyz);

    return box;
}

///////////////////////////////////////////////////////////////////////////////

Point3D JMeshGeometry :: getCenter() const
{
    JBoundingBox box =  getBoundingBox();
    return box.getCenter();
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshGeometry::explode( double alpha)
{
    logger->setInfo("Exploding mesh" );

    JMeshPtr msh = JMesh::newObject();

    int topDim = mesh->getTopology()->getDimension();

    if( topDim == 2 ) {
        size_t numfaces = mesh->getSize(2);
        for( size_t i  = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            JFacePtr fe = f->explode( alpha );
            msh->addObjects( fe->getNodes() );
            msh->addObject(fe);
        }
    }

    if( topDim == 3 ) {
        size_t numcells = mesh->getSize(3);
        for( size_t i  = 0; i < numcells; i++) {
            const JCellPtr &c  = mesh->getCellAt(i);
            JCellPtr ce = c->explode( alpha );
            msh->addObjects( ce->getNodes() );
            msh->addObject(ce);
        }
    }

    return msh;
}

////////////////////////////////////////////////////////////////////////////////

size_t
JMeshGeometry ::count_concave_faces() const
{
    logger->setInfo("Counting concave faces" );

    size_t numfaces = mesh->getSize(2);
    size_t ncount = 0;
    #pragma omp parallel for reduction(+:ncount)
    for (size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        if(f->isActive()  &&  !JFaceGeometry::isConvex(f) ) ncount++;
    }
    return ncount;
}
////////////////////////////////////////////////////////////////////////////////

size_t
JMeshGeometry ::getNumOfInvertedElements() const
{
    logger->setInfo("Counting inverted faces" );


    int topDim = mesh->getTopology()->getDimension();

    double val;

    JMeshQuality mq;
    size_t nCount = 0;
    if( topDim == 2 ) {
        size_t numfaces = mesh->getSize(2);
        nCount = 0;
        #pragma omp parallel private(mq) private(val) reduction(+:nCount)
        mq.setMesh(mesh);
        #pragma omp for
        for (size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                val = mq.getJacobian(face);
                if( val < 0.0) nCount++;
            }
        }
    }

    if( topDim == 3 ) {
        nCount = 0;
        size_t numCells = mesh->getSize(3);
        #pragma omp parallel private(mq) private(val) reduction(+:nCount)
        mq.setMesh(mesh);
        #pragma omp for
        for (size_t i = 0; i < numCells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val = mq.getJacobian(cell);
                if( val < 0.0) nCount++;
            }
        }
    }
    return nCount;
}

////////////////////////////////////////////////////////////////////////////////

int
JMeshGeometry::getCoordsArray( vector<double> &vcoords, vector<size_t> &l2g)
{
    if( mesh == nullptr ) return 1;

    logger->setInfo("Collecting node coordinates");

    size_t numnodes = mesh->getSize(0);

    vcoords.clear();
    l2g.clear();

    vcoords.reserve(3 * numnodes);
    l2g.reserve(numnodes);

    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) {
            l2g.push_back( v->getID() );
            const Point3D &xyz = v->getXYZCoords();
            vcoords.push_back( xyz[0] );
            vcoords.push_back( xyz[1] );
            vcoords.push_back( xyz[2] );
        }
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

int
JMeshGeometry::setCoordsArray(const vector<double> &vcoords, const vector<size_t> &l2g)
{
    logger->setInfo("Setting node coordinates" );

    size_t numnodes = mesh->getSize(0);

    Point3D xyz;
    int index = 0;
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            assert( vtx->getID() == l2g[index] );
            xyz[0] = vcoords[3*index+0];
            xyz[1] = vcoords[3*index+1];
            xyz[2] = vcoords[3*index+2];
            vtx->setXYZCoords(xyz);
            index++;
        }
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

double
JMeshGeometry ::getSurfaceArea() const
{
    logger->setInfo("Calculating surface area");

    double facearea, sumArea = 0.0;

    size_t numfaces = mesh->getSize(2);

    #pragma omp parallel for reduction(+:sumArea) private(facearea)
    for (size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            facearea = JFaceGeometry::getArea(face);
            sumArea += facearea;
        }
    }
    return sumArea;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeometry :: setNormal(const JFacePtr &face)
{
    if( face == nullptr ) return;
    if( !face->isActive() )  return;

    double x[100], y[100], z[100], nx, ny, nz;
    Point3D xyz;
    Vec3F normal;
    int nsize = face->getSize(0);

    for (int j = 0; j < nsize; j++) {
        const JNodePtr &vtx = face->getNodeAt(j);
        xyz = vtx->getXYZCoords();
        x[j] = xyz[0];
        y[j] = xyz[1];
        z[j] = xyz[2];
    }
    PolygonNormal3D(nsize, x, y, z, &nx, &ny, &nz);
    double val2 =  nx*nx + ny*ny + nz*nz;

    if( val2 < 1.0E-06) {
        cout << "Warning: Invalid face normal " << endl;
        normal[0] = nx;
        normal[1] = ny;
        normal[2] = nz;
    } else {
        double multby = 1.0/sqrt( val2 );

        normal[0] = nx*multby;
        normal[1] = ny*multby;
        normal[2] = nz*multby;

    }
    face->setAttribute("Normal", normal);
}

//////////////////////////////////////////////////////////////////////////////

void JMeshGeometry :: setNormal( const JNodePtr &vertex, int weight)
{
    if( !vertex->isActive() ) return ;

    Vec3F normal;
    JFaceSequence vfaces;
    JNode::getRelations(vertex, vfaces);
    if( vfaces.empty() ) {
        cout << "Warning: The vertex " << vertex->getID() << " is isolated for the normal calculation" << endl;
        normal[0] = 0.0;
        normal[1] = 0.0;
        normal[2] = 1.0;
        vertex->setAttribute("Normal", normal);
        return;
    }

    double nx = 0.0;
    double ny = 0.0;
    double nz = 0.0;
    double wsum = 0.0;
    double wj   = 1.0;

    for( size_t j = 0; j < vfaces.size(); j++) {
        if( !vfaces[j]->hasAttribute("Normal") )
            setNormal( vfaces[j] );
        vfaces[j]->getAttribute("Normal", normal);

        if( weight == NORMAL_AREA_WEIGHT )
            wj = fabs(JFaceGeometry::getArea( vfaces[j]));
        wsum += wj;
        nx += wj*normal[0];
        ny += wj*normal[1];
        nz += wj*normal[2];
    }

    if( weight == NORMAL_AREA_WEIGHT ) {
        cout << "Debug CSV Normal Area weight " << endl;
        nx /= wsum;
        ny /= wsum;
        nz /= wsum;
    }
    double val2 =  nx*nx + ny*ny + nz*nz;
    if( val2 < 1.0E-06) {
        cout << "Warning: Invalid node normal " << endl;
        normal[0] = nx;
        normal[1] = ny;
        normal[2] = nz;
    }
    else {
        double multby = 1.0/sqrt( val2 );
        normal[0] = nx*multby;
        normal[1] = ny*multby;
        normal[2] = nz*multby;
    }
    vertex->setAttribute("Normal", normal);
    return;
}

///////////////////////////////////////////////////////////////////////////////
double JFaceGeometry :: getAngleAt( const JFacePtr &face, int pos, int measure)
{
    int n      =  face->getSize(0) ;
    JNodePtr v0 =  face->getNodeAt(pos);
    JNodePtr v1 =  face->getNodeAt(pos+1);
    JNodePtr vn =  face->getNodeAt(pos+n-1);
    double theta;
    theta =  JMath::getTriAngle( v0->getXYZCoords(), v1->getXYZCoords(), vn->getXYZCoords(), measure );
    return theta;
}

///////////////////////////////////////////////////////////////////////////////
double JNodeGeometry :: getAngleDefect( const JNodePtr &vertex)
{
    if( !vertex->isActive() ) return 0;

    JFaceSequence vfaces;
    JNode::getRelations(vertex, vfaces);
    double sumangle = 0.0;
    for( size_t j = 0; j < vfaces.size(); j++) {
        int pos = vfaces[j]->getPosOf(vertex);
        sumangle += JFaceGeometry::getAngleAt(vfaces[j], pos, ANGLE_IN_DEGREES);
    }

    double defect = 360.0 - sumangle;

    return defect;
}
///////////////////////////////////////////////////////////////////////////////

void
JMeshGeometry ::setFacesNormal()
{
    logger->setInfo("Calculating face normals ");

    size_t nSize = mesh->getSize(2);

    #pragma omp parallel for
    for (size_t i = 0; i < nSize; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        setNormal(face);
    }
}

///////////////////////////////////////////////////////////////////////////////

size_t
JMeshGeometry ::setAngleDefects(double cutoff)
{
    logger->setInfo("Calculating angle defects ");

    mesh->deleteNodeAttribute("AngleDefect");

    int relexist2 = mesh->buildRelations(0,2);

    if( !mesh->getTopology()->isConsistent() )
        mesh->getTopology()->getConsistent();

    size_t ncount = 0;
    size_t nSize = mesh->getSize(0);
    for (size_t i = 0; i < nSize; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            float angle = JNodeGeometry::getAngleDefect(vertex);
            if( fabs(angle) >= cutoff ) {
                vertex->setAttribute("AngleDefect", angle);
                ncount++;
            }
        }
    }

    if( !relexist2 )  mesh->clearRelations(0,2);
    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshGeometry ::setNodesNormal(int weight)
{
    if( mesh == nullptr ) return;
    // You need surface to calculate the normals ...
    if( mesh->getTopology()->getDimension() < 2) return;

    logger->setInfo("Calculating nodes normal ");

    mesh->buildRelations(0,2);

    size_t nSize = mesh->getSize(0);
    #pragma omp parallel for
    for (size_t i = 0; i < nSize; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        setNormal(vertex, weight);
    }
}

///////////////////////////////////////////////////////////////////////////////

size_t
JMeshGeometry :: setSharpEdges( double creaseAngle)
{
    if( mesh == nullptr ) return 0;

    logger->setInfo("Calculating sharp edges in the model");

    setFacesNormal();

    size_t numEdges = mesh->getSize(1);

    Vec3F fn1, fn2;

    mesh->buildRelations(1,2);

    mesh->deleteEdgeAttribute("CreaseAngle");

    int err = 0;
    int ncount = 0;

    JFaceSequence efaces;
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if(edge->isActive()) {
            JEdge::getRelations(edge, efaces );
            if( efaces.size() == 2 ) {
                err = efaces[0]->getAttribute( "Normal", fn1);
                assert( err == 0);
                err = efaces[1]->getAttribute( "Normal", fn2);
                assert( err == 0);
                float  angle = JMath::getVecAngle(fn1, fn2, ANGLE_IN_DEGREES);
                if (angle <= 90 && angle >= creaseAngle) {
                    edge->setAttribute("CreaseAngle", angle);
                    ncount++;
                } else if (angle >= 90 && fabs(180 - angle) >= creaseAngle) {
                    edge->setAttribute("CreaseAngle", angle);
                    ncount++;
                }
            }
        }
    }
    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

double JNodeGeometry ::getFeatureLength(const JNodePtr &vertex)
{
    /*
        if (!vertex->isBoundary()) return MAXDOUBLE;

        JNodeSequence vneighs;
        JNode::getRelations( vertex, vneighs );

        assert(!vneighs.empty());

        double minlen = MAXDOUBLE;

        int  nSize = vneighs.size();
        for (int j = 0; j < nSize; j++) {
            if (vneighs[j]->isBoundary()) {
                const Point3D &p0 = vertex->getXYZCoords();
                const Point3D &p1 = vneighs[j]->getXYZCoords();
                minlen = min(minlen, JMath::length(p0, p1));
            }
        }
        return minlen;
    */
    cout << "Not implemendted " << endl;
    exit(0);
}

///////////////////////////////////////////////////////////////////////////////
JNodePtr JMeshGeometry :: getNearest( const Point3D &query)
{
    JNodePtr minNode = nullptr;
    if( mesh == nullptr ) return nullptr;

    double minval = std::numeric_limits<double>::max();
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        const Point3D &p3d = vtx->getXYZCoords();
        double dist = JMath::length2( p3d, query);
        if( dist < minval ) {
            minNode = vtx;
            minval  = dist;
        }
    }
    return minNode;

}
///////////////////////////////////////////////////////////////////////////////

JSphere JMeshGeometry :: getMinSphere()
{
    JSphere sph;
    logger->setInfo("Calculating the mesh minimum sphere " );

    double seed = 0;
    std::srand (seed);

    std::list<std::vector<double> > lp;
    size_t numnodes = mesh->getSize(0);
    for (size_t i=0; i< numnodes; ++i) {
        std::vector<double> p(3);
        const JNodePtr &vtx = mesh->getNodeAt(i);
        const Point3D &xyz = vtx->getXYZCoords();
        p[0] =  xyz[0];
        p[1] =  xyz[1];
        p[2] =  xyz[2];
        lp.push_back(p);
    }

    // define the types of iterators through the points and their coordinates
    // ----------------------------------------------------------------------
    typedef std::list<std::vector<double> >::const_iterator PointIterator;
    typedef std::vector<double>::const_iterator CoordIterator;

    // create an instance of Miniball
    // ------------------------------
    typedef Miniball::
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator>> MB;

    MB mb (3, lp.begin(), lp.end());

    const double*xyz = mb.center();

    Point3D center;
    center[0] = xyz[0];
    center[1] = xyz[1];
    center[2] = xyz[2];

    sph.setCenter( center );
    sph.setRadius( sqrt(mb.squared_radius()) );

    return sph;
}

////////////////////////////////////////////////////////////////////////////////

JCellPtr JMeshGeometry :: getMinBox( const vector<Point3D> &inpoints)
{
    size_t numpoints = inpoints.size();

    gdiam_real  *points;
    points = (gdiam_point)malloc( sizeof( gdiam_point_t ) * numpoints);
    assert( points != NULL );

    #pragma omp parallel for
    for  ( size_t i = 0; i < numpoints; i++ ) {
        points[ 3*i + 0 ] = inpoints[i][0];
        points[ 3*i + 1 ] = inpoints[i][1];
        points[ 3*i + 2 ] = inpoints[i][2];
    }

    GPointPair  pair;
    pair = gdiam_approx_diam_pair( (gdiam_real *)points, numpoints, 0.0 );

    gdiam_point  *pnt_arr;
    pnt_arr = gdiam_convert( (gdiam_real *)points, numpoints );

    gdiam_bbox  bb;
    bb = gdiam_approx_mvbb( pnt_arr, numpoints, 1.0E-03);

    double coord[3];
    JNodeSequence nodes(8);
    for( int i = 0; i < 8; i++)
        nodes[i] = JNode::newObject();

    bb.getCorner(  0, 0, 0, coord );
    nodes[0]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  1, 0, 0, coord );
    nodes[1]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  1, 1, 0, coord );
    nodes[2]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  0, 1, 0, coord );
    nodes[3]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  0, 0, 1, coord );
    nodes[4]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  1, 0, 1, coord );
    nodes[5]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  1, 1, 1, coord );
    nodes[6]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  0, 1, 1, coord );
    nodes[7]->setXYZCoords( coord[0], coord[1], coord[2] );

    JCellPtr hex = JHexahedron::newObject();
    hex->setNodes( nodes );

    free(pnt_arr);
    return hex;
}

//////////////////////////////////////////////////////////////////////////////

JFacePtr JMeshGeometry :: getMinRectangle( const vector<Point2D> &inpoints)
{
    size_t numpoints = inpoints.size();

    gdiam_real  *points;
    points = (gdiam_point)malloc( sizeof( gdiam_point_t ) * numpoints);
    assert( points != NULL );

    #pragma omp parallel for
    for  ( size_t i = 0; i < numpoints; i++ ) {
        points[ 3*i + 0 ] = inpoints[i][0];
        points[ 3*i + 1 ] = inpoints[i][1];
        points[ 3*i + 2 ] = 0.0;
    }

    GPointPair  pair;
    pair = gdiam_approx_diam_pair( (gdiam_real *)points, numpoints, 0.0 );

    gdiam_point  *pnt_arr;
    pnt_arr = gdiam_convert( (gdiam_real *)points, numpoints );

    gdiam_bbox  bb;
    bb = gdiam_approx_mvbb( pnt_arr, numpoints, 1.0E-03);

    JNodeSequence nodes(4);
    for( int i = 0; i < 4; i++)
        nodes[i] = JNode::newObject();

    double coord[3];
    bb.getCorner(  0, 0, 0, coord );
    nodes[0]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  1, 0, 0, coord );
    nodes[1]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  1, 1, 0, coord );
    nodes[2]->setXYZCoords( coord[0], coord[1], coord[2] );

    bb.getCorner(  0, 1, 0, coord );
    nodes[3]->setXYZCoords( coord[0], coord[1], coord[2] );

    JFacePtr quad = JQuadrilateral::newObject();
    quad->setNodes( nodes );

    free(pnt_arr);
    return quad;
}
//////////////////////////////////////////////////////////////////////////////
JFacePtr JMeshGeometry :: getMinRectangle( const JNodeSequence &uvNodes)
{
    size_t numPoints = uvNodes.size();
    if( numPoints < 2 ) return nullptr;

    vector<Point2D> inPoints(numPoints);
    Point3D xyz;
    size_t index = 0;
    for( const JNodePtr &vtx : uvNodes) {
         xyz = vtx->getXYZCoords();
         inPoints[index][0] = xyz[0];
         inPoints[index][1] = xyz[1];
         index++;
    }
    return getMinRectangle( inPoints );
}
//////////////////////////////////////////////////////////////////////////////

JCellPtr JMeshGeometry :: getMinBox( const JNodeSequence &nodes)
{
    if( nodes.empty() ) return nullptr;

    vector<Point3D> points;
    size_t numpoints = nodes.size();
    points.reserve( numpoints);

    for( const JNodePtr &vtx : nodes)
        if( vtx->isActive() ) points.push_back( vtx->getXYZCoords() );

    return getMinBox( points);
}

///////////////////////////////////////////////////////////////////////////////

JCellPtr JMeshGeometry :: getMinBox()
{
    if( mesh == nullptr) return nullptr;

    JNodeSequence nodes = mesh->getNodes();
    return getMinBox( nodes);
}

///////////////////////////////////////////////////////////////////////////////

JNodeSequence JMeshGeometry :: getConvexCorners(double  cutoff_angle) const
{
    double  convexAngle  = 180 - cutoff_angle;

    JNodeSequence nodes;
    JFaceSequence faces;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0;  i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() && vtx->isBoundary() ) {
            JNode::getRelations(vtx, faces);
            double sum = 0.0;
            for( const JFacePtr &face: faces)
                sum = sum + JFaceGeometry::getAngleAt(face,vtx, ANGLE_IN_DEGREES);
            if( sum < convexAngle ) nodes.push_back(vtx);
        }
    }
    return nodes;
}

///////////////////////////////////////////////////////////////////////////////

JNodeSequence JMeshGeometry :: getConcaveCorners(double  cutoff_angle) const
{
//    double  convexAngle  = 180 - cutoff_angle;
    double  concaveAngle = 180 + cutoff_angle;

    JNodeSequence nodes;
    JFaceSequence faces;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0;  i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() && vtx->isBoundary() ) {
            JNode::getRelations(vtx, faces);
            double sum = 0.0;
            for( const JFacePtr &face: faces)
                sum = sum + JFaceGeometry::getAngleAt(face,vtx, ANGLE_IN_DEGREES);
            if( sum > concaveAngle) nodes.push_back(vtx);
        }
    }
    return nodes;
}
///////////////////////////////////////////////////////////////////////////////
JNodeSequence JMeshGeometry :: getBoundaryCorners(double  cutoff_angle) const
{
    double  convexAngle  = 180 - cutoff_angle;
    double  concaveAngle = 180 + cutoff_angle;

    mesh->deleteNodeAttribute("CornerAngle");

    mesh->getTopology()->searchBoundary();
    mesh->buildRelations(0,2);

    JNodeSequence nodes;
    JFaceSequence faces;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0;  i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() && vtx->isBoundary() ) {
            JNode::getRelations(vtx, faces);
            double sum = 0.0;
            for( const JFacePtr &face: faces)
                sum = sum + JFaceGeometry::getAngleAt(face,vtx, ANGLE_IN_DEGREES);
            if( sum < convexAngle || sum > concaveAngle) {
                vtx->setAttribute("CornerAngle", sum);
                nodes.push_back(vtx);
            }
        }
    }
    return nodes;
}
///////////////////////////////////////////////////////////////////////////////

/*
JTinyStatistics JMeshGeometry :: getBoundaryEdgesLength() const
{
    logger->setInfo("Collecting minimum edge length");

     vector<double>  val;

    size_t numedges = mesh->getSize(1);
    for( size_t i  = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() && edge->isBoundary() ) {
            double l2  = JEdgeGeometry::getLength2(edge);
            val.push_back(l2);
        }
    }
    JTinyStatistics jstat;
    jstat.setData( val );
    return jstat;
}
//////////////////////////////////////////////////////////////////////////////

JTinyStatistics JMeshGeometry :: getEdgeLength() const
{
    logger->setInfo("Collecting minimum edge length");
    vector<double>  val;

    size_t numedges = mesh->getSize(1);
    for( size_t i  = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            double l2  = JEdgeGeometry::getLength2(edge);
            val.push_back(l2);
        }
    }
    JTinyStatistics jstat;
    jstat.setData( val );
    return jstat;
}
*/

