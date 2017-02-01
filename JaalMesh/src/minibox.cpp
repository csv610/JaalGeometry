#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <memory.h>
#include  <math.h>

#include "MeshGeometry.hpp"
#include "MinDiam.hpp"

Sphere JMeshGeometry :: getMinimumSphere()
{
    Sphere sph;
#ifdef HAVE_MINIBALL
    mesh->getLogger()->setInfo("Calculating the mesh minimum sphere " );

    double seed = 0;
    std::srand (seed);

    std::list<std::vector<double> > lp;
    size_t numnodes = mesh->getSize(0);
    for (size_t i=0; i< numnodes; ++i) {
        std::vector<double> p(3);
        JNodePtr vtx = mesh->getNodeAt(i);
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
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> >
    MB;
    MB mb (3, lp.begin(), lp.end());

    const double*xyz = mb.center();

    Point3D center;
    center[0] = xyz[0];
    center[1] = xyz[1];
    center[2] = xyz[2];

    sph.setCenter( center );
    sph.setRadius( sqrt(mb.squared_radius()) );

#endif
    return sph;
}

////////////////////////////////////////////////////////////////////////////////

JHexahedronPtr JMeshGeometry :: getMinimumBox()
{
    if( mesh == nullptr) return nullptr;

    JHexahedronPtr hex;
    mesh->getLogger()->setInfo("Calculating the mesh minimum box  " );

    vector<double> points;
    vector<size_t> l2g;
    getCoordsArray(points, l2g);

    int numpoints = points.size()/3;
    if( numpoints == 0) return nullptr;

    gdiam_point  *pnt_arr = gdiam_convert( &points[0], numpoints);

    gdiam_bbox bb = gdiam_approx_mvbb( pnt_arr, numpoints, 1.0E-03);

    JNodeSequence nodes(8);
    for( int i = 0; i < 8; i++)
        nodes[i] = JNode::newObject();

    double x, y, z;

    bb.get_vertex(  0, 0, 0, &x, &y, &z);
    nodes[0]->setXYZCoords( x, y, z);

    bb.get_vertex(  1, 0, 0, &x, &y, &z);
    nodes[1]->setXYZCoords(x, y, z);

    bb.get_vertex(  1, 1, 0, &x, &y, &z);
    nodes[2]->setXYZCoords( x, y, z);

    bb.get_vertex(  0, 1, 0, &x, &y, &z);
    nodes[3]->setXYZCoords( x, y, z);

    bb.get_vertex(  0, 0, 1, &x, &y, &z);
    nodes[4]->setXYZCoords( x, y, z);

    bb.get_vertex(  1, 0, 1, &x, &y, &z);
    nodes[5]->setXYZCoords( x, y, z);

    bb.get_vertex(  1, 1, 1, &x, &y, &z);
    nodes[6]->setXYZCoords( x, y, z);

    bb.get_vertex(  0, 1, 1, &x, &y, &z);
    nodes[7]->setXYZCoords( x, y, z);

    hex = Hexahedron::newObject();
    hex->setNodes( nodes );
    free(pnt_arr);
    return hex;
}

#ifdef CSV

JHexahedronPtr JMeshGeometry :: getMinimumBox( const vector<Point3D> &inpoints)
{
    int numpoints = inpoints.size();

    gdiam_real  *points;
    points = (gdiam_point)malloc( sizeof( gdiam_point_t ) * numpoints);
    assert( points != NULL );

    for  ( size_t i = 0; i < numpoints; i++ ) {
        points[ ind * 3 + 0 ] = inpoints[i][0];
        points[ ind * 3 + 1 ] = inpoints[i][1];
        points[ ind * 3 + 2 ] = inpoints[i][2];
    }

    GPointPair  pair;
    pair = gdiam_approx_diam_pair( (gdiam_real *)points, numpoints, 0.0 );

    gdiam_point  *pnt_arr;
    pnt_arr = gdiam_convert( (gdiam_real *)points, numpoints );

    gdiam_bbox   bb;
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

    JHexahedronPtr hex = Hexahedron::newObject();
    hex->setNodes( nodes );

    free(pnt_arr);
    return hex;
}

//////////////////////////////////////////////////////////////////////////////

JHexahedronPtr JMeshGeometry :: getMinimumBox( const JNodeSequence &nodes)
{
    if( nodes.empty() ) retunr nullptr;

    vector<Point3D> points;
    size_t numpoints = nodes.size();
    points.reserve( numpoints);

    for( const JNodePtr &vtx : nodes) {
        if( vtx->isActive() ) points.push_back( vtx->getXYZCoords() );

        return getMinimumBox( points);
    }

///////////////////////////////////////////////////////////////////////////////

    JHexahedronPtr JMeshGeometry :: getMinimumBox()
    {
        if( mesh == nullptr) return nullptr;

        JNodeSequence nodes;
        mesh->getNode(nodes);
        return getMinimumBox( nodes);
    }

///////////////////////////////////////////////////////////////////////////////
#endif
