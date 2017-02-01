#include "MeshSlicer.hpp"

void JMeshSlicer :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

#ifdef USE_CGAL
    JMeshIO mio;
    mio.saveAs(mesh, "tmp.off");

    std::ifstream input("tmp.off");
    input >> cgalMesh;
#endif
}

///////////////////////////////////////////////////////////////////////

#ifdef USE_CGAL
JEdgeSequence JMeshSlicer :: getEdges( const PolyPoints  &polypoints)
{
    JEdgeSequence contour;

    int numPoints = polypoints.size();
    JNodeSequence nodes(numPoints);
    contour.resize(numPoints);
    Point3D xyz;
    for( int i = 0; i < numPoints; i++) {
        xyz[0] = polypoints[i][0];
        xyz[1] = polypoints[i][1];
        xyz[2] = polypoints[i][2];
        nodes[i] = JNode::newObject();
        nodes[i]->setXYZCoords(xyz);
        nodes[i]->setID(i);
    }
    for( int i = 0; i < numPoints; i++) {
        const JNodePtr v0 = nodes[i];
        const JNodePtr v1 = nodes[(i+1)%numPoints];
        contour[i] = JEdge::newObject(v0,v1);
        contour[i]->setID(i);
    }
    return contour;
}
#endif

///////////////////////////////////////////////////////////////////////
vector<JEdgeSequence> JMeshSlicer :: getContours( const Vec3D &normal, const Point3D &passThru)
{
   vector<JEdgeSequence> contours;
#ifdef USE_CGAL
   double nx = normal[0];
   double ny = normal[1];
   double nz = normal[2];
   double nl = sqrt(nx*nx + ny*ny + nz*nz);
   
   double a = nx/nl;
   double b = ny/nl;
   double c = nz/nl;
   double d = -a*passThru[0] - b*passThru[1] - c*passThru[2];

    CGAL::Polygon_mesh_slicer<CGALMesh, K> slicer(cgalMesh);

    Polylines polylines;
    slicer(K::Plane_3(a, b, c, d), std::back_inserter(polylines));

    int nsize = polylines.size();

    contours.resize(nsize);
    int index = 0;
    for( const PolyPoints &p : polylines)
        contours[index++] = getEdges( p );
    return contours;
#endif
}
///////////////////////////////////////////////////////////////////////
