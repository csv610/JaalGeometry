#pragma once

#include "Mesh.hpp"
#include "GeomPredicates.hpp"
#include "circumcenter.hpp"
#include "SwapEdges.hpp"

#define REAL double
#define ANSI_DECLARATORS
extern "C" {
#include <triangle.h>
}

namespace Jaal {

struct JDelaunayMesh2D {
    // Check if the edge is a Delunay. If the Coordinates are (x,y,z = 0), circumcircle
    // otherwise circumsphere test is performed.
    // Return value:
    //     -1 :   Can not perform the test, as necessary data structure not available.
    //      0 :   Not a Delaunay, some point inside the circle(sphere);
    //      1 :   A Delaunay: empty circum circle(sphere);
    // Preconditions:
    //      Mesh must be counterclockwise oriented ( See the triangle document );
    //      Edge-Faces relations must be present.

    // Check if the  face is a Delunay. If the Coordinates are (x,y,z = 0), circumcircle
    // otherwise circumsphere test is performed.
    // Return value:
    //     -1 :   Can not perform the test, as necessary data structure not available.
    //      0 :   Not a Delaunay, some point inside the circle(sphere);
    //      1 :   A Delaunay: empty circum circle(sphere);
    //      Mesh must be counterclockwise oriented ( See the triangle document );
    //      Edge-Faces relations must be present.

    static bool isDelaunay(const JEdgePtr &edge);
    static bool isDelaunay(const JFacePtr  &f);

    JDelaunayMesh2D()
    {
        exactinit();
        boundarySplit = 0;
    }

    void  setMesh( const JMeshPtr &m) {
        mesh  = m;
    }

    size_t countNonDelaunayEdges();
    size_t countNonDelaunayFaces();

    // Check if the surface triangulation is Delaunay. Works only for
    // triangle mesh. We use Jonathan Shewchuk's predicates and all his
    // ideas. For 2D check with Circum-Circle and for 2D Equitorial Sphere.
    //
    //
    // If the surface triangulation is not Delaunay, use edge-flips to
    // tranform the mesh into Delaunay.
    //
    void  setCreaseAngle(double a) {
        creaseAngle = a;
    }
    void  setMinAngle( double minA)
    {
        minAngle = minA;
    }

    void setMaxArea( double area)
    {
        maxArea = area;
    }

    void setBoundarySplit( bool v ) {
        boundarySplit = v;
    }

    int   isDelaunay();
    void  addPoints( const vector<Point2D> &p );
    void  addSegments( JEdgeSequence  &e);
    void  addSegments( vector<JEdgeSequence> &e);
    int   retriangulate();

    JMeshPtr getConvexHull();
    JMeshPtr getSimpleMesh();
    JMeshPtr getQualityMesh();
    JMeshPtr getDual();
    JMeshPtr getMedialAxis();
    int  getRemeshed();
    int  getIntrinsicMesh();

private:
    JMeshPtr mesh, newmesh;
    string   options;
    double   minAngle = 34.0;
    double   maxArea  = 0.0;
    double   creaseAngle = 30.0;

    void    build( vector<double> &uvCoords, vector<int> &segments,
                   vector<double> &holeCoords);
    bool    boundarySplit;
    int     flipEdge( const JEdgePtr &e);
};

////////////////////////////////////////////////////////////////////////////////////////

struct JDelaunayMesh3D {
    // Check if the edge is a Delunay. If the Coordinates are (x,y,z = 0), circumcircle
    // otherwise circumsphere test is performed.
    // Return value:
    //     -1 :   Can not perform the test, as necessary data structure not available.
    //      0 :   Not a Delaunay, some point inside the circle(sphere);
    //      1 :   A Delaunay: empty circum circle(sphere);
    // Preconditions:
    //      Mesh must be counterclockwise oriented ( See the triangle document );
    //      Edge-Faces relations must be present.
    static int  isDelaunay(const JEdgePtr &edge);

    // Check if the  face is a Delunay. If the Coordinates are (x,y,z = 0), circumcircle
    // otherwise circumsphere test is performed.
    // Return value:
    //     -1 :   Can not perform the test, as necessary data structure not available.
    //      0 :   Not a Delaunay, some point inside the circle(sphere);
    //      1 :   A Delaunay: empty circum circle(sphere);
    //      Mesh must be counterclockwise oriented ( See the triangle document );
    //      Edge-Faces relations must be present.
    static int  isDelaunay3D(const JFacePtr  &f);
    static int  isDelaunay(const JCellPtr  &f);

    static size_t countNonDelaunayCells( const JMeshPtr &m);

    JDelaunayMesh3D()
    {
        simplex = 3;  // By default generate triangles in 2D and 3D (surface);
        exactinit();
    }

    // Check if the surface triangulation is Delaunay. Works only for
    // triangle mesh. We use Jonathan Shewchuk's predicates and all his
    // ideas. For 2D check with Circum-Circle and for 2D Equitorial Sphere.
    //
    int  isDelaunay(const JMeshPtr &m);
    //
    // If the surface triangulation is not Delaunay, use edge-flips to
    // tranform the mesh into Delaunay.
    //
    void  setDihedralAngle ( double a )
    {
        dihedralAngle = a;
    }

    int   retriangulate(const JMeshPtr &m);

    JMeshPtr addPoints( const vector<Point3D> &p );

    JMeshPtr addFacets( const JEdgeSequence    &e);
    JFaceSequence getConvexHull( const vector<Point3D> &p);

    JMeshPtr makeDelaunay( const vector<Point3D> &p );

    int   refine(const JMeshPtr &m);

    int   insert(const Point3D &p);
    int   remove(JNodePtr &v);

    JMeshPtr toTriVoronoi( const JMeshPtr &m );  // Convert Delaunay Mesh to Voronoi
    JMeshPtr fromTriVoronoi(const JMeshPtr &m ); // Convert Voronoi to Delaunay

private:
    JMeshPtr mesh;
    void   get_simple_triangulation( const vector<Point3D> &p);
    double dihedralAngle;
    int    simplex;    // Triangle, Tets etc...
    void   build( vector<double> &uvCoords, vector<int> &segments,
                  vector<double> &holeCoords);
};

}
