#pragma  once

/*

#include "Mesh.hpp"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

class MeshIntersection
{
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::FT FT;
    typedef K::Point_3 Point;
    typedef K::Vector_3 Vector;
    typedef K::Ray_3 Ray;
    typedef CGAL::Surface_mesh<Point> CGALMesh;
    typedef boost::graph_traits<CGALMesh>::face_descriptor face_descriptor;
    typedef boost::graph_traits<CGALMesh>::halfedge_descriptor halfedge_descriptor;
    typedef CGAL::AABB_face_graph_triangle_primitive<CGALMesh> Primitive;
    typedef CGAL::AABB_traits<K, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
public:
    void setMesh( const JMeshPtr &m);

    JFaceSequence getAllSegmentIntersections( const JEdgePtr &e);
    JFaceSequence getAllAntiPodalIntersections( const JFacePtr &e);
    JFacePtr      getFirstAntiPodalIntersction( const JFacePtr &e);
    JFacePtr      getFirstRayIntersection( const Vec3D  &r);
private:
    JMeshPtr   mesh;
    boost::scoped_ptr<Tree> tree;
};
*/



