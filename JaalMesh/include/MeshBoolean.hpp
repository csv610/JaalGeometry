#pragma once

#include "Mesh.hpp"
#include "MeshMatrix.hpp"
#include <Eigen/Core>

#ifdef USE_CGAL
#include <boost/property_map/property_map.hpp>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#ifdef USE_IGL
#include <igl/copyleft/cgal/mesh_boolean.h>
#endif
#endif

#ifdef USE_IGL
#include <igl/copyleft/cork/mesh_boolean.h>
#endif

class JMeshBoolean
{
#ifdef USE_CGAL
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::FT FT;
    typedef K::Point_3 Point;
    typedef K::Vector_3 Vector;
    typedef K::Segment_3 Segment;
    typedef CGAL::Surface_mesh<Point> CGALMesh;
    typedef boost::graph_traits<CGALMesh>::face_descriptor face_descriptor;
    typedef boost::graph_traits<CGALMesh>::halfedge_descriptor halfedge_descriptor;
    typedef CGAL::AABB_face_graph_triangle_primitive<CGALMesh> Primitive;
    typedef CGAL::AABB_traits<K, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef boost::optional<Tree::Intersection_and_primitive_id<Segment>::Type> Segmen_intersection;
#endif

public:
    static const int MESH_UNION = 0;
    static const int MESH_COMPLEMENT   = 1;
    static const int MESH_INTERSECTION = 2;
    static const int MESH_DIFFERENCE   = 3;
    static const int MESH_SYMMETRIC_DIFFERENCE   = 4;
    static const int MESH_SPLIT_AT_INTERSECATION = 5;

    static int  getOp( const string &op);

    void setMesh(const JMeshPtr &a, const JMeshPtr &b);
    void setMesh1(const JMeshPtr &a);
    void setMesh2(const JMeshPtr &b);

    bool anyIntersection();

    JMeshPtr getMesh(int op);
    int getSkinMesh( const JMeshPtr &m);

private:
    JMeshPtr aMesh, bMesh;

    Eigen::MatrixXd VA, VB, VC;
    Eigen::MatrixXi FA, FB, FC;

    void splitAlongIntersection();
};
