#pragma once

#include <fstream>
#include "Mesh.hpp"

#ifdef USE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polygon_mesh_slicer.h>
#endif

class JMeshSlicer
{
#ifdef USE_CGAL
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Surface_mesh<K::Point_3> CGALMesh;
    typedef std::vector<K::Point_3> PolyPoints;
    typedef std::list< PolyPoints> Polylines;
    typedef CGAL::AABB_halfedge_graph_segment_primitive<CGALMesh> HGSP;
    typedef CGAL::AABB_traits<K, HGSP>    AABB_traits;
    typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;
#endif

public:
    void setMesh( const JMeshPtr &m);
    vector<JEdgeSequence> getContours(const Vec3D  &normal, const Point3D &passThru);

private:
    JMeshPtr mesh;

#ifdef USE_CGAL
    CGALMesh cgalMesh;
    JEdgeSequence getEdges( const PolyPoints &p);
#endif

};

