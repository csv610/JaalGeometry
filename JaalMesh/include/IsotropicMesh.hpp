#pragma once

#include "Mesh.hpp"

#include <boost/property_map/property_map.hpp>

#ifdef USE_CGAL
#include <CGAL/algorithm.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Cartesian.h>
#endif

#include <boost/function_output_iterator.hpp>

class JIsotropicMesh
{
#ifdef USE_CGAL
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Surface_mesh<K::Point_3> Mesh;
    typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
    typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;

    typedef CGAL::Polyhedron_3<K>  CGALPolyhedron;
    typedef CGALPolyhedron::Vertex_handle   Vertex_handle;

    struct halfedge2edge
    {
        halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
            : m_mesh(m), m_edges(edges)
        {}
        void operator()(const halfedge_descriptor& h) const
        {
            m_edges.push_back(edge(h, m_mesh));
        }
        const Mesh& m_mesh;
        std::vector<edge_descriptor>& m_edges;
    };
#endif
public:
    void setMesh(const JMeshPtr &m) { inMesh = m; }
    void setBoundarySplit( bool b ) { allowBoundarySplit = b; }
    void setEdgeLength(double elen) { desiredEdgeLength = elen; }

    JMeshPtr  getMesh();
private:
    JMeshPtr  inMesh;
    JMeshPtr  newMesh;
    double    desiredEdgeLength = 0.1;
    bool      allowBoundarySplit = 0;

    void  refineEdge( const JEdgePtr &e);
    void  refineFace( const JFacePtr &e);
    void  getTriangles(double ien[3], int nodes[3], vector<int> &trinodes);
};

