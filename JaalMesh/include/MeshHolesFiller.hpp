#pragma once

#ifdef USE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>

#include "Mesh.hpp"
#include "MeshExporter.hpp"
#include "MeshImporter.hpp"

class JMeshHolesFiller
{
#ifdef USE_CGAL
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel>     CGALPolyhedron;
    typedef CGALPolyhedron::Halfedge_handle   Halfedge_handle;
    typedef CGALPolyhedron::Facet_handle      Facet_handle;
    typedef CGALPolyhedron::Vertex_handle     Vertex_handle;
    typedef Kernel::Point_3 CGALPoint;
    CGALPolyhedron poly;
#endif

public:
    void setMesh( const JMeshPtr &m);

    int  getNumOfHoles() const;
    void fill(const JEdgeSequence &e);
    JMeshPtr fillAll();
private:
    JMeshPtr mesh;
};
