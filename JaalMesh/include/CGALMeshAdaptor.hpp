#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include "Mesh.hpp"

class JCGALMeshAdaptor
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel> CGALPolyhedron;
public:
    CGALPolyhedron  getMesh( const JMeshPtr &m);
    JMeshPtr        getMesh( const CGALPolyhedron &m);
};
