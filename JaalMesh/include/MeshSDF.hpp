#pragma once

#include "Mesh.hpp"

#ifdef USE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
#endif

class JMeshSDF
{
#ifdef USE_CGAL
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel> CGALPolyhedron;
#endif
public:
    void setMesh(const JMeshPtr &m) {
        mesh = m;
    }

    int  setNumOfClusters(int n) { numClusters = n; }
    int  getNumOfClusters() const;
    int  execute();
private:
    JMeshPtr mesh;
    int numClusters = 10;
};
