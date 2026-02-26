#pragma once

#ifdef USE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/partition_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#endif

#include <cassert>
#include <list>

#include "Mesh.hpp"
#include "NearestNeighbours.hpp"
#include "polypartition.hpp"

class JCGALPolygonPartitioner
{
#ifdef USE_CGAL
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Partition_traits_2<K>  Traits;
    typedef CGAL::Is_convex_2<Traits>    Is_convex_2;
    typedef Traits::Polygon_2            Polygon_2;
    typedef Traits::Point_2              Point_2;
    typedef Polygon_2::Vertex_iterator   VertexIterator;
    typedef std::list<Polygon_2>         PolygonList;
#endif

public:
    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    JMeshPtr   getPolyMesh();

private:
    JMeshPtr   mesh;
    vector<JPolygonPtr>  polygons;

#ifdef USE_CGAL
    Polygon_2            polygon;
    PolygonList          partition;
    size_t     getID( const Point_2 &p);
    void  makePolygon(Polygon_2 &p);
    JMeshPtr   getCGALPolyMesh();
#endif
};

