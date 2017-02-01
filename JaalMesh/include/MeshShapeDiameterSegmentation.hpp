#pragma once

#include <algorithm>
#include <fstream>
#include <Eigen/Core>

#ifdef USE_CGAL
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#endif

#include "Mesh.hpp"
#include "AllTriMeshGenerator.hpp"

class JMeshShapeDiameterSegmentation
{
#ifdef USE_CGAL
    typedef CGAL::Simple_cartesian<double>                               Kernel;
    typedef Kernel::Point_3                                              Point;
    typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;
    typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
    typedef boost::graph_traits<Polyhedron>::halfedge_descriptor         halfedge_descriptor;
    typedef boost::graph_traits<Polyhedron>::face_descriptor             face_descriptor;
    typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron>        Skeletonization;
    typedef Skeletonization::Skeleton                                    Skeleton;
    typedef Skeleton::vertex_descriptor                                  Skeleton_vertex;

    struct SkelPolylines {
        const Skeleton& skeleton;
        std::ofstream& out;
        int polyline_size;
        std::stringstream sstr;
        SkelPolylines(const Skeleton& skeleton, std::ofstream& out)
            : skeleton(skeleton), out(out)
        {}
        void start_new_polyline() {
            polyline_size=0;
            sstr.str("");
            sstr.clear();
        }
        void add_node(Skeleton_vertex v) {
            ++polyline_size;
            sstr << " " << skeleton[v].point << endl;
        }
        void end_polyline()
        {
            out << polyline_size << endl;
            out << sstr.str()    << endl;
        }
    };
    template<class ValueType>
    struct Facet_with_id_pmap
        : public boost::put_get_helper<ValueType&,
          Facet_with_id_pmap<ValueType> >
    {
        typedef face_descriptor key_type;
        typedef ValueType value_type;
        typedef value_type& reference;
        typedef boost::lvalue_property_map_tag category;
        Facet_with_id_pmap(
            std::vector<ValueType>& internal_vector
        ) : internal_vector(internal_vector) { }
        reference operator[](key_type key) const
        {
            return internal_vector[key->id()];
        }
    private:
        std::vector<ValueType>& internal_vector;
    };
#endif

public:
    void setMesh(const JMeshPtr &m) {
        mesh = m;
    }
    void setNumCluster( int n ) { numClusters = n; }

    int getPartitions();

private:
    JMeshPtr mesh;
    int  numClusters = 10;
    void storeTriMesh();
};

