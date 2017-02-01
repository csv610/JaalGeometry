#include "MeshSDF.hpp"

int  JMeshSDF :: execute()
{
    if( mesh == nullptr) return 1;

#ifdef USE_CGAL
    CGALPolyhedron poly;
    JMeshIO meshio;
    meshio.saveAs(mesh, "tmp.off");
    std::ifstream input("tmp.off");
    input >> poly;

    // create a property-map for segment-ids
    typedef std::map<CGALPolyhedron::Facet_const_handle, std::size_t> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);
    // calculate SDF values and segment the mesh using default parameters.

    double coneAngle = 2.0*M_PI/3.0;
    int    numRays    = 25;
    std::size_t number_of_segments = CGAL::segmentation_via_sdf_values(poly, segment_property_map, 
                                                                coneAngle, numRays, numClusters);

    size_t numfaces = mesh->getSize(2);
    vector<int> segmentid(numfaces);
    size_t index = 0;

    for(CGALPolyhedron::Facet_const_iterator facet_it = poly.facets_begin();
            facet_it != poly.facets_end(); ++facet_it) {
        segmentid[index++] = segment_property_map[facet_it];
    }

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->setAttribute("Partition", segmentid[i]);
    }
    return 0;
#endif
    return 1;
}

