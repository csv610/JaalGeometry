#include "MeshShapeDiameterSegmentation.hpp"

///////////////////////////////////////////////////////////////////////////////

void JMeshShapeDiameterSegmentation :: storeTriMesh()
{
    if( mesh == nullptr) return;

    int dim = mesh->getTopology()->getDimension();
    JMeshPtr  surfMesh;

    // If the mesh have 3 cells, extract the boundary mesh ...
    if( dim == 3)
        surfMesh = mesh->getTopology()->getSurfaceMesh();
    else
        surfMesh = mesh;

    JMeshPtr  triMesh;
    AllTriMeshGenerator alltri;

    int elemType = surfMesh->getTopology()->getElementsType(2);
    if( elemType == JFace::QUADRILATERAL) {
        JMeshPtr tmpmesh = surfMesh->deepCopy();
        triMesh =  alltri.getFromQuadMesh(tmpmesh,4);
    } else if( elemType == JFace::TRIANGLE)
        triMesh = surfMesh->deepCopy();

    triMesh->enumerate(0);
    JMeshIO meshio;
    meshio.saveAs(triMesh, "tmp.off");
}

///////////////////////////////////////////////////////////////////////////////

int JMeshShapeDiameterSegmentation :: getPartitions()
{
    if( mesh == nullptr) return 0;

#ifdef USE_CGAL
    storeTriMesh();

    std::ifstream input( "tmp.off");
    Polyhedron tmesh;
    input >> tmesh;
    Skeleton skeleton;

    CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);
    std::ofstream output("skel.cgal");
    SkelPolylines display(skeleton, output);
    CGAL::split_graph_into_polylines(skeleton, display);
    output.close();

// init the polyhedron simplex indices
    CGAL::set_halfedgeds_items_id(tmesh);

//for each input vertex compute its distance to the skeleton
    std::vector<double> distances(num_vertices(tmesh));
    BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton) )
    {
        const Point& skel_pt = skeleton[v].point;
        BOOST_FOREACH(vertex_descriptor mesh_v, skeleton[v].vertices)
        {
            const Point& mesh_pt = mesh_v->point();
            distances[mesh_v->id()] = std::sqrt(CGAL::squared_distance(skel_pt, mesh_pt));
        }
    }

    // create a property-map for sdf values
    std::vector<double> sdf_values( num_faces(tmesh) );
    Facet_with_id_pmap<double> sdf_property_map(sdf_values);
    // compute sdf values with skeleton
    BOOST_FOREACH(face_descriptor f, faces(tmesh))
    {
        double dist = 0;
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(f, tmesh), tmesh))
        dist+=distances[target(hd, tmesh)->id()];
        sdf_property_map[f] = dist / 3.;
    }

    // post-process the sdf values
    CGAL::sdf_values_postprocessing(tmesh, sdf_property_map);
    // create a property-map for segment-ids (it is an adaptor for this case)
    std::vector<std::size_t> segment_ids( num_faces(tmesh) );
    Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);
    // segment the mesh using default parameters
    int nCount = CGAL::segmentation_from_sdf_values(tmesh, sdf_property_map, segment_property_map, numClusters);

    size_t numFaces =  mesh->getSize(2);
    assert( numFaces == segment_ids.size() );

    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        int pid = segment_ids[i];
        f->setAttribute("Partition", pid );
    }
    return 0;
#endif
    return 1;
}

