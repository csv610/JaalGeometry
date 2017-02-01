
#include "IsotropicMesh.hpp"

#ifdef USE_CGAL
namespace PMP = CGAL::Polygon_mesh_processing;
#endif

JMeshPtr JIsotropicMesh :: getMesh()
{
#ifdef USE_CGAL
    JMeshIO mio;
    mio.saveAs(inMesh, "tmp.off");
    std::ifstream input("tmp.off");
    CGALPolyhedron poly;
    if ( !input || !(input >> poly) || poly.empty() ) {
        std::cerr << "Not a valid off file." << std::endl;
        return nullptr;
    }
   
    std::vector<CGALPolyhedron::Facet_handle>  new_facets;
    std::vector<Vertex_handle> new_vertices;
    PMP::refine(poly, faces(poly),
                std::back_inserter(new_facets),
                std::back_inserter(new_vertices),
                PMP::parameters::density_control_factor(100.0));

    std::ofstream refined_off("iso.off");
    refined_off << poly;
    refined_off.close();

    Mesh mesh;
    if (!input || !(input >> mesh)) {
        std::cerr << "Not a valid off file." << std::endl;
        return nullptr;
    }

        unsigned int nb_iter = 5;
        std::vector<edge_descriptor> border;

        PMP::border_halfedges(faces(mesh), mesh,
                              boost::make_function_output_iterator(halfedge2edge(mesh, border)));

        PMP::split_long_edges(border, desiredEdgeLength, mesh);

        PMP::isotropic_remeshing( faces(mesh), desiredEdgeLength, mesh,
                                  PMP::parameters::number_of_iterations(nb_iter)
                                  .protect_constraints(true)//i.e. protect border, here
                                );

    std::ofstream outfile("iso.off");
    outfile << mesh;
    outfile.close();

    JMeshPtr isomesh = JMeshIO::readFile("iso.off");
    return isomesh;
#endif
    return nullptr;
}
