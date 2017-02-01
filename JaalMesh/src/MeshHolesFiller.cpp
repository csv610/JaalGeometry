#include "MeshHolesFiller.hpp"

void JMeshHolesFiller::setMesh( const JMeshPtr &m)
{
    mesh = m;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshHolesFiller :: fill( const JEdgeSequence &edges)
{

#ifdef USE_CGAL
    if( nullptr == mesh ) return;

    JNodeSequence nodes;
    JEdgeTopology::getChainNodes( edges, nodes);

    std::vector<CGALPoint> polyline;
    for( const JNodePtr &vtx : nodes) {
        const Point3D &xyz = vtx->getXYZCoords();
        polyline.push_back(CGALPoint( xyz[0], xyz[1], xyz[2]));
    }

    typedef CGAL::Triple<int, int, int> TriConnect;
    std::vector<TriConnect> patch;
    patch.reserve(polyline.size() -2); // there will be exactly n-2 triangles in the patch
    CGAL::Polygon_mesh_processing::triangulate_hole_polyline(
        polyline, std::back_inserter(patch));

    JNodeSequence tconn(3);
    for(std::size_t i = 0; i < patch.size(); ++i)
    {
        tconn[0] = nodes[ patch[i].first ];
        tconn[1] = nodes[ patch[i].second ];
        tconn[2] = nodes[ patch[i].third ];
        JFacePtr  face = Triangle::newObject( tconn );
        mesh->addObject(face);
    }
#endif
}
///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshHolesFiller::fillAll()
{
    if( nullptr == mesh ) return nullptr;

#ifdef USE_CGAL

    JMeshIO mio;
    mio.saveAs( mesh, "tmp.off");
    std::ifstream input("tmp.off");

    CGALPolyhedron poly;
    input >> poly;

    size_t nb_holes = 0;
    BOOST_FOREACH(Halfedge_handle h, halfedges(poly))
    {
        if(h->is_border())
        {
            std::vector<Facet_handle>  patch_facets;
            std::vector<Vertex_handle> patch_vertices;
            bool success = CGAL::cpp11::get<0>(
                               CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                                   poly,
                                   h,
                                   std::back_inserter(patch_facets),
                                   std::back_inserter(patch_vertices),
                                   CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, poly)).
                                   geom_traits(Kernel())) );
            std::cout << " Number of facets in constructed patch: " << patch_facets.size() << std::endl;
            std::cout << " Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
            std::cout << " Fairing : " << (success ? "succeeded" : "failed") << std::endl;
            ++nb_holes;
        }
    }
    std::cout << std::endl;
    std::cout << nb_holes << " holes have been filled" << std::endl;

    std::ofstream out("filled.off");
    out.precision(17);
    out << poly << std::endl;

    JMeshOFFImporter mimp;
    JMeshPtr newmesh = mimp.readFile("filled.off");
    return newmesh;
#endif
    return nullptr;
}

