#include "PolygonPartitioner.hpp"

//////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_CGAL
void JCGALPolygonPartitioner :: makePolygon(Polygon_2 &polygon)
{
    if( mesh == nullptr) return;

    vector<JEdgeSequence>  boundary;
    mesh->getTopology()->getBoundary(boundary);

    JNodeSequence boundnodes;
    for( size_t i = 0; i < boundary.size(); i++) {
        JEdgeTopology::getChainNodes( boundary[i], boundnodes);
        int ori = JEdgeGeometry::getOrientation( boundary[i] );
        if( ori < 0.0)  boost::reverse( boundnodes );

        for( const JNodePtr &vtx : boundnodes) {
            const Point3D &xyz = vtx->getXYZCoords();
            polygon.push_back(Point_2(xyz[0], xyz[1] ));
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////

JMeshPtr JCGALPolygonPartitioner :: getCGALPolyMesh()
{
    if( mesh == nullptr) return nullptr;

    Polygon_2     polygon;
    PolygonList   partitions;
    Traits        partition_traits;

    makePolygon(polygon);

    /*
        CGAL::approx_convex_partition_2(polygon.vertices_begin(), polygon.vertices_end(),
                                        std::back_inserter(partitions), partition_traits);
    */

    CGAL::optimal_convex_partition_2(polygon.vertices_begin(), polygon.vertices_end(),
                                     std::back_inserter(partitions), partition_traits);

    JMeshPtr  polymesh = JMesh::newObject();
    JNodeSequence nodes = mesh->getNodes();
    polymesh->addObjects( nodes);

    JNearestNeighbours jnn;
    jnn.setMesh( mesh );

    Point3D xyz;
    JNodeSequence facenodes;
    for( auto &poly :  partitions)
    {
        facenodes.clear();
        for (VertexIterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi) {
            xyz[0] = vi->x();
            xyz[1] = vi->y();
            xyz[2] = 0.0;
            JNodePtr vtx = jnn.getNearest( xyz );
            facenodes.push_back(vtx);
        }
        cout << endl;
        JFacePtr  newface = JFace::newObject( facenodes);
        polymesh->addObject(newface);
    }

    return polymesh;

}
//////////////////////////////////////////////////////////////////////////////////////
#endif
