#include "MeshGeodesicVoronoiRegions.hpp"


void JMeshGeodesicVoronoiRegions :: setMesh( const JMeshPtr &m)
{


}

void JMeshGeodesicVoronoiRegions ::   extracRegions( int numRegions )
{
    vector<double> dist(numnodes);

    // First collect the centers ...
    JNodeSequence  sources;
    size_t numnodes = mesh->getSize(0);
    source.resize(1);
    source[0] = mesh->getNodeAt(0);

    for( int i = 0; i < numRegions; i++) {
        meshGeodesic->getDistances(sources, dist);
        double maxdist = JMaths::max_element_at(dist);
        sources.push_back( mesh->getNodeAt(id));
    }
    meshGeodesic->getNearestSource(sources);
}

