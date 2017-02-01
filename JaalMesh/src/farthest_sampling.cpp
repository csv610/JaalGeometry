/* read mesh from file and - if one vertex is specified, for all vertices of the mesh print their distances to this vertex
 - if two vertices are specified, print the shortest path between these vertices

	Danil Kirsanov, 01/2008
*/
#include <iostream>
#include <fstream>
#include "geodesic_algorithm_exact.h"
#include "geodesic_algorithm_dijkstra.h"
#include <vector>
#include <map>
#include <string>

using namespace std;

std::vector<double>    points;
std::vector<unsigned>  faces;
std::vector<unsigned>  sampled;
int numSamples;
std::vector<int>       faceGroup;

int read_off_mesh( const string &file)
{
    std::ifstream infile( file.c_str(), ios::in);
    if( infile.fail() ) {
        exit(0);
    }

    int numnodes, numfaces, numedges;
    std::string str;
    infile >> str;
    assert( str == "OFF");
    infile >> numnodes >> numfaces >> numedges;
    points.resize(3*numnodes);
    double x, y, z;
    for( int i = 0; i < numnodes; i++) {
        infile >> x >> y >> z;
        points[3*i]    = x;
        points[3*i+1]  = y;
        points[3*i+2]  = z;
    }

    int nn;
    faces.resize(3*numfaces);
    int n1, n2, n3;
    for( int i = 0; i < numfaces; i++)  {
        infile >> nn >> n1 >> n2 >> n3;
        faces[3*i]   = n1;
        faces[3*i+1] = n2;
        faces[3*i+2] = n3;
        assert(nn == 3);
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void farthest_samples()
{
    geodesic::Mesh mesh;
    mesh.initialize_mesh_data(points, faces);

    sampled.clear();
    geodesic::GeodesicAlgorithmExact algorithm(&mesh);
    std::vector<geodesic::SurfacePoint> all_sources;

    int new_index  = 1;
    double distance, maxdist = 0;
    for( int j = 0; j < numSamples; j++) {
        sampled.push_back(new_index);
        geodesic::SurfacePoint newsource(&mesh.vertices()[new_index]);		//create source
        all_sources.push_back(newsource);
        algorithm.propagate(all_sources);	//cover the whole mesh
        maxdist = 0.0;
        for(unsigned i=0; i<mesh.vertices().size(); ++i)
        {
            geodesic::SurfacePoint p(&mesh.vertices()[i]);
            unsigned best_source = algorithm.best_source(p,distance);
            if( distance > maxdist) {
                maxdist = distance;
                new_index = i;
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////////
void voronoi_regions()
{
    std::vector<double>    newPoints;
    std::vector<unsigned>  newFaces;
    newPoints = points;

    size_t numfaces  = faces.size()/3;
    size_t numPoints = points.size()/3;
    int conn[3];
    for( size_t i = 0; i < numfaces; i++) {
        conn[0] = faces[3*i+0];
        conn[1] = faces[3*i+1];
        conn[2] = faces[3*i+2];
        double x = 0;
        double y = 0;
        double z = 0;
        for( int j = 0; j < 3; j++) {
            x += points[3*conn[j] + 0];
            y += points[3*conn[j] + 1];
            z += points[3*conn[j] + 2];
        }
        newPoints.push_back(x/3.0);
        newPoints.push_back(y/3.0);
        newPoints.push_back(z/3.0);
        int n4 =  numPoints + i;
        for( int j = 0; j < 3; j++) {
            newFaces.push_back(n4);
            newFaces.push_back(conn[(j+1)%3]);
            newFaces.push_back(conn[(j+2)%3]);
        }
    }

    geodesic::Mesh newMesh;
    newMesh.initialize_mesh_data(newPoints, newFaces);
    geodesic::GeodesicAlgorithmDijkstra algorithm(&newMesh);

    std::vector<geodesic::SurfacePoint> all_sources;
    for( int j = 0; j < sampled.size(); j++) {
        int new_index  = sampled[j];
        geodesic::SurfacePoint newsource(&newMesh.vertices()[new_index]);		//create source
        all_sources.push_back(newsource);
    }
    algorithm.propagate(all_sources);	//cover the whole mesh

    faceGroup.resize(numfaces);
    double distance;
    for(int i=0; i< numfaces; ++i)
    {
        int id = numPoints + i;
        geodesic::SurfacePoint p(&newMesh.vertices()[id]);
        unsigned best_source = algorithm.best_source(p,distance);
        faceGroup[i] = best_source;
    }
}
//////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if(argc != 3)
    {
        std::cout << "usage: mesh_file_name  #samples" << std::endl;
        return 0;
    }

    read_off_mesh(argv[1]);
    numSamples = atoi( argv[2] );

    farthest_samples();
    voronoi_regions();

    for( size_t i = 0; i < faceGroup.size(); i++)
        cout << i << " Group " << faceGroup[i] << endl;

    return 0;
}
