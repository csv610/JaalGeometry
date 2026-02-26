#include <iostream>
#include <vector>
#include "DelaunayMesh.hpp"
#include "MeshDualGraph.hpp"

using namespace Jaal;

/**
 * Example: Generating a 2D Voronoi Diagram
 * 
 * In Jaal, a Voronoi diagram is obtained by taking the dual graph
 * of a Delaunay triangulation.
 */
int main() {
    // 1. Create a set of 2D seed points for the Voronoi diagram
    std::vector<Point2D> points;
    points.push_back({0.1, 0.1});
    points.push_back({0.9, 0.1});
    points.push_back({0.9, 0.9});
    points.push_back({0.1, 0.9});
    points.push_back({0.5, 0.5});

    std::cout << "Step 1: Generating Delaunay Triangulation for " << points.size() << " seeds..." << std::endl;

    // 2. Generate the Delaunay Triangulation
    JDelaunayMesh2D delGen;
    delGen.addPoints(points);
    JMeshPtr delaunayMesh = delGen.getSimpleMesh();

    if (delaunayMesh) {
        std::cout << "Delaunay Mesh generated with " << delaunayMesh->getSize(2) << " triangles." << std::endl;

        // 3. Construct the Dual Graph (Voronoi Diagram)
        // By default, Jaal places dual nodes at centroids or circumcenters.
        JMeshDualGraph dualBuilder;
        dualBuilder.setMesh(delaunayMesh);
        dualBuilder.setBoundaryNodes(true); // Include Voronoi vertices on the boundary
        
        JMeshPtr voronoiDiagram = dualBuilder.getGraph();

        if (voronoiDiagram) {
            std::cout << "Step 2: Successfully generated Voronoi Diagram." << std::endl;
            std::cout << "Number of Voronoi Vertices: " << voronoiDiagram->getSize(0) << std::endl;
            std::cout << "Number of Voronoi Edges: " << voronoiDiagram->getSize(1) << std::endl;
            
            // Note: In 2D, the dual edges represent the boundaries between Voronoi cells.
        } else {
            std::cerr << "Failed to generate Dual Graph." << std::endl;
        }
    } else {
        std::cerr << "Failed to generate Delaunay Triangulation." << std::endl;
    }

    return 0;
}
