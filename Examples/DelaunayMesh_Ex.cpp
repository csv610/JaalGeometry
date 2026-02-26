#include <iostream>
#include <vector>
#include "DelaunayMesh.hpp"

using namespace Jaal;

/**
 * Example: Generating a 2D Delaunay Triangle Mesh from a Point Set
 */
int main() {
    // 1. Create a set of 2D points (random-ish coordinates)
    std::vector<Point2D> points;
    points.push_back({0.0, 0.0});
    points.push_back({1.0, 0.0});
    points.push_back({1.0, 1.0});
    points.push_back({0.0, 1.0});
    points.push_back({0.5, 0.5}); // Center point
    points.push_back({0.2, 0.8});
    points.push_back({0.8, 0.2});

    std::cout << "Creating Delaunay Mesh from " << points.size() << " points..." << std::endl;

    // 2. Initialize the Delaunay Mesh generator
    JDelaunayMesh2D delGen;
    
    // 3. Add the points to the generator
    delGen.addPoints(points);

    // 4. Generate the simple Delaunay triangulation (Convex Hull by default)
    JMeshPtr mesh = delGen.getSimpleMesh();

    if (mesh) {
        std::cout << "Successfully generated Delaunay Mesh." << std::endl;
        std::cout << "Number of Nodes: " << mesh->getSize(0) << std::endl;
        std::cout << "Number of Triangles: " << mesh->getSize(2) << std::endl;

        // 5. Example: Perform quality-constrained Delaunay refinement
        // You can set min angle and max area before calling getQualityMesh()
        delGen.setMinAngle(30.0);
        JMeshPtr qualityMesh = delGen.getQualityMesh();
        
        if (qualityMesh) {
            std::cout << "Quality Refinement complete." << std::endl;
            std::cout << "Refined Mesh Triangles: " << qualityMesh->getSize(2) << std::endl;
        }
    } else {
        std::cerr << "Failed to generate Delaunay Mesh." << std::endl;
    }

    return 0;
}
