#include <iostream>
#include <vector>
#include "AllTetMeshGenerator.hpp"

using namespace Jaal;

/**
 * Example: Generating a 3D Tetrahedral Mesh from a Point Set
 */
int main() {
    // 1. Create a set of 3D points
    JNodeSequence nodes;
    
    // Define 8 corners of a unit cube
    double coords[8][3] = {
        {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
        {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
    };

    for(int i = 0; i < 8; ++i) {
        JNodePtr v = JNode::newObject();
        v->setXYZCoords(coords[i][0], coords[i][1], coords[i][2]);
        nodes.push_back(v);
    }

    // Add a center point to make it more interesting
    JNodePtr vCenter = JNode::newObject();
    vCenter->setXYZCoords(0.5, 0.5, 0.5);
    nodes.push_back(vCenter);

    std::cout << "Creating Tet Mesh from " << nodes.size() << " points..." << std::endl;

    // 2. Initialize the Tetrahedral Mesh Generator
    AllTetMeshGenerator tetGen;

    // 3. Generate the tetrahedral mesh (Convex Hull of points)
    // internally calls tetgen to perform Delaunay tetrahedralization
    JMeshPtr tetMesh = tetGen.getConvexHull(nodes);

    if (tetMesh) {
        std::cout << "Successfully generated Tetrahedral Mesh." << std::endl;
        std::cout << "Number of Nodes: " << tetMesh->getSize(0) << std::endl;
        std::cout << "Number of Tetrahedrons: " << tetMesh->getSize(3) << std::endl;
    } else {
        std::cerr << "Failed to generate Tetrahedral Mesh." << std::endl;
    }

    return 0;
}
