#include <iostream>
#include "AlphaMSTQuadMesh.hpp"
#include "AllQuadMeshGenerator.hpp"

using namespace Jaal;

/**
 * Example: Local Quad Mesh Remeshing using JAlphaMSTQuadMesh
 */
int main() {
    // 1. Generate a large structured quad mesh (10x10 cells)
    int dim[2] = {11, 11};
    JMeshPtr mesh = AllQuadMeshGenerator::getStructuredMesh(dim);
    
    if (!mesh) {
        std::cerr << "Failed to generate initial mesh." << std::endl;
        return 1;
    }
    std::cout << "Original Mesh: " << mesh->getSize(0) << " nodes, " << mesh->getSize(2) << " faces." << std::endl;

    // 2. Initialize the AlphaMST remesher
    JAlphaMSTQuadMesh remesher;
    remesher.setMesh(mesh);

    // 3. Define a circular region to remesh in the middle of the 1.0x1.0 square
    JCircle circle;
    circle.setCenter({0.5, 0.5, 0.0});
    circle.setRadius(0.25);
    remesher.setCircle(circle);

    // 4. Build the patch (identifies faces/nodes within the circle)
    remesher.buildPatch();
    
    if (remesher.isEmpty()) {
        std::cout << "No patch found in the specified region." << std::endl;
    } else {
        std::cout << "Patch identified with " << remesher.getFaces().size() << " faces." << std::endl;
        std::cout << "Singularities in patch before remeshing: " << remesher.getNumSingularities() << std::endl;

        // 5. Perform the local remeshing
        // This replaces the old patch with a new quad template and smooths it.
        remesher.remeshPatch();

        std::cout << "Remeshing complete." << std::endl;
        std::cout << "Final Mesh: " << mesh->getSize(0) << " nodes, " << mesh->getSize(2) << " faces." << std::endl;
        std::cout << "New singularities in remeshed patch: " << remesher.getNumSingularities() << std::endl;
    }

    return 0;
}
