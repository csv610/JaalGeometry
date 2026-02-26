#include <iostream>
#include "AllQuadMeshGenerator.hpp"

using namespace Jaal;

/**
 * Example: Structured Quadrilateral Mesh Generation and Simple Property Queries
 */
int main() {
    // 1. Define dimensions for a 10x10 grid of quad cells (11x11 vertices)
    int dim[2] = {11, 11};

    // 2. Generate the structured quad mesh
    // dim[0] * dim[1] vertices (11*11 = 121)
    // results in (dim[0]-1) * (dim[1]-1) quad faces (10*10 = 100)
    JMeshPtr quadMesh = AllQuadMeshGenerator::getStructuredMesh(dim);

    if (quadMesh) {
        // 3. Print basic mesh information
        std::cout << "Successfully generated a structured Quad Mesh." << std::endl;
        std::cout << "Number of Nodes: " << quadMesh->getSize(0) << std::endl; // Should be 121
        std::cout << "Number of Quadrilaterals: " << quadMesh->getSize(2) << std::endl; // Should be 100

        // 4. Access individual faces (quads)
        JFacePtr firstFace = quadMesh->getFaceAt(0);
        if (firstFace) {
            std::cout << "First Face has " << firstFace->getSize(0) << " vertices." << std::endl;
        }

        // 5. Example: Schneider's Pyramid (a common quad-based test shape)
        JMeshPtr pyramid = AllQuadMeshGenerator::SchneiderPyramid();
        if (pyramid) {
            std::cout << "Generated Schneider Pyramid with " 
                      << pyramid->getSize(2) << " quad faces." << std::endl;
        }
    } else {
        std::cerr << "Failed to generate Quad Mesh." << std::endl;
    }

    return 0;
}
