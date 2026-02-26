#include <iostream>
#include "AllTriMeshGenerator.hpp"

/**
 * Example: Structured Triangle Mesh Generation and Simple Property Queries
 */
int main() {
    // 1. Define dimensions for a 10x10 structured triangle mesh (grid of vertices)
    int nx = 10;
    int ny = 10;

    // 2. Generate the structured triangle mesh
    // A grid of 10x10 vertices results in 9x9 quad cells, each split into 2 triangles.
    // Total triangles expected: 9 * 9 * 2 = 162
    JMeshPtr triMesh = AllTriMeshGenerator::getStructuredMesh(nx, ny);

    if (triMesh) {
        // 3. Print basic mesh information
        std::cout << "Successfully generated a structured Triangle Mesh." << std::endl;
        std::cout << "Number of Nodes: " << triMesh->getSize(0) << std::endl; // Should be 100
        std::cout << "Number of Triangles: " << triMesh->getSize(2) << std::endl; // Should be 162

        // 4. Access individual faces (triangles)
        JFacePtr firstFace = triMesh->getFaceAt(0);
        if (firstFace) {
            std::cout << "First Face has " << firstFace->getSize(0) << " vertices." << std::endl;
        }

        // 5. Example: Get a fractal-like Sierpinski triangle
        JMeshPtr sierpinski = AllTriMeshGenerator::getSierpinski(3);
        if (sierpinski) {
            std::cout << "Generated Sierpinski triangle (Level 3) with " 
                      << sierpinski->getSize(2) << " triangles." << std::endl;
        }
    } else {
        std::cerr << "Failed to generate Triangle Mesh." << std::endl;
    }

    return 0;
}
