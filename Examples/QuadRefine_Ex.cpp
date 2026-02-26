#include <iostream>
#include "AllQuadMeshGenerator.hpp"
#include "MeshRefine.hpp"

using namespace Jaal;

/**
 * Example: Quadrilateral Mesh Refinement using JQuadRefiner
 * 
 * JQuadRefiner provides various topological refinement schemes for quad meshes.
 * This example demonstrates the standard 1-to-4 refinement, where each quad
 * is split into four smaller quadrilaterals.
 */
int main() {
    // 1. Generate a simple structured 2x2 quad mesh (3x3 vertices)
    int dim[2] = {3, 3};
    JMeshPtr mesh = AllQuadMeshGenerator::getStructuredMesh(dim);

    if (!mesh) {
        std::cerr << "Failed to generate initial Quad Mesh." << std::endl;
        return 1;
    }
    std::cout << "Initial Mesh: " << mesh->getSize(0) << " nodes, " << mesh->getSize(2) << " quad faces." << std::endl;

    // 2. Initialize the Quad Refiner
    JQuadRefiner refiner;
    refiner.setMesh(mesh);

    // 3. Perform 1-to-4 refinement on all faces in the mesh
    // Scheme 14 or 4 depending on implementation details, usually refineAll(14) or similar.
    // In this API, we can refine all faces using a specific scheme.
    std::cout << "Refining all quad faces (1-to-4 split)..." << std::endl;
    refiner.refineAll(JQuadRefiner::QUAD14); 

    // 4. Print results
    // For a 2x2 grid (4 faces), each face splits into 4 -> 4 * 4 = 16 faces.
    std::cout << "Refined Mesh: " << mesh->getSize(0) << " nodes, " << mesh->getSize(2) << " quad faces." << std::endl;

    return 0;
}
