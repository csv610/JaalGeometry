#include <iostream>
#include "AllTetMeshGenerator.hpp"
#include "AllHexMeshGenerator.hpp"

/**
 * Example: Tetrahedral Mesh Generation from a Hexahedral (Structured) Input
 */
int main() {
    // 1. First, create a 2x2x2 hex mesh (1 hex cell, 8 vertices) as input
    int hexdim[3] = {2, 2, 2};
    JMeshPtr hexMesh = AllHexMeshGenerator::getStructuredMesh(hexdim);

    if (hexMesh) {
        std::cout << "Step 1: Generated input hex mesh with " << hexMesh->getSize(3) << " cell." << std::endl;

        // 2. Initialize the Tetrahedral Mesh Generator
        AllTetMeshGenerator tetGen;

        // 3. Convert the hex mesh into a tetrahedral mesh
        // Each hex cell is split into 5 or 6 tetrahedral cells.
        JMeshPtr tetMesh = tetGen.fromHexMesh(hexMesh);

        if (tetMesh) {
            // 4. Print basic tetrahedral mesh information
            std::cout << "Step 2: Successfully converted hex mesh to Tet Mesh." << std::endl;
            std::cout << "Number of Nodes: " << tetMesh->getSize(0) << std::endl; // Should still be 8
            std::cout << "Number of Tetrahedrons: " << tetMesh->getSize(3) << std::endl; // Expect 5 or 6

            // 5. Example: Generate a 3D fractal-like Sierpinski tetrahedron
            JMeshPtr sierpinski = tetGen.getSierpinski(2);
            if (sierpinski) {
                std::cout << "Generated Sierpinski tetrahedron (Level 2) with " 
                          << sierpinski->getSize(3) << " tet cells." << std::endl;
            }
        } else {
            std::cerr << "Failed to generate Tetrahedral Mesh from hex input." << std::endl;
        }
    } else {
        std::cerr << "Failed to generate initial hex mesh for tet conversion." << std::endl;
    }

    return 0;
}
