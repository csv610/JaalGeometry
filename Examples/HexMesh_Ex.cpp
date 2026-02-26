#include <iostream>
#include "AllHexMeshGenerator.hpp"

/**
 * Example: 3D Structured Hexahedral Mesh Generation and Property Queries
 */
int main() {
    // 1. Define dimensions for a 5x5x5 grid of hex cells (6x6x6 vertices)
    int gridim[3] = {6, 6, 6};

    // 2. Generate the structured hexahedral mesh
    // results in (6*6*6) = 216 nodes and (5*5*5) = 125 hex cells
    JMeshPtr hexMesh = AllHexMeshGenerator::getStructuredMesh(gridim);

    if (hexMesh) {
        // 3. Print basic hexahedral mesh information
        std::cout << "Successfully generated a structured Hex Mesh." << std::endl;
        std::cout << "Number of Nodes: " << hexMesh->getSize(0) << std::endl; // Should be 216
        std::cout << "Number of Hexahedrons: " << hexMesh->getSize(3) << std::endl; // Should be 125

        // 4. Access individual cells (hexahedrons)
        JCellPtr firstCell = hexMesh->getCellAt(0);
        if (firstCell) {
            std::cout << "First Cell has " << firstCell->getSize(0) << " vertices." << std::endl;
        }

        // 5. Example: Topological query - check for "singlets"
        // A singlet is a vertex shared by only one cell in the mesh (boundary or isolated)
        size_t singletCount = 0;
        size_t numNodes = hexMesh->getSize(0);
        for (size_t i = 0; i < numNodes; ++i) {
            if (AllHexMeshGenerator::isSinglet(hexMesh->getNodeAt(i))) {
                singletCount++;
            }
        }
        std::cout << "Number of 'singlet' nodes detected: " << singletCount << std::endl;

    } else {
        std::cerr << "Failed to generate Hexahedral Mesh." << std::endl;
    }

    return 0;
}
