#include <iostream>
#include <cassert>
#include "AllHexMeshGenerator.hpp"

void test_structured_hex_mesh() {
    std::cout << "Testing AllHexMeshGenerator::getStructuredMesh..." << std::endl;
    
    // Create a 2x2x2 grid of vertices, which means 1x1x1 cells.
    int dim[3] = {2, 2, 2};
    JMeshPtr mesh = AllHexMeshGenerator::getStructuredMesh(dim);
    
    assert(mesh != nullptr);
    assert(mesh->getSize(0) == 8); // 8 vertices
    assert(mesh->getSize(3) == 1); // 1 hex cell
    
    std::cout << "AllHexMeshGenerator::getStructuredMesh passed!" << std::endl;
}

int main() {
    test_structured_hex_mesh();
    std::cout << "All HexMesh tests passed!" << std::endl;
    return 0;
}
