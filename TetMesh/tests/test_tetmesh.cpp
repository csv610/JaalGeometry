#include <iostream>
#include <cassert>
#include "AllTetMeshGenerator.hpp"
#include "AllHexMeshGenerator.hpp"

void test_tet_from_hex_mesh() {
    std::cout << "Testing AllTetMeshGenerator::fromHexMesh..." << std::endl;
    
    // Create a 2x2x2 hex mesh (1 hex cell).
    int hexdim[3] = {2, 2, 2};
    JMeshPtr hexmesh = AllHexMeshGenerator::getStructuredMesh(hexdim);
    assert(hexmesh != nullptr);
    assert(hexmesh->getSize(3) == 1);
    
    // Convert hex cell into tet cells. A hex is usually split into 5 or 6 tets.
    AllTetMeshGenerator tet_gen;
    JMeshPtr tetmesh = tet_gen.fromHexMesh(hexmesh);
    
    assert(tetmesh != nullptr);
    assert(tetmesh->getSize(0) == 8); // 8 vertices
    assert(tetmesh->getSize(3) >= 5); // At least 5 tets
    
    std::cout << "AllTetMeshGenerator::fromHexMesh passed!" << std::endl;
}

int main() {
    test_tet_from_hex_mesh();
    std::cout << "All TetMesh tests passed!" << std::endl;
    return 0;
}
