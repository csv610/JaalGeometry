#include <iostream>
#include <cassert>
#include "AllTriMeshGenerator.hpp"

void test_structured_tri_mesh() {
    std::cout << "Testing AllTriMeshGenerator::getStructuredMesh..." << std::endl;
    
    // Create a 2x2 grid of vertices, which means 1x1 cells.
    // Each quad cell in a structured grid is typically split into 2 triangles.
    // So a 1x1 grid should have 2 triangles.
    JMeshPtr mesh = AllTriMeshGenerator::getStructuredMesh(2, 2);
    
    assert(mesh != nullptr);
    assert(mesh->getSize(0) == 4); // 4 vertices
    assert(mesh->getSize(2) == 2); // 2 triangles
    
    std::cout << "AllTriMeshGenerator::getStructuredMesh passed!" << std::endl;
}

int main() {
    test_structured_tri_mesh();
    std::cout << "All TriangleMesh tests passed!" << std::endl;
    return 0;
}
