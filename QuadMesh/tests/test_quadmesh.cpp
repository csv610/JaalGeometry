#include <iostream>
#include <cassert>
#include "AllQuadMeshGenerator.hpp"

using namespace Jaal;

void test_structured_quad_mesh() {
    std::cout << "Testing AllQuadMeshGenerator::getStructuredMesh..." << std::endl;
    
    // Create a 2x2 grid of vertices, which means 1x1 cells.
    int dim[2] = {2, 2};
    JMeshPtr mesh = AllQuadMeshGenerator::getStructuredMesh(dim);
    
    assert(mesh != nullptr);
    assert(mesh->getSize(0) == 4); // 4 vertices
    assert(mesh->getSize(2) == 1); // 1 quad
    
    std::cout << "AllQuadMeshGenerator::getStructuredMesh passed!" << std::endl;
}

int main() {
    test_structured_quad_mesh();
    std::cout << "All QuadMesh tests passed!" << std::endl;
    return 0;
}
