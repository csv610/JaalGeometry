#include <iostream>
#include "AllTriMeshGenerator.hpp"
#include "AllQuadMeshGenerator.hpp"

using namespace Jaal;

/**
 * Example: Converting a Triangle Mesh to a Quadrilateral Mesh
 */
int main() {
    // 1. Generate an initial Triangle Mesh (10x10 vertices)
    JMeshPtr triMesh = AllTriMeshGenerator::getStructuredMesh(10, 10);
    if (!triMesh) {
        std::cerr << "Failed to generate Triangle Mesh." << std::endl;
        return 1;
    }
    std::cout << "Original Triangle Mesh: " << triMesh->getSize(2) << " triangles." << std::endl;

    // 2. Initialize the Quad Mesh Generator with the triangle mesh
    AllQuadMeshGenerator quadGen;
    quadGen.setMesh(triMesh);

    // --- Option A: Triangle Matching (Pairs adjacent triangles) ---
    // This reduces the number of elements but may leave some triangles.
    // Algorithms: GREEDY_MATCHING, EDMONDS_MATCHING, BINARY_TREE_MATCHING
    JMeshPtr matchedQuadMesh = quadGen.getTrianglesMatching(AllQuadMeshGenerator::BINARY_TREE_MATCHING);
    
    if (matchedQuadMesh) {
        std::cout << "
After Triangle Matching (Binary Tree):" << std::endl;
        std::cout << "Number of Quads: " << matchedQuadMesh->getSize(2) << std::endl;
        // Note: Some triangles might still remain if the matching wasn't perfect.
    }

    // --- Option B: Catmull-Clark Subdivision (All-Quad Mesh) ---
    // Every triangle is split into exactly 3 quads. This guarantees an all-quad mesh.
    quadGen.setMesh(triMesh); // Reset to original triangle mesh
    JMeshPtr allQuadMesh = quadGen.getCatmullClarkMesh();

    if (allQuadMesh) {
        std::cout << "
After Catmull-Clark Subdivision:" << std::endl;
        std::cout << "Number of Quads: " << allQuadMesh->getSize(2) << std::endl;
        // Every triangle (162) becomes 3 quads -> 162 * 3 = 486 quads.
    }

    return 0;
}
