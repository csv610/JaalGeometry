#include <iostream>
#include "AllQuadMeshGenerator.hpp"
#include "QuadChord.hpp"

using namespace Jaal;

/**
 * Example: Extracting a Quad Chord from a Quadrilateral Mesh
 * 
 * A Quad Chord is a topological "strip" of quadrilaterals that flows through
 * the mesh, defined by starting from a "seed edge" and following opposite edges.
 */
int main() {
    // 1. Generate a structured 10x10 quad mesh (11x11 vertices)
    int dim[2] = {11, 11};
    JMeshPtr mesh = AllQuadMeshGenerator::getStructuredMesh(dim);

    if (!mesh) {
        std::cerr << "Failed to generate initial Quad Mesh." << std::endl;
        return 1;
    }
    std::cout << "Generated 10x10 Quad Mesh." << std::endl;

    // 2. Pick a "seed edge" to start the chord.
    // Let's pick an edge from the first face.
    JFacePtr startFace = mesh->getFaceAt(0);
    JEdgePtr seedEdge  = startFace->getEdgeAt(0);

    // 3. Initialize the Quad Chord extractor
    JQuadChord chord;
    chord.setMesh(mesh);
    chord.setSeed(seedEdge);

    // 4. Extract and print information about the chord
    JFaceSequence faces = chord.getFaces();
    JEdgeSequence edges = chord.getEdges();

    std::cout << "Chord extracted starting from edge " << seedEdge->getID() << ":" << std::endl;
    std::cout << "Number of faces in the chord: " << faces.size() << std::endl; // For a 10x10 grid, expect 10.
    std::cout << "Number of edges in the chord: " << edges.size() << std::endl;

    // 5. Check if the chord is cyclic (forms a loop)
    if (chord.isCyclic()) {
        std::cout << "The extracted chord is cyclic (forms a closed loop)." << std::endl;
    } else {
        std::cout << "The extracted chord is a boundary-to-boundary strip." << std::endl;
    }

    return 0;
}
