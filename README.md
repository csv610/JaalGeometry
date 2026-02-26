# JaalGeometry

JaalGeometry is a comprehensive C++ library suite for geometry processing and mesh generation. It focuses on advanced algorithms for quad and hex meshing, topological optimization, and discrete differential geometry.

## Project Structure

- **JaalMesh**: The core engine containing data structures and algorithms for mesh generation, optimization, and analysis. This is the primary active library.
- **JaalGraphics (Legacy)**: A Qt-based graphical user interface for visualizing and interacting with meshes. *Note: This component is currently deprecated as the project moves away from Qt dependencies.*
- **Demo**: A collection of video demonstrations showing the meshing algorithms (e.g., MST, quad meshing) in action.

## Core Capabilities (JaalMesh)

- **Mesh Generation**: Support for Hex, Quad, Tet, and Tri meshes.
- **Specialized Meshing**: MST-based quad meshing, Sphere-Hex meshing, and Advancing Front techniques.
- **Topological Operations**: Dual graph construction, chordal analysis, and singularity management.
- **Discrete Differential Geometry**: Implementations for geodesics, curvature flow, and harmonic maps.
- **Mesh Optimization**: Untangling, smoothing, and quality-driven refinement.

## Examples

Illustrative C++ examples are provided in the `Examples/` directory to demonstrate the usage of each module:

- **`CompGeom_Ex.cpp`**: Basic geometric queries (point-in-polygon, area calculation).
- **`TriangleMesh_Ex.cpp`**: Generation of structured and fractal-like (Sierpinski) triangle meshes.
- **`QuadMesh_Ex.cpp`**: Generation of structured quadrilateral meshes and common test shapes (Schneider's Pyramid).
- **`TetMesh_Ex.cpp`**: Tetrahedral mesh generation from structured hexahedral inputs.
- **`HexMesh_Ex.cpp`**: 3D structured hexahedral mesh generation and basic topological queries.
- **`DelaunayMesh_Ex.cpp`**: 2D Delaunay triangulation from a custom set of points.
- **`AlphaMSTQuadMesh_Ex.cpp`**: Local quadrilateral remeshing within a user-defined circular region.
- **`Tri2Quad_Ex.cpp`**: Conversion of a triangle mesh into a quadrilateral mesh using matching and subdivision.
- **`QuadChord_Ex.cpp`**: Extracting topological "strips" (chords) from a quadrilateral mesh starting from a seed edge.
- **`QuadRefine_Ex.cpp`**: Topological refinement of quadrilateral meshes (e.g., 1-to-4 splitting).

## Getting Started

Refer to the `README.md` in the `JaalMesh` directory for technical details on building and using the library.
