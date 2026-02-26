# JaalMesh

JaalMesh is the core geometry processing engine of the JaalGeometry suite. It provides a robust set of C++ templates and classes for mesh data structures, generation, optimization, and topological analysis.

## Key Modules

### 1. Mesh Data Structures and Infrastructure
- **Mesh/MeshEntity**: Core classes for handling vertices, edges, faces, and cells.
- **Attributes**: Generic attribute system for attaching data to mesh entities.
- **MeshImporter/MeshExporter**: Support for various mesh formats (OFF, OBJ, STL, PLY, etc.).
- **JaalMoabConverter**: Integration with the MOAB (Mesh-Oriented datABase) framework.

### 2. Mesh Generation
- **Standard Generators**: `AllHexMeshGenerator`, `AllQuadMeshGenerator`, `AllTetMeshGenerator`, `AllTriMeshGenerator`.
- **Specialized Meshing**:
  - `MSTQuadMesher`: Minimum Spanning Tree-based quad meshing.
  - `DiskQuadMesher` and `RingQuadMesh`: Focused on specific topologies.
  - `SphereHexMesher`: Generators for spherical hex meshes.
- **Advanced Techniques**: `MarchingTriangles`, `DelaunayMesh`, `ConvexHull`.

### 3. Mesh Optimization and Repair
- **Untangling**: `MeshUntangle` and `MeshOptBoundaryLayer` for fixing inverted elements.
- **Smoothing**: `LloydOptimizer`, `LaplacianMeshDeformation`, and `NormalSmoothing`.
- **Quality**: `MeshQuality` for evaluating element aspect ratio, Jacobian, and other metrics.
- **Cleanup**: `QuadCleanUp` and `TriDecimator` for simplifying and refining topology.

### 4. Topology and Graph Analysis
- **Dual Operations**: `MeshDual`, `QuadDual`, and `HexDual` for dual graph construction.
- **Singularity Management**: `MeshSingularityGraph` and `MeshMinSingularity`.
- **Connectivity**: `MeshTopology`, `MeshTraversal`, and `SwapEdges`.

### 5. Discrete Differential Geometry (DDG)
- Extensive support for DDG algorithms, including:
  - `DDG_MeshGeodesics`: Geodesic distance computation.
  - `DDG_MeshFlatten`: Mesh parameterization and flattening.
  - `DDG_DiscreteExteriorCalculus`: Discrete operators (grad, div, curl) on meshes.
  - `DDG_HarmonicBases` and `HarmonicField`.

### 6. Geometry Processing Algorithms
- **Parameterization**: `SurfaceParameterization` and `HarmonicMap`.
- **Segmentation**: `MeshSegmentation` and `MeshShapeDiameterSegmentation`.
- **Boolean Operations**: `MeshBoolean` and `PolyBoolean`.
- **Deformation**: `MeshMeanCurvatureFlow` and `LaplacianMeshDeformation`.

## Dependencies
- **Eigen**: For linear algebra operations.
- **CGAL (Optional)**: Used for certain advanced geometric predicates and mesh boolean operations.
- **MOAB (Optional)**: For large-scale mesh data handling.
- **ShapeOp**: Integrated for geometry optimization.

## Usage
Most of the library is header-only or follows a standard `include/src` structure. To use JaalMesh, include `JaalHeaders.hpp` and link against the compiled library.
