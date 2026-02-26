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

## Getting Started

Refer to the `README.md` in the `JaalMesh` directory for technical details on building and using the library.
