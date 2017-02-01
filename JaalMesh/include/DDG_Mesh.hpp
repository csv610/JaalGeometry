// -----------------------------------------------------------------------------
// libDDG -- Mesh.h
// -----------------------------------------------------------------------------
//
// Mesh represents a polygonal surface mesh using the halfedge data structure.
// It is essentially a large collection of disjoint vertices, edges, and faces
// that are ``glued together'' by halfedges which encode connectivity (see
// the documentation for an illustration).  By construction, the halfedge data
// structure cannot represent nonorientable surfaces or meshes with nonmanifold
// edges.
//
// Mesh elements are referenced using iterators -- common usage of these
// iterators is to either traverse an entire vector of mesh elements:
//
//    // visit all vertices
//    for( VertexIter i = vertices.begin(); i != vertices.end(); i++ )
//    {
//       //...
//    }
//
// or to perform a local traversal over the neighborhood of some mesh element:
//
//    // visit both halfedges of edge e
//    HalfEdgeIter he = e->he;
//    do
//    {
//       // ...
//
//       he = he->flip;
//    }
//    while( he != e->he );
//
// (See Types.h for an explicit definition of iterator types.)
//
// Meshes with boundary are handled by creating an additional face for each
// boundary loop (the method Face::isBoundary() determines whether a given
// face is a boundary loop).  Isolated vertices (i.e., vertiecs not contained
// in any edge or face) reference a dummy halfedge and can be checked via
// the method Vertex::isIsolated().
//

#pragma once

#include <vector>
#include <string>

#include "DDG_HalfEdge.hpp"
#include "DDG_Vertex.hpp"
#include "DDG_Edge.hpp"
#include "DDG_Face.hpp"
#include "DDG_SparseMatrix.hpp"

namespace DDG
{
class Mesh
{
public:

    Mesh();

    Mesh( const Mesh& mesh );

    const Mesh& operator=( const Mesh& mesh );

    // reads a mesh from a Wavefront OBJ file; return value is nonzero
    // only if there was an error
    int read( const std::string& filename );

    // writes a mesh to a Wavefront OBJ file; return value is nonzero
    // only if there was an error
    int write( const std::string& filename ) const;

    // reloads a mesh from disk using the most recent input filename
    bool reload( void );

    // centers around the origin and rescales to have unit radius
    void normalize( void );

    // returns total mesh area
    double area( void ) const;

    // returns mean edge lenght
    double meanEdgeLength( void  ) const;

    // returns Euler Characteristic Number
    int getEulerCharacteristicNumber( void ) const;

    void clear() {
        halfedges.clear();
        vertices.clear();
        edges.clear();
        faces.clear();
        boundaries.clear();
    }

    std::vector<HalfEdge> halfedges;
    std::vector<Vertex>   vertices;
    std::vector<Edge>     edges;
    std::vector<Face>     faces;
    std::vector<Face>     boundaries;
    // storage for mesh elements

protected:
    std::string inputFilename;

    // assigns a unique, 0-based index to each mesh element
    void indexElements( void );
};
}

