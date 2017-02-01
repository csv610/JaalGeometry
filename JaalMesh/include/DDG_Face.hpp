// -----------------------------------------------------------------------------
// libDDG -- Face.h
// -----------------------------------------------------------------------------
//
// Face stores attributes associated with a mesh edge.  The iterator he points
// to one of its associated halfedges.  (See the documentation for a more
// in-depth discussion of the halfedge data structure.)
//

#pragma once

#include "DDG_Vector.hpp"
#include "DDG_Types.hpp"

#include <string>
#include <boost/any.hpp>

namespace DDG
{
class Face
{
public:
    HalfEdgeIter he;
    // points to one of the halfedges associated with this face

    int index;
    // unique integer ID in the range 0, ..., nFaces-1

//   Vector vector;
    // unit vector tangent to triangle

//   bool tag;
    // auxiliary tag

//  FaceIter parent;
    // parent in tree-cotree decomposition

    Face() : index(0) { }

    bool isBoundary( void ) const;
    // returns true if this face corresponds to a
    // boundary loop; false otherwise

    double area( void ) const;
    // returns the triangle area

    Vector normal( void ) const;
    // returns the unit normal associated with this face; normal
    // orientation is determined by the circulation order of halfedges

    Vector circumcenter( void ) const;
    // returns triangle circumcenter

    Vector barycenter( void ) const;
    // returns triangle barycenter

    Vector dualPoint( void ) const;
    // returns dual vertex from weighted triangulation

    std::map<std::string, boost::any>  attrib;
    Vector vector;

};
}
