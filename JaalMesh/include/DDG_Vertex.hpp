// -----------------------------------------------------------------------------
// libDDG -- Vertex.h
// -----------------------------------------------------------------------------
//
// Vertex stores attributes associated with a mesh edge.  The iterator he
// points to its "outgoing" halfedge.  (See the documentation for a more
// in-depth discussion of the halfedge data structure.)
//

#pragma once

#include "DDG_Vector.hpp"
#include "DDG_Types.hpp"

#include <map>
#include <string>
#include <boost/any.hpp>

namespace DDG
{
class Vertex
{
public:
    HalfEdgeIter he;
    // points to the "outgoing" halfedge

    Vector position;
    // location of vertex in Euclidean 3-space

    int index;
    // unique integer ID in the range 0, ..., nVertices-1

    bool tag;
    // true if vertex is selected by the user; false otherwise

    Vertex() : index(0), tag(false)
    { }

    // returns true if vertex is at boundary
    bool onBoundary( void ) const;

    double area() const;
    // returns the barycentric area associated with this vertex

    Vector normal() const;
    // returns the vertex normal

    bool isIsolated() const;
    // returns true if the vertex is not contained in any face or edge; false otherwise

    int valence() const;
    // returns the number of incident faces / edges

    void toggleTag();

    std::map<std::string, boost::any> attrib;

    // toggle vertex tag

    // returns sum_tip_angles
    double theta() const;

    // zero-form that controls curvature
//   double potential;

    // singularity index
//    double singularity;

//    VertexIter parent;
    double weight;
    double distance;

};
}


