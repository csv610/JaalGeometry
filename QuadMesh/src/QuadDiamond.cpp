#include "Diamond.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////
bool
JDiamond::isDiamond(const JFacePtr &face, int &pos, int type)
{
    pos = -1;
    const JNodePtr &v0 = face->getNodeAt(0);
    const JNodePtr &v1 = face->getNodeAt(1);
    const JNodePtr &v2 = face->getNodeAt(2);
    const JNodePtr &v3 = face->getNodeAt(3);

    int d0 = v0->getNumRelations(2);
    int d1 = v1->getNumRelations(2);
    int d2 = v2->getNumRelations(2);
    int d3 = v3->getNumRelations(2);

    assert(d0 > 0 && d1 > 0 && d2 > 0 && d3 > 0);

    /*
         // Boundary Cases ...
         if (v0->isBoundary() || v2->isBoundary()) {
              if( d0 <= v0->get_ideal_face_degree( Face::QUADRILATERAL ) ) return 0;
              if( d2 <= v2->get_ideal_face_degree( Face::QUADRILATERAL ) ) return 0;
              if (!v1->isBoundary() && !v3->isBoundary()) {
                   if (d1 == 3 && d3 == 3) {
                        pos = 1;
                        return 1;
                   }
              }
         }

         if (v1->isBoundary() || v3->isBoundary()) {
              if( d1 <= v1->get_ideal_face_degree( Face::QUADRILATERAL ) ) return 0;
              if( d3 <= v3->get_ideal_face_degree( Face::QUADRILATERAL ) ) return 0;

              if (!v0->isBoundary() && !v2->isBoundary()) {
                   if ((d0 == 3 && d2 == 3)) {
                        pos = 0;
                        return 1;
                   }
              }
         }
    */

    if (v0->isBoundary() || v1->isBoundary() || v2->isBoundary() || v3->isBoundary() ) return 0;

    if( type == 0) {
        if( d0 != 4) return 1;
        if( d1 != 4) return 1;
        if( d2 != 4) return 1;
        if( d3 != 4) return 1;
    }

    if( type == 33 ) {
        if ((d0 == 3 && d2 == 3)) {
            pos = 0;
            return 1;
        }

        if (d1 == 3 && d3 == 3) {
            pos = 1;
            return 1;
        }
        return 0;
    }

    if( type == 34 ) {
        if ( (d0 == 3 && d2 == 4) || ( d0 == 4 && d2 == 3)) {
            if( d1 < 4 || d3 < 4)  return 0;
            pos = 0;
            return 1;
        }

        if ( (d1 == 4 && d3 == 3) || (d1 == 3 && d3 == 4)) {
            if( d0 < 4 || d2 < 4)  return 0;
            pos = 1;
            return 1;
        }
        return 0;
    }

    if( type == 55 ) {
        if (d1 == 5 && d3 == 5) {
            if ((d0 == 3 && d2 == 4) || (d0 == 4 && d2 == 3)) {
                pos = 0;
                return 1;
            }
        }

        if (d0 == 5 && d2 == 5) {
            if ((d1 == 3 && d3 == 4) || (d1 == 4 && d3 == 3)) {
                pos = 1;
                return 1;
            }
        }
        return 0;
    }

    if( type == 3333 ) {
        if (d0 == 3 && d1 == 3 && d2 == 3 && d3 == 3 ) return 1;
    }

    if( type == 5555 ) {
        if (d0 == 5 && d1 == 5 && d2 == 5 && d3 == 5 ) return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JFaceSequence
JDiamond::getDiamonds(int type)
{
    JFaceSequence diamonds;
    if( mesh == nullptr) return diamonds;

    mesh->buildRelations(0, 2);

    int pos;
    size_t numfaces = mesh->getSize(2);
    for (size_t iface = 0; iface < numfaces; iface++) {
        JFacePtr face = mesh->getFaceAt(iface);
        if ( face->isActive() && JDiamond::isDiamond(face, pos, type)) {
            diamonds.push_back(face );
        }
    }
    return diamonds;
}

///////////////////////////////////////////////////////////////////////////////

bool
JFaceClose::isSafe() const
{
    if (!face->isActive()) return 0;

    if (!vertex0->isActive()) return 0;
    if (!vertex2->isActive()) return 0;

    if (vertex0->isBoundary()) return 0;
    if (vertex2->isBoundary()) return 0;

    int vpos = face->getPosOf(vertex0);
    assert(vpos >= 0);

    if (face->getNodeAt(vpos + 2) != vertex2) {
        cout << "Warning: Face-open requires opposite vertices " << endl;
        return 0;
    }

    // Make sure that any neigbour do not have both the vertex0 and vertex2
    JFaceSequence neighs;
    JNode::getRelations(vertex0, neighs );
    int nSize = neighs.size();
    for (int  j = 0; j < nSize; j++) {
        if (neighs[j] != face) {
            int val0 = neighs[j]->hasNode(vertex0);
            int val1 = neighs[j]->hasNode(vertex2);
            if (val0 + val1 == 2) return 0;
        }
    }

    JNode::getRelations( vertex2, neighs );
    nSize = neighs.size();
    for (int j = 0; j < nSize; j++) {
        if (neighs[j] != face) {
            int val0 = neighs[j]->hasNode(vertex0);
            int val1 = neighs[j]->hasNode(vertex2);
            if (val0 + val1 == 2) return 0;
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
JFaceClose::build()
{
    replaceNode = nullptr;

    if (!isSafe()) return 1;

    JNodePtr v0 = vertex0;
    JNodePtr v2 = vertex2;

    // Back up Coordinates ...
    const Point3D &p0 = v0->getXYZCoords();
    const Point3D &p2 = v2->getXYZCoords();

    // Temporary assign the new coordinate and check what happen to
    // other neighbouring faces, if someone get inverted, then we
    // should apply the face close operation. We may visit the same
    // face again, after applying shape optimization. and hope that
    // it will not invert then.

    Point3D p3d;

    switch( node_placement_policy) {
    case MID_POINT_POSITION:
        p3d = JNodeGeometry::getMidPoint(v0, v2);
        vertex0->setXYZCoords(p3d);
        vertex2->setXYZCoords(p3d);
        break;
    case FIRST_POINT_POSITION:
        vertex0->setXYZCoords(p0);
        vertex2->setXYZCoords(p0);
        break;
    case SECOND_POINT_POSITION:
        vertex0->setXYZCoords(p2);
        vertex2->setXYZCoords(p2);
        break;
    }

    int pass = 1;
    JFaceSequence vfaces;

    if( check_inversion ) {
        if( pass ) {
            JNode::getRelations( vertex0, vfaces );
            for( size_t i = 0; i < vfaces.size(); i++) {
                if( vfaces[i] !=  face ) {
                    if (!JFaceGeometry::isSimple( vfaces[i])) {
                        pass = 0;
                        break;
                    }
                }
            }
        }

        if( pass ) {
            JNode::getRelations( vertex2, vfaces );
            for( size_t i = 0; i < vfaces.size(); i++) {
                if( vfaces[i] !=  face ) {
                    if (!JFaceGeometry::isSimple( vfaces[i]) ) {
                        pass = 0;
                        break;
                    }
                }
            }
        }

        if( !pass ) {
            vertex0->setXYZCoords( p0 );
            vertex2->setXYZCoords( p2 );
            return 1;
        }
    }

    replaceNode = JNode::newObject();
    replaceNode->setXYZCoords(p3d);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
JFaceClose::commit()
{

#ifdef CSV
    if (replaceNode == nullptr) return 1;

    mesh->addObject(replaceNode);

    JFaceSequence vfaces;
    JNode::getRelations( vertex0, vfaces );
    size_t nSize = vfaces.size();

    for (size_t i = 0; i < nSize; i++) {
        Face *fn = vfaces[i];
        if ( fn != face) {
            mesh->deactivate(fn);
            fn->replace(vertex0, replaceNode);
            mesh->reactivate(fn);
        }
    }

    JNode::getRelations( vertex2, vfaces );
    nSize = vfaces.size();
    for (size_t i = 0; i < vfaces.size(); i++) {
        Face *fn = vfaces[i];
        if ( fn != face) {
            mesh->deactivate(fn);
            fn->replace(vertex2, replaceNode);
            mesh->reactivate(fn);
        }
    }

    // Two nodes and face go away from the mesh..
    mesh->remove(face);
    mesh->remove(vertex0);
    mesh->remove(vertex2);
#endif

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
