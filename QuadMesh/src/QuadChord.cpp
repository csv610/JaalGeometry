#include "QuadDual.hpp"

/////////////////////////////////////////////////////////////////////////////////////
void JQuadChord :: clear()
{
    chordEdges.clear();
    chordFaces.clear();
}
/////////////////////////////////////////////////////////////////////////////////////

bool JQuadChord :: isCyclic() const
{
    if( chordEdges.empty() ) return 0;
    if( chordType == CLOSED_CHORD ) return 1;
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////

void JQuadChord :: setSeed(const JEdgePtr &edge)
{
    if( mesh == nullptr) {
        cout << "Warning : A null pointer to mesh passed " << endl;
        return;
    }

    if( edge == nullptr ) return;
    if( !edge->isActive() ) return;

    JFaceSequence faceneighs;
    JEdge::getRelations(edge, faceneighs);

    if( faceneighs.empty() ) {
        cout << "Warning: edge-face relationis are  absent " << endl;
        return ;
    }

    if( faceneighs.size() > 2 ) {
        cout << "Warning: A Manifold edge cann't detect the chord " << endl;
        return;
    }

    chordType = -1;

    seedEdge = edge;
    expand(seedEdge, faceneighs[0]);

    if( faceneighs.size() == 1) return;

    if( chordType == CLOSED_CHORD) return;

    boost::reverse(chordEdges);
    boost::reverse(chordFaces);

    chordEdges.pop_back();   // Because it wull be reinserted by the next expand.

    seedEdge = edge;
    expand(seedEdge, faceneighs[1]);

    assert( chordType == BOUNDARY_CHORD);
}

///////////////////////////////////////////////////////////////////////////////

void JQuadChord :: expand(const JEdgePtr &thisEdge, const JFacePtr &thisFace)
{
    if( thisFace == nullptr || thisEdge == nullptr ) return;

    JQuadrilateralPtr quad = JQuadrilateral::down_cast(thisFace);
    assert( quad );

    // Get the opposite edge in the given face ...
    JEdgePtr nextEdge = JQuadrilateral::getOppositeEdge( thisFace, thisEdge );
    assert( nextEdge );

    chordEdges.push_back(thisEdge);
    chordFaces.push_back(thisFace);

    if(nextEdge == seedEdge ) {
        chordType = CLOSED_CHORD;
        return;
    }

    JFacePtr nextFace;
    JFaceSequence neighs;
    JEdge::getRelations( nextEdge, neighs );
    if( neighs.size() == 2 ) {
        if( neighs[0] == thisFace ) nextFace = neighs[1];
        if( neighs[1] == thisFace ) nextFace = neighs[0];
    }

    if( nextEdge == chordEdges[0] ) return;

    // We have hit the boundary. No more expansion ...
    if( nextFace == nullptr ) {
        chordType = BOUNDARY_CHORD;
        assert(nextEdge != nullptr);
        chordEdges.push_back(nextEdge);
        return;
    }

    // Recursively make progress..
    expand( nextEdge, nextFace);
}

///////////////////////////////////////////////////////////////////////////////


#ifdef CSV
int JQuadChord :: remove_self_intersecting_cycle( Face *face)
{
    unused_parameter(face);
    JNodeSequence newnodes, newfaces;
    /*
          Quadrilateral::refine(face, dim, newnodes, newfaces);
             assert( newnodes.size() == 5 );
             assert( newfaces.size() == 4 );
             Edge *edge0 = face->getEdgeAt(0);
             Edge *edge1 = face->getEdgeAt(1);
             dice(edge0, 2);
             dice(edge1, 2);
         */

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void JDualChord :: shrink_parallel_edges()
{
    Point3D p0, p1, pm;
    Vertex *vm, *v0, *v1, *v2, *v3;
    Edge *edge;

    size_t nSize = dElements.size();
    for( size_t i = 0; i < nSize; i++) {
        edge = dElements[i].parEdge1;
        if( edge )  {
            v0 = edge->getNodeAt(0);
            v1 = edge->getNodeAt(1);
            Vertex::getMidPoint( v0, v1, p0, 0.25 );
            Vertex::getMidPoint( v0, v1, p1, 0.75 );
            v0->setXYZCoords(p0);
            v1->setXYZCoords(p1);
        }
        edge = dElements[i].parEdge2;
        if( edge )  {
            v0 = edge->getNodeAt(0);
            v1 = edge->getNodeAt(1);
            Vertex::getMidPoint( v0, v1, p0, 0.25 );
            Vertex::getMidPoint( v0, v1, p1, 0.75 );
            v0->setXYZCoords(p0);
            v1->setXYZCoords(p1);
        }
    }

    for( size_t i = 0; i < nSize; i++) {
        edge = dElements[i].parEdge1;
        if( edge )  {
            if( edge->hasAttribute("Steiner")) {
                edge->getAttribute("Steiner", vm);
                if( vm ) {
                    v0 = edge->getNodeAt(0);
                    v1 = edge->getNodeAt(1);
                    Vertex::getMidPoint( v0, v1, pm);
                    vm->setXYZCoords(pm);
                }
            }
        }

        edge = dElements[i].parEdge2;
        if( edge )  {
            if( edge->hasAttribute("Steiner")) {
                edge->getAttribute("Steiner", vm);
                if( vm ) {
                    v0 = edge->getNodeAt(0);
                    v1 = edge->getNodeAt(1);
                    Vertex::getMidPoint( v0, v1, pm);
                    vm->setXYZCoords(pm);
                }
            }
        }

        Face *face = dElements[i].face;
        if( face ) {
            if( face->hasAttribute("Steiner")) {
                face->getAttribute("Steiner", vm);
                if( vm ) {
                    v0 = face->getNodeAt(0);
                    v1 = face->getNodeAt(1);
                    v2 = face->getNodeAt(2);
                    v3 = face->getNodeAt(3);
                    FaceGeometry::getCentroid(v0,v1,v2,v3, pm);
                    vm->setXYZCoords(pm);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

int JDualChord :: delete_parallel_edges()
{
    if( mesh == nullptr ) return 1;

    map<Vertex*, Vertex*> newnode;

    JEdgeSequence oldedges;
    getEdges(oldedges);

    size_t nsize = oldedges.size();
    for(size_t i = 0; i < nsize; i++) {
        Vertex *v0 = oldedges[i]->getNodeAt(0);
        Vertex *v1 = oldedges[i]->getNodeAt(1);
        Vertex *vm = Vertex::getMidNode(v0,v1);
        newnode[v0] = vm;
        newnode[v1] = vm;
        mesh->addObject(vm);
    }

    JFaceSequence oldfaces;
    getFaces(oldfaces);

    nsize = oldfaces.size();
    for( size_t i = 0; i < nsize; i++) {
        oldfaces[i]->setStatus(MeshEntity::REMOVE);
    }

    mesh->buildRelations(0,2);
    set<Face*> updatefaces;
    JFaceSequence faceneighs;

    nsize = oldedges.size();
    for(size_t i = 0; i < nsize; i++) {
        Vertex *v0 = oldedges[i]->getNodeAt(0);
        Vertex *v1 = oldedges[i]->getNodeAt(1);
        JVertex::getRelations(v0, faceneighs);
        copy(faceneighs.begin(), faceneighs.end(), inserter(updatefaces, updatefaces.end()));
        JVertex::getRelations(v1, faceneighs);
        copy(faceneighs.begin(), faceneighs.end(), inserter(updatefaces, updatefaces.end()));
        oldedges[i]->setStatus(MeshEntity::REMOVE);
        v0->setStatus(MeshEntity::REMOVE);
        v1->setStatus(MeshEntity::REMOVE);
    }
    JNodeSequence qnodes(4);

    Face *face;
    foreach_(face, updatefaces) {
        if( face->isActive() ) {
            Vertex *v0 = face->getNodeAt(0);
            Vertex *v1 = face->getNodeAt(1);
            Vertex *v2 = face->getNodeAt(2);
            Vertex *v3 = face->getNodeAt(3);

            qnodes[0] = v0;
            if( newnode.find(v0) != newnode.end() )
                qnodes[0] = newnode[v0];

            qnodes[1] = v1;
            if( newnode.find(v1) != newnode.end() )
                qnodes[1] = newnode[v1];

            qnodes[2] = v2;
            if( newnode.find(v2) != newnode.end() )
                qnodes[2] = newnode[v2];

            qnodes[3] = v3;
            if( newnode.find(v3) != newnode.end() )
                qnodes[3] = newnode[v3];
        }
        face->setStatus(MeshEntity::REMOVE);
        face = Quadrilateral::newObject(qnodes);
        mesh->addObject(face);
    }
    mesh->pruneAll();
    delete_dual_segments();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void JDualChord :: shrinkColumn()
{
    int  nSize = collapseNodesPair.size()/2;
    Point3D xyz;
    for( int i = 0; i < nSize; i++) {
        Vertex *v0 = collapseNodesPair[2*i];
        Vertex *v1 = collapseNodesPair[2*i+1];
        Vertex::getMidPoint( v0, v1, xyz, 0.75 );
        v1->setXYZCoords(xyz);
    }
    /*
         size_t nSize = dElements.size();
         Point3D p0, p1;
         for( size_t i = 0; i < nSize; i++) {
              Edge *edge = dElements[i].edge;
              Vertex *v0 = edge->getNodeAt(0);
              Vertex *v1 = edge->getNodeAt(1);
              Vertex::getMidPoint( v0, v1, p0, 0.25 );
              Vertex::getMidPoint( v0, v1, p1, 0.75 );
              v0->setXYZCoords(p0);
         }
    */
}

///////////////////////////////////////////////////////////////////////////////

#endif

