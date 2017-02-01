#include "MeshDual.hpp"
#include "MeshRefine.hpp"

///////////////////////////////////////////////////////////////////////////////

int DualChord :: remove_self_intersecting_cycle(Face *face)
{
    JNodeSequence newnodes, newfaces;

    int dim[] = {3,3};
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

void DualChord :: reverse()
{
    size_t nSize = dElements.size();
    for( int i = nSize-1; i > 0; i--) {
        dElements[i].cell = dElements[i-1].cell;
        dElements[i].face = dElements[i-1].face;
        dElements[i].segments[0] = dElements[i-1].segments[1];
        dElements[i].segments[1] = dElements[i-1].segments[0];
        dElements[i].segments[0]->reverse();
        dElements[i].segments[1]->reverse();
    }
    std::reverse( dElements.begin(), dElements.end() );
    dElements.pop_back();   // Otherwise there will be duplication of the last
}

///////////////////////////////////////////////////////////////////////////////

int DualChord :: initialize()
{
    string attribname  = "Steiner";
    topDim             = 0;
    ring               = 0;
    end_at_boundary    = 0;
    simple             = 0;
    self_touching      = 0;
    self_intersecting  = 0;
    dElements.clear();

    selfTouchingNodes.clear();
    selfTouchingEdges.clear();
    selfTouchingFaces.clear();

    selfIntersectingFaces.clear();
    selfIntersectingCells.clear();

    /*
     if( edge->hasAttribute( attribname ) )  {
         cout << "Warning: Edge attribute " << attribname << " may be overwritten "  << endl;
         return 1;
     }

     if( face->hasAttribute( attribname ) ) {
         cout << "Warning: Face attribute " << attribname << " may be overwritten "  << endl;
         return 1;
     }
    */
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

void DualChord :: print()
{
    for( size_t i = 0;  i < dElements.size(); i++)  {
        dElements[i].edge->print();
    }
}
///////////////////////////////////////////////////////////////////////////////

void DualChord :: Delete()
{
    for( size_t i = 0;  i < dElements.size(); i++) {
        dElements[i].Delete();
    }
    initialize();
}

///////////////////////////////////////////////////////////////////////////////

int DualChord :: diceQuads(Mesh *primalMesh, int npieces)
{
    size_t nSize;
    JNodeSequence newnodes;
    JFaceSequence newfaces, faceneighs;

    JEdgeSequence edges;
    getEdges(edges);

    nSize = edges.size();
    for( size_t i = 0; i < nSize; i++) {
        Edge::generate_linear_nodes( edges[i], npieces, primalMesh);
    }

    JFaceSequence faces;
    getFaces(faces);
    nSize = faces.size();
    for( size_t i = 0; i < nSize; i++) {
        QuadRefiner::refine( faces[i], primalMesh);
    }

    nSize = edges.size();
    for( size_t i = 0; i < nSize; i++) {
        edges[i]->deleteAttribute("Steiner");
        edges[i]->setStatus( MeshEntity::REMOVE );
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void DualChord :: shrinkEdges()
{

    Point3D p0, p1, pm;
    Vertex *vm, *v0, *v1, *v2, *v3;

    size_t nSize = dElements.size();
    for( size_t i = 0; i < nSize; i++) {
        Edge *edge = dElements[i].edge;
        if( edge )  {
            v0 = edge->getNodeAt(0);
            v1 = edge->getNodeAt(1);
            Vertex::mid_point( v0, v1, p0, 0.25 );
            Vertex::mid_point( v0, v1, p1, 0.75 );
            v0->setXYZCoords(p0);
            v1->setXYZCoords(p1);
        }
    }

    for( size_t i = 0; i < nSize; i++) {
        Edge *edge = dElements[i].edge;
        if( edge )  {
            if( edge->hasAttribute("Steiner")) {
                edge->getAttribute("Steiner", vm);
                if( vm ) {
                    v0 = edge->getNodeAt(0);
                    v1 = edge->getNodeAt(1);
                    Vertex::mid_point( v0, v1, pm);
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
                    Face::centroid(v0,v1,v2,v3, pm);
                    vm->setXYZCoords(pm);
                }
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

void DualChord :: shrinkColumn()
{
    int  nSize = collapseNodesPair.size()/2;
    Point3D xyz;
    for( int i = 0; i < nSize; i++) {
        Vertex *v0 = collapseNodesPair[2*i];
        Vertex *v1 = collapseNodesPair[2*i+1];
        Vertex::mid_point( v0, v1, xyz, 0.75 );
        v1->setXYZCoords(xyz);
    }
    /*
         size_t nSize = dElements.size();
         Point3D p0, p1;
         for( size_t i = 0; i < nSize; i++) {
              Edge *edge = dElements[i].edge;
              Vertex *v0 = edge->getNodeAt(0);
              Vertex *v1 = edge->getNodeAt(1);
              Vertex::mid_point( v0, v1, p0, 0.25 );
              Vertex::mid_point( v0, v1, p1, 0.75 );
              v0->setXYZCoords(p0);
         }
    */
}

///////////////////////////////////////////////////////////////////////////////

void  DualChord :: addElements( Edge *edge, Face *face)
{
    assert( edge != NULL );

    topDim = 2;

    Point3D xyz;
    Vertex *v0 = NULL, *v1 = NULL, *v2 = NULL;
    if( !edge->hasAttribute("Steiner") ) {
        v0 = Vertex ::newObject();
        edge->getAvgXYZ(xyz);
        v0->setXYZCoords(xyz);
        edge->setAttribute("Steiner", v0);
    }
    edge->getAttribute("Steiner", v0);

    Element delem(edge, face);

    if( face ) {
        if( !face->hasAttribute("Steiner") ) {
            v1 = Vertex ::newObject();
            face->getAvgXYZ(xyz);
            v1->setXYZCoords(xyz);
            face->setAttribute("Steiner", v1);
        }
        face->getAttribute("Steiner", v1);

        Quadrilateral *qface = Quadrilateral::down_cast(face);
        Edge *oedge = qface->getOppositeEdge( edge);
        assert( oedge);

        if( !oedge->hasAttribute("Steiner") ) {
            v2 = Vertex ::newObject();
            oedge->getAvgXYZ(xyz);
            v2->setXYZCoords(xyz);
            oedge->setAttribute("Steiner", v2);
        }
        oedge->getAttribute("Steiner", v2);
        Edge *seg1  = Edge::newObject(v0, v1);
        Edge *seg2  = Edge::newObject(v1, v2);
        delem.segments[0] = seg1;
        delem.segments[1] = seg2;
    }

    dElements.push_back( delem );
}

///////////////////////////////////////////////////////////////////////////////

void  DualChord :: addElements( Face *face, Cell *cell)
{
    assert( face != NULL );
    topDim = 3;

    Vertex *v0, *v1, *v2;
    Point3D p3d;
    if( !face->hasAttribute("Steiner") ) {
        v0 = Vertex ::newObject();
        face->getAvgXYZ(p3d);
        v0->setXYZCoords(p3d);
        face->setAttribute("Steiner", v0);
    }
    face->getAttribute("Steiner", v0);
    assert( v0 );

    Edge *seg1 = NULL, *seg2 = NULL;

    if( cell  ) {
        if( !cell->hasAttribute("Steiner") ) {
            v1 = Vertex ::newObject();
            cell->getAvgXYZ(p3d);
            v1->setXYZCoords(p3d);
            cell->setAttribute("Steiner", v1);
        }
        cell->getAttribute("Steiner", v1);
        assert( v1 );

        Hexahedron  *hex = Hexahedron::down_cast(cell);
        Face *oface = hex->getOppositeFace( face );
        assert( oface);

        if( !oface->hasAttribute("Steiner") ) {
            v2 = Vertex ::newObject();
            oface->getAvgXYZ(p3d);
            v2->setXYZCoords(p3d);
            oface->setAttribute("Steiner", v2);
        }
        oface->getAttribute("Steiner", v2);
        assert( v2 );

        seg1  = Edge::newObject(v0, v1);
        seg2  = Edge::newObject(v1, v2);
    }

    Element delem(face, cell);
    delem.segments[0] = seg1;
    delem.segments[1] = seg2;

    dElements.push_back( delem );
}

///////////////////////////////////////////////////////////////////////////////

void DualChord :: getEdges( JEdgeSequence &seq) const
{
    seq.clear();

    if( dElements.empty() ) return;

    size_t nSize = dElements.size();
    seq.reserve( nSize );
    for( size_t i = 0; i < nSize; i++) {
        Edge *e = dElements[i].edge;
        if( e ) seq.push_back(e);
    }
}

///////////////////////////////////////////////////////////////////////////////

void DualChord :: getFaces( JFaceSequence &seq) const
{
    seq.clear();
    if( dElements.empty() ) return;
    size_t nSize = dElements.size();
    seq.reserve( nSize );
    for( size_t i = 0; i < nSize; i++) {
        Face *f = dElements[i].face;
        if( f ) seq.push_back(f);
    }
}

///////////////////////////////////////////////////////////////////////////////
bool DualChord :: isTouching() const
{
    if( dElements.empty() ) return 0;

    Element elmlast = dElements.back();
    Face  *flast = elmlast.face;
    Edge  *elast = elmlast.edge;

    int pos = flast->getPosOf(elast);
    assert(pos >= 0);

    Edge *side1 = flast->getEdgeAt(pos+1);
    Edge *side2 = flast->getEdgeAt(pos+3);

    int nSize = dElements.size();
    for( int i = 0; i < nSize-1; i++) {
        Face *fpast = dElements[i].face;
        Edge *epast = dElements[i].edge;

        pos = fpast->getPosOf(epast);
        Edge *side3 = fpast->getEdgeAt(pos+1);
        Edge *side4 = fpast->getEdgeAt(pos+3);
        if( side1 == side3 || side1 == side4 ) return 1;
        if( side2 == side3 || side2 == side4 ) return 1;
    }

    return 0;
}


void DualChord :: getCells( JCellSequence &seq) const
{
    seq.clear();

    if( dElements.empty() ) return;

    size_t nSize = dElements.size();
    seq.reserve( nSize );
    for( size_t i = 0; i < nSize; i++) {
        Cell *c = dElements[i].cell;
        if( c ) seq.push_back(c);
    }
}

///////////////////////////////////////////////////////////////////////////////

int DualChord :: isMergeable( Vertex  *v1, Vertex *v2) const
{
    //////////////////////////////////////////////////////////////////////////
    //  Case#     v1           v2            Action              Return value  Priority
    //
    //  1   Internal       Internal       Always merge             1
    //  2   boundary       internal       move towards v1          1
    //  3   Intenral       boundary       move towards v2         -1
    //  4   Rigid          Bound(face)    move towards v1          1
    //  5   Bound(face)    Rigid          move towards v2         -1
    //  6   Rigid          Internal       move towards v1          1
    //  7   Internal       Rigid          move towards v2         -1
    //  8   Default                       Don't move               0
    ///////////////////////////////////////////////////////////////////////////

    JNoImpl();

    // Case 1:
    if( !v1->isBoundary() && !v2->isBoundary() ) return 1;

    int v1rigid = v1->hasAttribute("Rigid");
    int v2rigid = v2->hasAttribute("Rigid");

    if( v1rigid && v2rigid ) return 0;
    exit(0);
    /*
          int gid1, gid2;

         int v1bound = v1->isBoundary();
         int v2bound = v2->isBoundary();
         int geo1    = v1->hasAttribute("GeoDimension");
         int geo2    = v2->hasAttribute("GeoDimension");

         // case 2:
         if( v1rigid && !v2rigid && v2bound && geo2)  {
              v2->getAttribute("GeoDimension", gid2);
              if( gid2 == 2) return  1;
         }

         // case 3:
         if( !v1rigid && v2rigid && v1bound && geo1)  {
              v1->getAttribute("GeoDimension", gid1);
              if( gid1 == 2) return  -1;
         }

         // case 4:
         if( v1rigid && !v2bound) return 1;

         // case 5:
         if( !v1bound && v2rigid) return -1;

         if( v1bound && !v2bound) return  1;

         if( !v1bound && v2bound) return -1;
    */

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

bool DualChord :: isRemovable()
{
    collapseNodesPair.clear();
    if( !isSimple() ) return 0;

    Face *firstFace = getSeedFace();

    Vertex *v0 = NULL, *v1 = NULL, *v2 = NULL, *v3 = NULL;

    v0 = firstFace->getNodeAt(0);
    v1 = firstFace->getNodeAt(2);

    int dir;
    dir = isMergeable(v0,v1);
    if( dir == -1)  std::swap(v0,v1);

    if( dir == 0) {
        v0  = firstFace->getNodeAt(1);
        v1  = firstFace->getNodeAt(3);
        dir = isMergeable(v0,v1);
        if( dir == -1)  std::swap(v0,v1);
        if( dir == 0) {
            collapseNodesPair.clear();
            return 0;
        }
    }

    JCellSequence column;
    getCells(column);

    size_t numFaces;
    if( ring )
        numFaces = column.size();
    else
        numFaces = column.size()+1;

    Face *currface = firstFace;
    for( size_t i = 0; i < numFaces; i++) {
        int dir = isMergeable( v0, v1) ;
        if( dir == 0) {
            collapseNodesPair.clear();
            return 0;
        } else {
            if( dir > 0 ) {
                collapseNodesPair.push_back(v0);
                collapseNodesPair.push_back(v1);
            } else {
                collapseNodesPair.push_back(v1);
                collapseNodesPair.push_back(v0);
            }
        }

        if( i < column.size() ) {
            Hexahedron *hex = Hexahedron::down_cast( column[i] );
            hex->get_parallel_face_diagonal(v0, v1, v2, v3);
            currface = hex->getOppositeFace(currface);
            assert( currface->hasNode(v2) );
            assert( currface->hasNode(v3) );
            v0 = v2;
            v1 = v3;
        }
    }

    return 1;
}


///////////////////////////////////////////////////////////////////////////////
void DualChord :: append(const DualChord  &other)
{
    int nSize = other.dElements.size();
    for( int i = 0; i < nSize; i++)
        dElements.push_back( other.dElements[i] );
}

///////////////////////////////////////////////////////////////////////////////
void QuadDual :: NoIntersection_Refinement( Face *face, JNodeSequence &newnodes,
        JFaceSequence &newfaces)
{
    newnodes.clear();
    newfaces.clear();

    Vertex *v[17];

    Point3D xyz;
    v[0] = Vertex ::newObject();
    face->getAvgXYZ(xyz);
    v[0]->setXYZCoords(xyz);
    newnodes.push_back( v[0] );

    v[1] = face->getNodeAt(0);
    v[2] = face->getNodeAt(1);
    v[3] = face->getNodeAt(2);
    v[4] = face->getNodeAt(3);

    Edge *edge;

    edge = face->getEdgeAt(0);
    if( !edge->hasAttribute("Steiner") ) {
        v[5] = Vertex ::newObject();
        edge->getAvgXYZ(xyz);
        v[5]->setXYZCoords(xyz);
        edge->setAttribute("Steiner", v[5]);
        newnodes.push_back( v[5] );
    }
    edge->getAttribute("Steiner", v[5] );

    edge = face->getEdgeAt(1);
    if( !edge->hasAttribute("Steiner") ) {
        v[6] = Vertex ::newObject();
        edge->getAvgXYZ(xyz);
        v[6]->setXYZCoords(xyz);
        edge->setAttribute("Steiner", v[6]);
        newnodes.push_back( v[6] );
    }
    edge->getAttribute("Steiner", v[6] );

    edge = face->getEdgeAt(2);
    if( !edge->hasAttribute("Steiner") ) {
        v[7] = Vertex ::newObject();
        edge->getAvgXYZ(xyz);
        v[7]->setXYZCoords(xyz);
        edge->setAttribute("Steiner", v[7]);
        newnodes.push_back( v[7] );
    }
    edge->getAttribute("Steiner", v[7] );

    edge = face->getEdgeAt(3);
    if( !edge->hasAttribute("Steiner") ) {
        v[8] = Vertex ::newObject();
        edge->getAvgXYZ(xyz);
        v[8]->setXYZCoords(xyz);
        edge->setAttribute("Steiner", v[8]);
        newnodes.push_back( v[8] );

    }
    edge->getAttribute("Steiner", v[8] );

    v[9]  = Vertex::mid_node(v[0], v[1]);
    v[10] = Vertex::mid_node(v[0], v[2]);
    v[11] = Vertex::mid_node(v[0], v[3]);
    v[12] = Vertex::mid_node(v[0], v[4]);

    v[13]  = Face::centroid(v[0], v[1], v[2] );
    v[14]  = Face::centroid(v[0], v[2], v[3] );
    v[15]  = Face::centroid(v[0], v[3], v[4] );
    v[16]  = Face::centroid(v[0], v[4], v[1] );

    newnodes.push_back( v[9] );
    newnodes.push_back( v[10] );
    newnodes.push_back( v[11] );
    newnodes.push_back( v[12] );
    newnodes.push_back( v[13] );
    newnodes.push_back( v[14] );
    newnodes.push_back( v[15] );
    newnodes.push_back( v[16] );

    newfaces.resize(12);
    newfaces[0] = Quadrilateral::newObject( v[0], v[9], v[13], v[10] );
    newfaces[1] = Quadrilateral::newObject( v[1], v[5], v[13], v[9] );
    newfaces[2] = Quadrilateral::newObject( v[2], v[10], v[13], v[5] );

    newfaces[3] = Quadrilateral::newObject( v[0], v[10], v[14], v[11] );
    newfaces[4] = Quadrilateral::newObject( v[2], v[6], v[14], v[10] );
    newfaces[5] = Quadrilateral::newObject( v[3], v[11], v[14], v[6] );

    newfaces[6] = Quadrilateral::newObject( v[0], v[11], v[15], v[12] );
    newfaces[7] = Quadrilateral::newObject( v[3], v[7], v[15], v[11] );
    newfaces[8] = Quadrilateral::newObject( v[4], v[12], v[15], v[7] );

    newfaces[9]  = Quadrilateral::newObject( v[0], v[12], v[16], v[9] );
    newfaces[10] = Quadrilateral::newObject( v[4], v[8], v[16], v[12] );
    newfaces[11] = Quadrilateral::newObject( v[1], v[9], v[16], v[8] );

}

void QuadDual :: initialize()
{
    assert( pmesh );
    pmesh->getTopology()->collect_edges();
    pmesh->buildRelations(0,1);
    pmesh->buildRelations(0,2);
    pmesh->buildRelations(1,2);
    pmesh->getTopology()->search_boundary();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int  QuadDual :: dice( Edge *edge, int npieces)
{
    DualChord dc;
    getChord( edge, dc );
    if( dc.classify() == DualChord::SIMPLE_CHORD) {
        dc.diceQuads( pmesh, npieces);
        return 0;
    }
    return 1;
}
///////////////////////////////////////////////////////////////////////////////

void QuadDual :: expandChord( Face *thisFace, Edge *thisEdge, DualChord &chord)
{
    if( thisFace == NULL || thisEdge == NULL ) return;

    Edge *nextEdge = NULL;
    Face *nextFace = NULL;
    Quadrilateral *quad = Quadrilateral::down_cast(thisFace);

    // Get the opposite edge in the given face ...
    nextEdge = quad->getOppositeEdge( thisEdge );
    assert( nextEdge );

    if( chord.hasFace(thisFace) ) {
        chord.addIntersection( thisFace );
//        self_intersecting = 1;
//        chord.selfIntersectingFaces.push_back(thisFace);
    }

    chord.addElements( thisEdge, thisFace);

    /*
         if( chord.isTouching() ) {
              chord.self_touching = 1;
         }
    */

    JFaceSequence neighs;
    nextEdge->getRelations( neighs );
    if( neighs.size() == 2 ) {
        if( neighs[0] == thisFace ) nextFace = neighs[1];
        if( neighs[1] == thisFace ) nextFace = neighs[0];
    }

    // WHen there is cyclic chord, we should stop further...
    if( chord.hasEdge( nextEdge) ) {
//        chord.ring =  1;
        return;
    }

    // We have hit the boundary. No more expansion ...
    if( nextFace == NULL ) {
//        chord.end_at_boundary  = 1;
        chord.addElements( nextEdge );
        return;
    }

    // Recursively make progress..
    expandChord( nextFace, nextEdge, chord);
}
///////////////////////////////////////////////////////////////////////////////

int QuadDual :: getChord(Edge *edge, DualChord  &chord)
{
    nCounter = 0;
    chord.initialize();

    if( edge == NULL ) return 1;

    if( !edge->isActive() ) {
        cout << "Warning : Module enterered in inactive edge " << endl;
        return 1;
    }

    assert( pmesh  );  // Primal mesh must be present...

    if( pmesh->getAdjTable(1,2) == 0)
        pmesh->buildRelations(1,2);

    JFaceSequence faceneighs;
    edge->getRelations(faceneighs);

    if( faceneighs.empty() ) {
        cout << "Warning: edge-face relationis are  absent " << endl;
        return 1;
    }

    if( faceneighs.size() > 2 ) {
        cout << "Warning: A Manifold edge cann't detect the chord " << endl;
        return 1;
    }

    DualChord backwardchord;

    Edge *startedge = edge;

    expandChord(faceneighs[0], startedge, chord);

    if( faceneighs.size() == 1 && !chord.isBoundary()) {
        cout << "Error: A boundary chord must end at the boundary " << endl;
        exit(0);
    }

    // We started with an internal edge and while progressing, we hit the
    // boundary, so we must start again in the backward direction. This chord
    // must be a boundary chord.
    if( faceneighs.size() == 2 && chord.isBoundary()) {
        chord.reverse();
        startedge = edge;   // From the beginning edge
        expandChord(faceneighs[1], startedge, chord);
    }

    chord.finalize();

    return 0;
}

//////////////////////////////////////////////////////////////////////////////

void QuadDual :: getChords( JEdgeSequence &edges, vector<DualChord> &chords)
{
    for( size_t i = 0; i < chords.size(); i++)
        chords[i].Delete();
    chords.clear();

    size_t nSize = edges.size();
    if( nSize == 0) return;

    chords.reserve( nSize );

    DualChord ch;
    for( size_t i = 0; i < nSize; i++) {
        int stat = getChord( edges[i], ch);
        if( stat == 0)  chords.push_back(ch);
    }
}
//////////////////////////////////////////////////////////////////////////////

int QuadDual :: getChords(Face *s, vector<DualChord> &chords)
{
    seedFace = s;
    chords.resize(2);
    getChord(seedFace->getEdgeAt(0), chords[0] );
    getChord(seedFace->getEdgeAt(1), chords[1] );
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int QuadDual :: getAllChords( vector<DualChord> &dchords)
{
    if( !dchords.empty() ) {
        LOG4CXX_INFO( pmesh->logger, "Deleting all previously chords in the mesh ");

        for( size_t i = 0; i < dchords.size(); i++)
            dchords[i].Delete();
        dchords.clear();
    }

    if( pmesh->getTopology()->getDimension() != 2 ) {
        LOG4CXX_WARN( pmesh->logger, "Quadrilateral Dual only for surface mesh ");
        return 1;
    }

    LOG4CXX_INFO( pmesh->logger, "Searching all chords in the quad mesh ");

    Edge *seedEdge, *edge = NULL;
    int val = 0;
    size_t numedges = pmesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        edge = pmesh->getEdgeAt(i);
        if( edge->isActive() ) edge->setAttribute("ChordID", val );
    }

    DualChord ch;
    JEdgeSequence chedges;

    int  currID = 1;
    while(1) {
        numedges = pmesh->getSize(1);
        seedEdge = NULL;
        for( size_t i = 0; i < numedges; i++) {
            edge = pmesh->getEdgeAt(i);
            if( edge->isActive() ) {
                edge->getAttribute("ChordID", val );
                if( val == 0) {
                    seedEdge = edge;
                    break;
                }
            }
        }
        if( seedEdge == NULL ) break;
        getChord(seedEdge, ch);
        ch.getEdges( chedges);
        numedges = chedges.size();
        for( size_t i = 0; i < numedges; i++) {
            chedges[i]->setAttribute("ChordID", currID);
        }
        currID++;
        if( size_t(currID) == pmesh->getSize(1) ) break;
        dchords.push_back(ch);
    }

    pmesh->delete_edge_attribute("ChordID");

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void HexDual :: initialize()
{
    assert( pmesh );

    pmesh->getTopology()->collect_faces();
    pmesh->getTopology()->collect_edges();
    /*
    pmesh->buildRelations(0,1);
    pmesh->buildRelations(0,2);
    pmesh->buildRelations(0,3);
    pmesh->buildRelations(1,2);
        */
    pmesh->buildRelations(2,3);
    pmesh->getTopology()->search_boundary();
}
///////////////////////////////////////////////////////////////////////////////

/*
void HexDual::joinSheets( const DualSheet &sheet1, const DualSheet &sheet2,
                     DualSheet &sheet )
{
     sheet = sheet1;
     if( sheet2.faces.size() < 3 ) return;

     Face *face = sheet2.faces[2];
     if(find( sheet1.faces.begin(), sheet1.faces.end(), face ) != sheet1.faces.end() ) {
          return;
     }

     for( size_t i = 3; i < sheet2.points.size(); i++)
          sheet.points.push_front( sheet2.points[i] );

     for( size_t i = 2; i < sheet2.faces.size(); i++)
          sheet.faces.push_front( sheet2.faces[i] );

     for( size_t i = 1; i < sheet2.cells.size(); i++)
          sheet.cells.push_front( sheet2.cells[i] );

     for( size_t i = 1; i < sheet2.sheetfaces.size(); i++)
          sheet.sheetfaces.push_front( sheet2.sheetfaces[i] );
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
void HexDual::joinChords( const DualChord &chord1, const DualChord &chord2,
                     DualChord &chord )
{
     chord = chord1;
     if( chord2.faces.size() < 2 ) return;

     for( size_t i = 1; i < chord2.points.size(); i++)
          chord.points.push_front( chord2.points[i] );

     for( size_t i = 1; i < chord2.faces.size(); i++)
          chord.faces.push_front( chord2.faces[i] );

     for( size_t i = 0; i < chord2.cells.size(); i++)
          chord.cells.push_front( chord2.cells[i] );
}
*/


///////////////////////////////////////////////////////////////////////////////

void  DualSheet::getAffectedNeighbors( JCellSequence &affectedCells )
{
    affectedCells.clear();

    JCellSet neighSet;
    neighSet.clear();

    JCellSequence cellneighs;
    JEdgeSet::const_iterator it;
    for( it = parEdges.begin(); it != parEdges.end(); ++it) {
        Edge *edge = *it;
        Vertex *v0 = edge->getNodeAt(0);
        v0->getRelations( cellneighs );
        if( cellneighs.empty() ) {
            cout << "Warning: Vertex-Cell Relations absent" << endl;
        }
        for( size_t j = 0; j < cellneighs.size(); j++)
            neighSet.insert( cellneighs[j] );

        Vertex *v1 = edge->getNodeAt(1);
        v1->getRelations( cellneighs );
        if( cellneighs.empty() ) {
            cout << "Warning: Vertex-Cell Relations absent" << endl;
        }
        for( size_t j = 0; j < cellneighs.size(); j++)
            neighSet.insert( cellneighs[j] );
    }

    // Eliminate cells which are part of the sheets, as they anyhow will
    // be processed...
    for( size_t i = 0; i < dElements.size(); i++)
        neighSet.erase( dElements[i].cell );

    JCellSet::iterator cit;
    for( cit = neighSet.begin(); cit != neighSet.end(); ++cit)
        affectedCells.push_back(*cit);
}

///////////////////////////////////////////////////////////////////////////////

void DualSheet :: shrink()
{
    //
    // Shrink all the parallel edges by some amount. This is good for
    // visual debugging as it can tell the affect of sheet removal.
    //

    int s0, s1;
    Point3D p0, p1;
    set<Edge*>:: const_iterator it;

    for( it = parEdges.begin(); it != parEdges.end(); ++it) {
        Edge *edge = *it;
        Vertex *v0 = edge->getNodeAt(0);
        Vertex *v1 = edge->getNodeAt(1);
        int    b0  = v0->isBoundary();
        int    b1  = v1->isBoundary();

        // Both are internal nodes...
        if( b0 == 0  && b1 == 0 ) {
            Vertex::mid_point(v0, v1, p0, 0.25);
            Vertex::mid_point(v0, v1, p1, 0.75);
            v0->setXYZCoords(p0);
            v1->setXYZCoords(p1);
        }

        // v0 is free node and v1 is bounded, Shift the free node only.
        if( b0 == 0 && b1 == 1) {
            Vertex::mid_point(v0, v1, p0, 0.25);
            v0->setXYZCoords(p0);
        }

        // v0 is bounded node and v0 is free, Shift the free node only.
        if( b0 == 1 && b1 == 0) {
            Vertex::mid_point(v1, v0, p0, 0.25);
            v1->setXYZCoords(p0);
        }

        // Both are in the boundary, dangerous case....
        if( b0 == 1 && b1 == 1) {
            bool special = 0;
            if( v0->hasAttribute("Steiner") && v1->hasAttribute("Steiner") ) {
                v0->getAttribute("Steiner", s0);
                v1->getAttribute("Steiner", s1);
                if( s0 == 1 && s1 == 2) {
                    Vertex::mid_point(v0, v1, p0, 0.65);
                    v1->setXYZCoords(p0);
                    special = 1;
                }
                if( s0 == 2 && s1 == 1) {
                    Vertex::mid_point(v1, v0, p0, 0.65);
                    v0->setXYZCoords(p0);
                    special = 1;
                }
            }

            if( v0->hasAttribute("Steiner") && !v1->hasAttribute("Steiner") ) {
                Vertex::mid_point(v0, v1, p0, 0.25);
                v0->setXYZCoords(p0);
                special = 1;
            }

            if( !v0->hasAttribute("Steiner") && v1->hasAttribute("Steiner") ) {
                Vertex::mid_point(v1, v0, p0, 0.25);
                v1->setXYZCoords(p0);
                special = 1;
            }
            if( !special ) {
                Vertex::mid_point(v0, v1, p0, 0.50);
                v0->setXYZCoords(p0);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
int HexDual :: removeSheet(  DualSheet &sheet)
{
    size_t numedges, numfaces, numcells;

    assert( pmesh->getAdjTable(0,1) );
    assert( pmesh->getAdjTable(0,2) );
    assert( pmesh->getAdjTable(0,3) );

    JCellSequence cells, cellneighs;
    JFaceSequence faceneighs;
    JEdgeSequence edgeneighs;

    sheet.getCells(cells);
    numcells = cells.size();
    for(size_t i = 0; i < numcells; i++) {
        Cell *cell = cells[i];
        assert( cell->isActive() );
        cell->setStatus( MeshEntity::REMOVE);
    }

    pmesh->buildRelations(0,3);

    JEdgeSequence paredges;
    sheet.getEdges(paredges);
    int s0, s1;

    size_t numParedges = paredges.size();

    size_t nindex =   pmesh->getSize(0);
    for(size_t i = 0; i < numParedges; i++) {
        Vertex *v0 = paredges[i]->getNodeAt(0);
        Vertex *v1 = paredges[i]->getNodeAt(1);

        Vertex *vm = NULL;
        assert( v0->isActive() );
        assert( v1->isActive() );

        int    b0  = v0->isBoundary();
        int    b1  = v1->isBoundary();

        // Both are internal nodes...
        if( b0 == 0  && b1 == 0 ) vm = Vertex::mid_node(v0, v1);

        // v0 is free node and v1 is bounded, Shift the free node only.
        if( b0 == 0 && b1 == 1) vm = v1->getClone();

        // v0 is bounded node and v0 is free, Shift the free node only.
        if( b0 == 1 && b1 == 0) vm = v0->getClone();

        // Both are in the boundary, dangerous case....
        if( b0 == 1 && b1 == 1) {
            bool special = 0;
            if( v0->hasAttribute("Steiner") && v1->hasAttribute("Steiner") ) {
                v0->getAttribute("Steiner", s0);
                v1->getAttribute("Steiner", s1);
                if( s0 == 1 && s1 == 2) {
                    vm =  v0->getClone();
                    special = 1;
                }
                if( s0 == 2 && s1 == 1) {
                    vm =  v1->getClone();
                    special = 1;
                }
            }

            if( v0->hasAttribute("Steiner") && !v1->hasAttribute("Steiner") ) {
                vm = v1->getClone();
                special = 1;
            }

            if( !v0->hasAttribute("Steiner") && v1->hasAttribute("Steiner") ) {
                vm = v0->getClone();
                special = 1;
            }
            if( !special ) {
                vm = Vertex::mid_node(v0, v1);
            }
            vm->setBoundaryMark( paredges[i]->getBoundaryMark() );
        }

        if( vm ) {
            vm->setID( nindex++);
            pmesh->addObject(vm);

            v0->getRelations(cellneighs);
            assert( v0->isActive() );
            numcells = cellneighs.size();
            for( size_t j = 0; j < numcells; j++) {
                Cell *c = cellneighs[j];
                if( c->isActive() ) c->replace(v0,vm);
                assert( !c->hasNode(v0) );
                assert(  c->hasNode(vm) );
                vm->addRelation(c);
            }

            v1->getRelations(cellneighs);
            assert( v1->isActive() );
            numcells = cellneighs.size();
            for( size_t j = 0; j < numcells; j++) {
                Cell *c = cellneighs[j];
                if( c->isActive() ) c->replace(v1,vm);
                assert( !c->hasNode(v0) );
                assert(  c->hasNode(vm) );
                vm->addRelation(c);
            }

            paredges[i]->setStatus( MeshEntity::REMOVE);

            v0->getRelations(faceneighs);
            numfaces = faceneighs.size();
            for( size_t j = 0; j < numfaces; j++)
                faceneighs[j]->setStatus(MeshEntity::REMOVE);

            v0->getRelations(edgeneighs);
            numedges = edgeneighs.size();
            for( size_t j = 0; j < numedges; j++)
                edgeneighs[j]->setStatus(MeshEntity::REMOVE);

            v1->getRelations(faceneighs);
            numfaces = faceneighs.size();
            for( size_t j = 0; j < numfaces; j++)
                faceneighs[j]->setStatus(MeshEntity::REMOVE);

            v1->getRelations(edgeneighs);
            numedges = edgeneighs.size();
            for( size_t j = 0; j < numedges; j++)
                edgeneighs[j]->setStatus(MeshEntity::REMOVE);

            v0->setStatus( MeshEntity::REMOVE);
            v1->setStatus( MeshEntity::REMOVE);
        }
    }

    JFaceSequence faces;
    sheet.getSheetFaces( faces );
    for(size_t i = 0; i < faces.size(); i++)
        faces[i]->setStatus( MeshEntity::REMOVE );
    sheet.clear();

    numcells = pmesh->getSize(3);
    for( size_t i = 0; i < numcells; i++) {
        Cell *c = pmesh->getCellAt(i);
        if( c->isActive() ) {
            c->getEdges(edgeneighs);
            for( size_t i = 0; i < edgeneighs.size(); i++)
                assert( edgeneighs[i]->isActive() );

            c->getFaces(faceneighs);
            for( size_t i = 0; i < faceneighs.size(); i++)
                assert( faceneighs[i]->isActive() );
        }
    }

    initialize();

    cout << "Mesh After : " << endl;
    cout << "Nodes " << pmesh->getSize(0) << endl;
    cout << "Edges " << pmesh->getSize(1) << endl;
    cout << "Faces " << pmesh->getSize(2) << endl;
    cout << "Cell  " << pmesh->getSize(3) << endl;

    cout << "Active Mesh After : " << endl;
    cout << "Nodes " << pmesh->getActiveSize(0) << endl;
    cout << "Edges " << pmesh->getActiveSize(1) << endl;
    cout << "Faces " << pmesh->getActiveSize(2) << endl;
    cout << "Cell  " << pmesh->getActiveSize(3) << endl;

    numcells = pmesh->getSize(3);

    JNodeSequence nodeneighs;
    for( size_t i = 0; i < numcells; i++) {
        Cell *c = pmesh->getCellAt(i);
        if( c->isActive() ) {

            nodeneighs = c->getNodes();
            for( size_t k = 0; k < nodeneighs.size(); k++) {
                assert( nodeneighs[k]->isActive() );
            }

            c->getFaces( faceneighs );
            c->getEdges( edgeneighs );

            for( size_t j = 0; j < faceneighs.size(); j++) {
                assert( faceneighs[j]->isActive() );
                nodeneighs = faceneighs[j]->getNodes();
                for( size_t k = 0; k < nodeneighs.size(); k++) {
                    assert( nodeneighs[k]->isActive() );
                }
            }

            for( size_t j = 0; j < edgeneighs.size(); j++) {
                assert( edgeneighs[j]->isActive() );
                nodeneighs = edgeneighs[j]->getNodes();
                for( size_t k = 0; k < nodeneighs.size(); k++) {
                    assert( nodeneighs[k]->isActive() );
                }
            }
        }
    }

    /*
         numfaces = pmesh->getSize(2);
         for( int i = 0; i < numfaces; i++) {
              Face *f = pmesh->getFaceAt(i);
              if( f->isActive() ) {
                  f->getEdges( edgenighs );
                  for( int j = 0; j < edgeneighs.size(); j++) {
                       assert( edgeneighs[j]->isActive() );
                       nodeneighs = edgeneighs[j]->getNodes();
                       for( int k = 0; k < nodeneighs.size(); k++)
                            assert( nodeneighs[j]->isActive() );
                  }
                 nodeneighs = f->getNodes();
                 for( int k = 0; k < nodeneighs.size(); k++)
                     assert( nodeneighs[j]->isActive() );
         }
    */


    cout << "Done " << endl;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void HexDual :: expandSheet( Cell *thisCell, Face *thisFace, Edge *thisEdge, DualSheet &sheet)
{
    if( thisFace == NULL || thisEdge == NULL ) return;

    if( sheet.hasFace(thisFace) ) {
        sheet.ring =  1;
        return;
    }

    sheet.addFace(thisFace);
    if( thisCell == NULL  || !thisCell->isActive() ) return;

    Hexahedron *hex = Hexahedron::down_cast(thisCell);
    hex->get_topological_parallel_edges(thisEdge, par_edges);

    if( sheet.hasCell( thisCell ) ) {
        sheet.self_intersecting  =  1;
    } else {
        sheet.addElements(thisCell, thisEdge, par_edges[0], par_edges[1], par_edges[2] );
    }

    // Enter the opposite face ...
    Face *nextFace = hex->getOppositeFace( thisFace );
    Cell *nextCell = NULL;
    JCellSequence cellneighs;
    assert( nextFace );

    nextFace->getRelations( cellneighs );
    if( cellneighs.size() == 2 ) {
        if( cellneighs[0] == thisCell ) nextCell = cellneighs[1];
        if( cellneighs[1] == thisCell ) nextCell = cellneighs[0];
    }

    expandSheet( nextCell, nextFace, par_edges[1], sheet);

    // Enter the first side face ...
    JFaceSequence faceneighs;
    hex->getNeighbors( thisEdge, faceneighs);
    assert( faceneighs.size() == 2 );

    nextFace = NULL;
    if( faceneighs[0] == thisFace) nextFace = faceneighs[1];
    if( faceneighs[1] == thisFace) nextFace = faceneighs[0];
    assert( nextFace );

    nextCell = NULL;
    nextFace->getRelations( cellneighs );
    if( cellneighs.size() == 2 ) {
        if( cellneighs[0] == thisCell ) nextCell = cellneighs[1];
        if( cellneighs[1] == thisCell ) nextCell = cellneighs[0];
    }

    expandSheet( nextCell, nextFace, thisEdge, sheet);

    // Enter Second Side face ...
    Quadrilateral *quad = Quadrilateral::down_cast( thisFace );
    thisEdge = quad->getOppositeEdge( thisEdge );

    hex->getNeighbors( thisEdge, faceneighs);
    assert( faceneighs.size() == 2 );

    if( faceneighs[0] == thisFace) nextFace = faceneighs[1];
    if( faceneighs[1] == thisFace) nextFace = faceneighs[0];

    nextFace->getRelations( cellneighs );
    nextCell = NULL;
    if( cellneighs.size() == 2 ) {
        if( cellneighs[0] == thisCell ) nextCell = cellneighs[1];
        if( cellneighs[1] == thisCell ) nextCell = cellneighs[0];
    }

    expandSheet( nextCell, nextFace, thisEdge,  sheet);
}

////////////////////////////////////////////////////////////////////////////////

int HexDual :: getSheets(Face *f, vector<DualSheet> &sheets)
{
    if( !f->isActive() ) return 1;

    if( pmesh->getAdjTable(2,3) == 0 )
        pmesh->buildRelations(2,3);

    seedFace = f;
    sheets.resize(3);  // Keep it three not two. because 3 is general case.

    JCellSequence cellneighs;
    seedFace->getRelations(cellneighs);
    if( cellneighs.empty() ) return 1;

    sheets[0].clear();
    Edge *startedge = f->getEdgeAt(0);
    expandSheet( cellneighs[0], seedFace, startedge, sheets[0]);

    sheets[1].clear();
    startedge = f->getEdgeAt(1);
    expandSheet( cellneighs[0], seedFace, startedge, sheets[1]);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int HexDual :: getSheets(Cell *s, vector<DualSheet> &sheets)
{
    seedCell = s;
    sheets.clear();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void HexDual :: expandChord( Cell *thisCell, Face *thisFace, DualChord &chord)
{
    if( thisCell == NULL || thisFace == NULL ) return;

    Face *nextFace = NULL;
    Hexahedron *hex = Hexahedron::down_cast(thisCell);
    JCellSequence neighs;

    nextFace = hex->getOppositeFace( thisFace );
    assert( nextFace ) ;

    if( chord.hasFace(thisFace) ) {
//       chord.ring  =  1;
        return;
    }

    if( chord.hasCell( thisCell ) ) {
        chord.addIntersection( thisCell );
    }

    chord.addElements(thisFace, thisCell);

    nextFace->getRelations( neighs );
    Cell *nextCell = NULL;
    if( neighs.size() == 2 ) {
        if( neighs[0] == thisCell ) nextCell = neighs[1];
        if( neighs[1] == thisCell ) nextCell = neighs[0];
    }

    if( nextCell == NULL) {
//        chord.end_at_boundary  = 1;
        chord.addElements( nextFace );
        return;
    }

    expandChord( nextCell, nextFace, chord);
}

////////////////////////////////////////////////////////////////////////////////////

int HexDual :: getChord( Face *seedFace, DualChord &chord)
{
    chord.initialize();
    if( seedFace == NULL ) return 1;

    if( pmesh->getAdjTable(2,3) == 0)
        pmesh->buildRelations(2,3);

    JCellSequence cellneighs;
    seedFace->getRelations( cellneighs );
    if( cellneighs.empty() ) return 1;
    if( cellneighs.size() > 2 ) return 1;

    Face *startface = seedFace;
    expandChord( cellneighs[0], seedFace, chord);

    DualChord backwardchord;
    if( cellneighs.size() == 2 && chord.isBoundary() ) {
        chord.reverse();
        expandChord( cellneighs[1], startface, backwardchord);
        chord.append( backwardchord );
    }

    chord.finalize();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////

int HexDual :: getChords( JFaceSequence &faces, vector<DualChord> &chords)
{
    for( size_t i = 0; i < chords.size(); i++)
        chords[i].Delete();
    chords.clear();

    size_t nSize = faces.size();
    if( nSize == 0) return 1;

    chords.reserve( nSize );

    DualChord ch;
    for( size_t i = 0; i < nSize; i++) {
        int stat = getChord( faces[i], ch);
        if( stat == 0)  chords.push_back(ch);
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////

/*
int HexDual :: getChords( Edge *seedEdge, vector<DualChord> &chords)
{
     chords.clear();

     if( !seedEdge->isBoundary() ) {
          cout << "Warning: Edge must be boundary " << endl;
          return 1;
     }

     JFaceSequence faces;
     seedEdge->getRelations(faces);
     DualChord ch;
     for( size_t i = 0; i < faces.size(); i++) {
          if( faces[i]->isBoundary() ) {
               getChord(faces[i], ch);
               chords.push_back(ch);
          }
     }
     return 0;
}
*/
//////////////////////////////////////////////////////////////////////////////////

int HexDual :: getAllChords( vector<DualChord> &dchords)
{
    if( !dchords.empty() ) {
        LOG4CXX_INFO( pmesh->logger, "Deleting all previously chords in the mesh ");

        for( size_t i = 0; i < dchords.size(); i++)
            dchords[i].Delete();
        dchords.clear();
    }

    LOG4CXX_INFO( pmesh->logger, "Searching all chords in the hex mesh ");

    Face *seedFace, *face;

    size_t numfaces = pmesh->getSize(2);
    int val = 0;
    for( size_t i = 0 ; i < numfaces; i++) {
        face = pmesh->getFaceAt(i);
        if( face->isActive() ) face->setAttribute("ChordID", val);
    }

    int currID = 1;

    DualChord ch;
    JFaceSequence chfaces;
    while(1) {
        numfaces = pmesh->getSize(2);
        seedFace = NULL;
        for(size_t i = 0; i < numfaces; i++) {
            face = pmesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("ChordID", val);
                if( val == 0) {
                    seedFace = face;
                    break;
                }
            }
        }
        if( seedFace == NULL) break;

        getChord( seedFace, ch);
        ch.getFaces( chfaces);
        numfaces = chfaces.size();
        for( size_t i = 0; i < numfaces; i++)
            chfaces[i]->setAttribute("ChordID", currID);
        currID++;
        if( currID >= (int)pmesh->getSize(2) ) break;
        dchords.push_back( ch );
    }

    pmesh->delete_face_attribute("ChordID");

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int HexDual :: getAllSheets( vector<DualSheet> &sheets)
{
    sheets.clear();

    Face *seedFace, *face;

    size_t numfaces = pmesh->getSize(2);
    int val = 0;
    for( size_t i = 0 ; i < numfaces; i++) {
        face = pmesh->getFaceAt(i);
        if( face->isActive() ) face->setAttribute("SheetID", val);
    }

    int currID = 1;

    vector<DualSheet> sh;
    JFaceSequence chfaces;
    while(1) {
        numfaces = pmesh->getSize(2);
        seedFace = NULL;
        for(size_t i = 0; i < numfaces; i++) {
            face = pmesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("SheetID", val);
                if( val == 0) {
                    seedFace = face;
                    break;
                }
            }
        }
        if( seedFace == NULL) break;

        getSheets( seedFace, sh);
        for( size_t j = 0; j < sh.size(); j++ ) {
            sh[j].getFaces( chfaces);
            numfaces = chfaces.size();
            for( size_t i = 0; i < numfaces; i++)
                chfaces[i]->setAttribute("SheetID", currID);
            currID++;
            if( currID >= (int)pmesh->getSize(2) ) break;
            sheets.push_back( sh[j] );
        }
    }

    pmesh->delete_face_attribute("SheetID");

    return 0;
}

////////////////////////////////////////////////////////////////////////////////


int HexDual:: getColumn( Face *seedFace,
                         vector<DualSheet> &sheet1, vector<DualSheet> &sheet2)
{
    /*
       Cell  *seedCell;
       Edge * seedEdge;
       DualSheet hingSheet;
       expandSheet( seedCell, seedFace, seedEdge, hingSheet);

      for( size_t i = 0; i  < hingSheet.cells.size(); i++) {
           expandSheet( hingeSheet.cells[i], seedFace, seedEdge, dual);
      }
      for( size_t i = 0; i  < hingSheet.cells.size(); i++) {
           expandSheet( hingeSheet.cells[i], seedFace, seedEdge, dual);
      }
    */
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int HexDual :: getColumn(Face *seedFace, JCellSequence &column)
{
    column.clear();

    JCellSequence cellneighs;
    seedFace->getRelations( cellneighs );

    if( cellneighs.empty() ) return 1;

    DualChord forwardChord;
    expandChord( cellneighs[0], seedFace, forwardChord);

    if( cellneighs.size() > 1  && !forwardChord.isCyclic() ) {
        assert( cellneighs.size() == 2 );
        DualChord backwardChord;
        expandChord( cellneighs[1], seedFace, backwardChord);
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
int HexDual :: removeColumn(DualChord &chord)
{
    if( !chord.isRemovable() ) return 1;

    JCellSequence cells;
    chord.getCells( cells );

    JFaceSequence cellfaces;
    JEdgeSequence celledges;
    for( size_t i = 0; i < cells.size(); i++) {
        cells[i]->setStatus( MeshEntity::REMOVE );
        cells[i]->getFaces( cellfaces );
        for( size_t j = 0; j < cellfaces.size(); j++)
            cellfaces[j]->setStatus( MeshEntity::REMOVE );
        cells[i]->getEdges( celledges );
        for( size_t j = 0; j < celledges.size(); j++)
            celledges[j]->setStatus( MeshEntity::REMOVE );
    }

    for( size_t i = 0; i < chord.collapseNodesPair.size()/2; i++) {
        Vertex *v0 = chord.collapseNodesPair[2*i];
        Vertex *v1 = chord.collapseNodesPair[2*i+1];

        v1->getRelations( cells );
        for(  size_t j = 0; j < cells.size(); j++) {
            cells[j]->getFaces( cellfaces );
            for( size_t k = 0; k < cellfaces.size(); k++)
                cellfaces[k]->setStatus( MeshEntity::REMOVE );
            cells[j]->getEdges( celledges );
            for( size_t k = 0; k < celledges.size(); k++)
                celledges[k]->setStatus( MeshEntity::REMOVE );
        }

        for( size_t j = 0; j < cells.size(); j++) {
            int stat = cells[j]->replace( v1, v0);
            assert( stat == 0);
        }
        v1->setStatus( MeshEntity::REMOVE );
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int DualSheaf :: swapEdge(Vertex *newedgenode)
{
    if( mesh == NULL ) return 1;

    if( !isTight() ) return 1;

    int dir = 0;
    if( newedgenode == NULL )
        dir = 1;

    int nSize = chord1.getSize();

    JFaceSequence faceSeq1, faceSeq2;
    JCellSequence cellSeq1, cellSeq2;

    chord1.getFaces( faceSeq1 );
    chord2.getFaces( faceSeq2 );

    chord1.getCells( cellSeq1 );
    chord2.getCells( cellSeq2 );

    JEdgeSequence comm_edges;
    JFaceSequence comm_faces;

    JNodeSequence nodes(12);
    JNodeSequence hexnodes(8);

    Hexahedron *hex1, *hex2;

    for( int i = 0; i < nSize; i++) {
        Face *face1 = faceSeq1[i];
        Face *face2 = faceSeq2[i];
        Face::get_shared_entities( face1, face2, comm_edges);
        if(comm_edges.size() != 1 ) return 2;

        Hexahedron *cell1 = Hexahedron::down_cast( cellSeq1[i] );
        Hexahedron *cell2 = Hexahedron::down_cast( cellSeq2[i] );

        Cell::get_shared_entities( cell1, cell2, comm_faces);
        if(comm_faces.size() != 1 ) return 2;

        Quadrilateral *comm_quad = Quadrilateral::down_cast( comm_faces[0] );

        nodes[0] = comm_edges[0]->getNodeAt(0);
        nodes[1] = comm_edges[0]->getNodeAt(1);
        nodes[2] = comm_quad->getDiagonalNode( nodes[0] );
        nodes[3] = comm_quad->getDiagonalNode( nodes[1] );

        nodes[4] = cell1->getDiagonalNode( nodes[2] );
        nodes[5] = cell1->getDiagonalNode( nodes[3] );
        nodes[6] = cell1->getDiagonalNode( nodes[0] );
        nodes[7] = cell1->getDiagonalNode( nodes[1] );

        nodes[8]  = cell2->getDiagonalNode( nodes[2] );
        nodes[9]  = cell2->getDiagonalNode( nodes[3] );
        nodes[10] = cell2->getDiagonalNode( nodes[0] );
        nodes[11] = cell2->getDiagonalNode( nodes[1] );

        switch( dir ) {
        case 1:
            hexnodes[0] = nodes[7];
            hexnodes[1] = nodes[3];
            hexnodes[2] = nodes[11];
            hexnodes[3] = nodes[6];
            hexnodes[4] = nodes[4];
            hexnodes[5] = nodes[0];
            hexnodes[6] = nodes[8];
            hexnodes[7] = nodes[5];
            hex1 = Hexahedron::newObject(hexnodes);

            hexnodes[0] = nodes[6];
            hexnodes[1] = nodes[11];
            hexnodes[2] = nodes[10];
            hexnodes[3] = nodes[2];

            hexnodes[4] = nodes[5];
            hexnodes[5] = nodes[8];
            hexnodes[6] = nodes[9];
            hexnodes[7] = nodes[1];
            hex2 = Hexahedron::newObject(hexnodes);
            cell1->setStatus( MeshEntity::REMOVE );
            cell2->setStatus( MeshEntity::REMOVE );
            mesh->addObject( hex1 );
            mesh->addObject( hex2 );
            break;
        case 2:
            hexnodes[0] = nodes[7];
            hexnodes[1] = nodes[3];
            hexnodes[2] = nodes[11];
            hexnodes[3] = nodes[10];
            hexnodes[4] = nodes[4];
            hexnodes[5] = nodes[0];
            hexnodes[6] = nodes[8];
            hexnodes[7] = nodes[9];
            hex1 = Hexahedron::newObject(hexnodes);

            hexnodes[0] = nodes[7];
            hexnodes[1] = nodes[10];
            hexnodes[2] = nodes[2];
            hexnodes[3] = nodes[6];

            hexnodes[4] = nodes[4];
            hexnodes[5] = nodes[9];
            hexnodes[6] = nodes[1];
            hexnodes[7] = nodes[7];
            hex2 = Hexahedron::newObject(hexnodes);
            cell1->setStatus( MeshEntity::REMOVE );
            cell2->setStatus( MeshEntity::REMOVE );
            mesh->addObject( hex1 );
            mesh->addObject( hex2 );
            break;
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////


bool DualSheaf :: isTight()
{
    if(chord1.getSize() != chord2.getSize()  ) return 0;

    if( !chord1.isSimple() ) return 0;
    if( !chord2.isSimple() ) return 0;

    int nSize = chord1.getSize();

    JCellSequence cellSeq1, cellSeq2;
    chord1.getCells( cellSeq1 );
    chord2.getCells( cellSeq2 );

    JFaceSequence comm_faces;
    for( int i = 0; i < nSize; i++) {
        Cell *cell1 = cellSeq1[i];
        Cell *cell2 = cellSeq2[i];
        Cell::get_shared_entities( cell1, cell2, comm_faces);
        if(comm_faces.size() != 1 ) return 0;
    }

    return 1;
}
////////////////////////////////////////////////////////////////////////////////



