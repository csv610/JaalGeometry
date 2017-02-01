#include "MeshDual.hpp"
#include "MeshRefine.hpp"

///////////////////////////////////////////////////////////////////////////////

#ifdef DFDFD
int JDualChord :: remove_self_intersecting_cycle( const JFacePtr &face)
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


void DualChord :: reverse()
{
    size_t nSize = dElements.size();

    switch(topDim) {
    case 3:
        for( int i = nSize-1; i > 0; i--)
            dElements[i].cell = dElements[i-1].cell;
        break;
    case 2:
        for( int i = nSize-1; i > 0; i--)
            dElements[i].face = dElements[i-1].face;
        break;
    }

    for( int i = nSize-1; i > 0; i--) {
        dElements[i].segments[0] = dElements[i-1].segments[1];
        dElements[i].segments[1] = dElements[i-1].segments[0];
        dElements[i].segments[0]->reverse();
        dElements[i].segments[1]->reverse();
    }

    std::reverse( dElements.begin(), dElements.end() );
    dElements.pop_back();   // Otherwise there will be duplication of the last

    if( end_at_boundary ) {
        start_at_boundary = 1;
        end_at_boundary   = 0;
    }
}

int JDualChord :: initialize()
{
    string attribname  = "Steiner";

    ring               = 0;
    topDim             = 0;
    end_at_boundary    = 0;
    start_at_boundary  = 0;
    simple             = 0;
    self_touching      = 0;
    self_intersecting  = 0;

    mesh     = NULL;
    seedEdge = NULL;
    seedFace = NULL;

    dElements.clear();

    selfTouchingNodes.clear();
    selfTouchingEdges.clear();
    selfTouchingFaces.clear();

    selfIntersectingFaces.clear();
    selfIntersectingCells.clear();
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

void JDualChord :: print()
{
    for( size_t i = 0;  i < dElements.size(); i++)  {
//        dElements[i].edge->print();
        dElements[i].face->print();
    }
}
///////////////////////////////////////////////////////////////////////////////

void JDualChord :: delete_dual_segments()
{
    for( size_t i = 0;  i < dElements.size(); i++) {
        Edge *e1 = dElements[i].segments[0];
        Edge *e2 = dElements[i].segments[1];
        if( e1 ) delete e1;
        if( e2 ) delete e2;
        dElements[i].segments[0] = NULL;
        dElements[i].segments[1] = NULL;
    }
    initialize();
}

///////////////////////////////////////////////////////////////////////////////

int DualChord :: diceQuads(int npieces)
{
    unused_parameter(npieces);
    /*
         if( mesh == NULL ) return 1;

         size_t nSize;
         JNodeSequence newnodes;
         JFaceSequence newfaces, faceneighs;

         JEdgeSequence edges;
         getEdges(edges);

         //We shall create new nodes on some subset of edges.
         mesh->delete_edge_attribute("Steiner");

         nSize = edges.size();
         for( size_t i = 0; i < nSize; i++) {
              Edge::generate_linear_nodes( edges[i], npieces, mesh);
         }

         JFaceSequence faces;
         getFaces(faces);
         nSize = faces.size();
         for( size_t i = 0; i < nSize; i++) {
              QuadRefiner::refine( faces[i], mesh);
         }

         nSize = edges.size();
         for( size_t i = 0; i < nSize; i++) {
              edges[i]->deleteAttribute("Steiner");
              edges[i]->setStatus( MeshEntity::REMOVE );
         }
         delete_dual_segments();
    */

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void JDualChord :: shrink_parallel_edges()
{
    Point3D p0, p1, pm;
    JNodePtr vm, v0, v1, v2, v3;
    JEdgePtr edge;

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

        JFacePtr face = dElements[i].face;
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
    if( mesh == NULL ) return 1;

    map<JNodePtr, JNodePtr> newnode;

    JEdgeSequence oldedges;
    getEdges(oldedges);

    size_t nsize = oldedges.size();
    for(size_t i = 0; i < nsize; i++) {
        JNodePtr v0 = oldedges[i]->getNodeAt(0);
        JNodePtr v1 = oldedges[i]->getNodeAt(1);
        JNodeptr vm = Vertex::getMidNode(v0,v1);
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
    JFaceSet updatefaces;
    JFaceSequence faceneighs;

    nsize = oldedges.size();
    for(size_t i = 0; i < nsize; i++) {
        JNodePtr v0 = oldedges[i]->getNodeAt(0);
        JNodePtr v1 = oldedges[i]->getNodeAt(1);
        v0->getRelations(faceneighs);
        boost::copy(faceneighs. inserter(updatefaces, updatefaces.end()));
        v1->getRelations(faceneighs);
        boost::copy(faceneighs. inserter(updatefaces, updatefaces.end()));
        oldedges[i]->setStatus(MeshEntity::REMOVE);
        v0->setStatus(MeshEntity::REMOVE);
        v1->setStatus(MeshEntity::REMOVE);
    }
    JNodeSequence qnodes(4);

    JFacePtr face;
    foreach_(face, updatefaces) {
        if( face->isActive() ) {
            JNodePtr v0 = face->getNodeAt(0);
            JNodePtr v1 = face->getNodeAt(1);
            JNodePtr v2 = face->getNodeAt(2);
            JNodePtr v3 = face->getNodeAt(3);

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
        JNodePtr v0 = collapseNodesPair[2*i];
        JNodePtr v1 = collapseNodesPair[2*i+1];
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

void JDualChord :: addElements( JEdgePtr edge, JFacePtr face)
{
    assert( edge != nullptr );
    if( seedEdge == nullptr ) seedEdge = edge;

    if( face == NULL ) end_at_boundary = 1;

    topDim = 2;

    Point3D xyz;
    JNodePtr v0, v1, v2;
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

        JQuadrilateralPtr qface = Quadrilateral::down_cast(face);
        JEdgePtr oedge = qface->getOppositeEdge(edge);
        assert( oedge);

        if( !oedge->hasAttribute("Steiner") ) {
            v2 = Vertex ::newObject();
            oedge->getAvgXYZ(xyz);
            v2->setXYZCoords(xyz);
            oedge->setAttribute("Steiner", v2);
        }
        oedge->getAttribute("Steiner", v2);

        JEdgePtr seg1  = Edge::newObject(v0, v1);
        JEdgePtr seg2  = Edge::newObject(v1, v2);
        delem.segments[0] = seg1;
        delem.segments[1] = seg2;
    }

    dElements.push_back( delem );
}

///////////////////////////////////////////////////////////////////////////////

void JDualChord :: addElements( JFacePtr f, JEdgePtr edge1, JCellPtr cell)
{
    JQuadrilateralPtr qface = Quadrilateral::down_cast(f);
    assert( qface );

    topDim = 3;

    if( cell == nullptr) end_at_boundary = 1;

    JNodePtr v0, v1, v2;

    Point3D p3d;
    if( !qface->hasAttribute("Steiner") ) {
        v0 = Vertex ::newObject();
        qface->getAvgXYZ(p3d);
        v0->setXYZCoords(p3d);
        qface->setAttribute("Steiner", v0);
    }
    qface->getAttribute("Steiner", v0);
    assert( v0 );

    JEdgePtr edge2;
    if( edge1 )  {
        if( !edge1->hasAttribute("Steiner") ) {
            v0 = Vertex ::newObject();
            edge1->getAvgXYZ(p3d);
            v0->setXYZCoords(p3d);
            edge1->setAttribute("Steiner", v0);
        }

        edge2 = qface->getOppositeEdge(edge1);
        if( !edge2->hasAttribute("Steiner") ) {
            v0 = Vertex ::newObject();
            edge2->getAvgXYZ(p3d);
            v0->setXYZCoords(p3d);
            edge2->setAttribute("Steiner", v0);
        }
    }

    JEdgePtr seg1, seg2;
    if( cell  ) {
        if( !cell->hasAttribute("Steiner") ) {
            v1 = Vertex ::newObject();
            cell->getAvgXYZ(p3d);
            v1->setXYZCoords(p3d);
            cell->setAttribute("Steiner", v1);
        }
        cell->getAttribute("Steiner", v1);
        assert( v1 );

        JHexahedronPtr hex = Hexahedron::down_cast(cell);
        JFacePtr oface = hex->getOppositeFace( qface );
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


    Element delem(qface, edge1, edge2, cell);
    delem.segments[0] = seg1;
    delem.segments[1] = seg2;

    dElements.push_back( delem );
}

///////////////////////////////////////////////////////////////////////////////

void JDualChord :: getNodes( JNodeSequence &seq) const
{
    seq.clear();

    if( dElements.empty() ) return;

    JNodeSet ndset;

    size_t nSize = dElements.size();
    JNodePtr vtx;
    for( size_t i = 0; i < nSize; i++) {
        JEdgePtr e1  = dElements[i].parEdge1;
        if( e1 )  {
            e1->getAttribute("Steiner", vtx);
            ndset.insert(vtx);
        }

        JEdgePtr e2  = dElements[i].parEdge2;
        if( e2 )  {
            e2->getAttribute("Steiner", vtx);
            ndset.insert(vtx);
        }

        JFacePtr f  = dElements[i].face;
        if( f )  {
            f->getAttribute("Steiner", vtx);
            ndset.insert(vtx);
        }

        JCellPtr c  = dElements[i].cell;
        if( c )  {
            c->getAttribute("Steiner", vtx);
            ndset.insert(vtx);
        }
    }

    if( !ndset.empty() ) {
        seq.resize( ndset.size() );
        boost::copy( ndset.seq.begin() ) ;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JDualChord :: get_parallel_nodes( JNodeSequence &seq) const
{
    seq.clear();

    if( dElements.empty() ) return;

    JNodeSet ndset;

    size_t nSize = dElements.size();
    for( size_t i = 0; i < nSize; i++) {
        JEdgePtr e1  = dElements[i].parEdge1;
        if( e1 )  {
            ndset.insert( e1->getNodeAt(0) );
            ndset.insert( e1->getNodeAt(1) );
        }

        JEdgePtr e2  = dElements[i].parEdge2;
        if( e2 )  {
            ndset.insert( e2->getNodeAt(0) );
            ndset.insert( e2->getNodeAt(1) );
        }
    }

    if( !ndset.empty() ) {
        seq.resize( ndset.size() );
        boost::copy( ndset, seq.begin() ) ;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JDualChord :: getEdges( JEdgeSequence &seq) const
{
    seq.clear();

    if( dElements.empty() ) return;

    size_t nSize = dElements.size();
    seq.reserve( nSize );
    for( size_t i = 0; i < nSize; i++) {
        JEdgePtr e1 = dElements[i].parEdge1;
        if( e1 ) seq.push_back(e1);

        JEdgePtr e2 = dElements[i].parEdge2;
        if( e2 ) seq.push_back(e2);
    }
}

///////////////////////////////////////////////////////////////////////////////
void JDualChord :: getFaces( JFaceSequence &seq) const
{
    /*
        seq.clear();

        if( dElements.empty() ) return;
        size_t nSize = dElements.size();
        seq.reserve( nSize );
        for( size_t i = 0; i < nSize; i++) {
            JFacePtr f = dElements[i].face;
            if( f ) seq.push_back(f);
        }
    */
}
///////////////////////////////////////////////////////////////////////////////

void JDualChord :: getCells( JCellSequence &seq) const
{
    /*
        seq.clear();

        if( dElements.empty() ) return;

        size_t nSize = dElements.size();
        seq.reserve( nSize );
        for( size_t i = 0; i < nSize; i++) {
            JCellPtr c = dElements[i].cell;
            if( c ) seq.push_back(c);
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

bool JDualChord :: isTouching() const
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
#endif

///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
int JDualChord :: isMergeable( const JNodePtr &v1, const JNodePtr &v2) const
{
    /*
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

        NoImpl();

        // Case 1:
        if( !v1->isBoundary() && !v2->isBoundary() ) return 1;

        int v1rigid = v1->hasAttribute("Rigid");
        int v2rigid = v2->hasAttribute("Rigid");

        if( v1rigid && v2rigid ) return 0;
        exit(0);
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
#endif
///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
bool JDualChord :: isRemovable()
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
#endif


///////////////////////////////////////////////////////////////////////////////

