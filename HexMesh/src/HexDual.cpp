#include "MeshDual.hpp"
#include "MeshRefine.hpp"

///////////////////////////////////////////////////////////////////////////////
#ifdef CSV

void HexDual :: initialize()
{
    if( pmesh == NULL ) return;

//    pmesh->getTopology()->collect_faces();
//    pmesh->getTopology()->collect_edges();

    pmesh->buildRelations(2,3);
    pmesh->getTopology()->search_boundary();
}

///////////////////////////////////////////////////////////////////////////////

void HexDual :: expandChord( Cell *thisCell, Face *thisFace, Edge *thisEdge, DualChord *chord)
{
    if( thisCell == NULL || thisFace == NULL ) return;

    Face *nextFace = NULL;
    Edge *nextEdge = NULL;

    Hexahedron *hex = Hexahedron::down_cast(thisCell);

    // Get the next oppsite face through which you will enter next ...
    nextFace = hex->getOppositeFace( thisFace );
    assert( nextFace ) ;

    if( thisEdge )  {
        nextEdge = hex->getDiagonalEdge( thisEdge );
        assert( nextEdge ) ;
    }

    // Check if the chord is closed.
    if( chord->hasFace(thisFace) ) {
        return;
    }

    // Check if the chord is intersecting  ...
    if( chord->hasCell( thisCell ) ) {
        chord->addIntersection( thisCell );
    }

    chord->addElements(thisFace, thisEdge, thisCell);

    JCellSequence neighs;
    nextFace->getRelations( neighs );
    Cell *nextCell = NULL;
    if( neighs.size() == 2 ) {
        if( neighs[0] == thisCell ) nextCell = neighs[1];
        if( neighs[1] == thisCell ) nextCell = neighs[0];
    }

    if( nextCell == NULL) {
        chord->addElements( nextFace, nextEdge );
        return;
    }

    expandChord( nextCell, nextFace, nextEdge, chord);
}

////////////////////////////////////////////////////////////////////////////////////


DualChord* HexDual :: getChord(Face *thisFace, Edge *thisEdge)
{

#ifdef CSV
    DualChord *chord = new DualChord;
    chord->initialize();

    if( thisFace == NULL ) return NULL;

    if( pmesh->getAdjTable(2,3) == 0)
        pmesh->buildRelations(2,3);

    JCellSequence cellneighs;
    thisFace->getRelations( cellneighs );

    if( cellneighs.empty() ) return NULL;
    if( cellneighs.size() > 2 ) return NULL;

    Face *startFace = thisFace;
    Edge *startEdge = thisEdge;
    expandChord( cellneighs[0], startFace, startEdge, chord);

    if( cellneighs.size() == 2 && chord->isBoundary() ) {
        chord->reverse();
        startFace = thisFace;
        startEdge = thisEdge;
        expandChord( cellneighs[1], startFace, startEdge, chord);
    }

    chord->finalize();
#endif

    return chord;
}

////////////////////////////////////////////////////////////////////////////////////

int HexDual :: getChords( JFaceSequence &faces, vector<DualChord*> &chords)
{
    chords.clear();

    size_t nSize = faces.size();
    if( nSize == 0) return 1;

    chords.reserve( nSize );

    DualChord ch;
    for( size_t i = 0; i < nSize; i++) {
        DualChord *ch = getChord( faces[i]);
        if( ch )  chords.push_back(ch);
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////

int HexDual :: getAllChords( vector<DualChord*> &dchords)
{
    dchords.clear();

    pmesh->getLogger()->setInfo("Searching all chords in the hex mesh");

    Face *face;

    size_t numfaces = pmesh->getSize(2);
    int val = 0;
    for( size_t i = 0 ; i < numfaces; i++) {
        face = pmesh->getFaceAt(i);
        if( face->isActive() ) face->setAttribute("ChordID", val);
    }

    int currID = 1;

    JFaceSequence chfaces;
    while(1) {
        numfaces = pmesh->getSize(2);
        Face *seedFace = NULL;
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

        DualChord *ch = getChord( seedFace);

        ch->getFaces( chfaces);
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

///////////////////////////////////////////////////////////////////////////////

int HexDual :: getColumn(Face *face, JCellSequence &columnCells)
{
    columnCells.clear();

    DualChord *chord = getChord( face );

    chord->getCells(columnCells);
    chord->delete_dual_segments();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
/*
int HexDual :: getSheet( Face *face, vector<DualChord*> &chords)
{
     chords.clear();

     DualChord *commChord = getChord(face);
     chords.push_back(commChord);

     JFaceSequence commFaces;
     commChord->getFaces(commFaces);

     JCellSequence commCells;
     commChord->getCells(commCells);

     int nSize = commChord->getSize();

     DualChord *ch;
     for(int i = 0; i < nSize-1; i++) {
          Face *face = commFaces[i];
          Cell *cell = commCells[i];
          int pos = cell->getPosOf(face);
          if( pos == 0 || pos == 1 ) {
               ch = getChord( cell->getFaceAt(2));
               if( ch ) chords.push_back(ch);
               ch = getChord( cell->getFaceAt(4));
               if( ch ) chords.push_back(ch);
          }
          if( pos == 2 || pos == 3 ) {
               ch = getChord( cell->getFaceAt(0));
               if( ch ) chords.push_back(ch);
               ch = getChord( cell->getFaceAt(4));
               if( ch ) chords.push_back(ch);
          }

          if( pos == 4 || pos == 5 ) {
               ch = getChord( cell->getFaceAt(0));
               if( ch ) chords.push_back(ch);
               ch = getChord( cell->getFaceAt(2));
               if( ch ) chords.push_back(ch);
          }
     }
     return 0;
}
*/

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
            Vertex::getMidPoint(v0, v1, p0, 0.25);
            Vertex::getMidPoint(v0, v1, p1, 0.75);
            v0->setXYZCoords(p0);
            v1->setXYZCoords(p1);
        }

        // v0 is free node and v1 is bounded, Shift the free node only.
        if( b0 == 0 && b1 == 1) {
            Vertex::getMidPoint(v0, v1, p0, 0.25);
            v0->setXYZCoords(p0);
        }

        // v0 is bounded node and v0 is free, Shift the free node only.
        if( b0 == 1 && b1 == 0) {
            Vertex::getMidPoint(v1, v0, p0, 0.25);
            v1->setXYZCoords(p0);
        }

        // Both are in the boundary, dangerous case....
        if( b0 == 1 && b1 == 1) {
            bool special = 0;
            if( v0->hasAttribute("Steiner") && v1->hasAttribute("Steiner") ) {
                v0->getAttribute("Steiner", s0);
                v1->getAttribute("Steiner", s1);
                if( s0 == 1 && s1 == 2) {
                    Vertex::getMidPoint(v0, v1, p0, 0.65);
                    v1->setXYZCoords(p0);
                    special = 1;
                }
                if( s0 == 2 && s1 == 1) {
                    Vertex::getMidPoint(v1, v0, p0, 0.65);
                    v0->setXYZCoords(p0);
                    special = 1;
                }
            }

            if( v0->hasAttribute("Steiner") && !v1->hasAttribute("Steiner") ) {
                Vertex::getMidPoint(v0, v1, p0, 0.25);
                v0->setXYZCoords(p0);
                special = 1;
            }

            if( !v0->hasAttribute("Steiner") && v1->hasAttribute("Steiner") ) {
                Vertex::getMidPoint(v1, v0, p0, 0.25);
                v1->setXYZCoords(p0);
                special = 1;
            }
            if( !special ) {
                Vertex::getMidPoint(v0, v1, p0, 0.50);
                v0->setXYZCoords(p0);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
int HexDual :: removeSheet(DualSheet *sheet)
{
#ifdef CSV
    size_t numedges, numfaces, numcells;

    assert( pmesh->getAdjTable(0,1) );
    assert( pmesh->getAdjTable(0,2) );
    assert( pmesh->getAdjTable(0,3) );

    JCellSequence cells, cellneighs;
    JFaceSequence faceneighs;
    JEdgeSequence edgeneighs;

    sheet->getCells(cells);
    numcells = cells.size();
    for(size_t i = 0; i < numcells; i++) {
        Cell *cell = cells[i];
        assert( cell->isActive() );
        cell->setStatus( MeshEntity::REMOVE);
    }

    pmesh->buildRelations(0,3);

    JEdgeSequence paredges;
    sheet->getEdges(paredges);
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
        if( b0 == 0  && b1 == 0 ) vm = Vertex::getMidNode(v0, v1);

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
                vm = Vertex::getMidNode(v0, v1);
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
    sheet->getSheetFaces( faces );
    for(size_t i = 0; i < faces.size(); i++)
        faces[i]->setStatus( MeshEntity::REMOVE );
    sheet->clear();

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
#endif

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


    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void HexDual :: expandSheet( Cell *thisCell, Face *thisFace, Edge *thisEdge, DualSheet *sheet)
{
    if( thisFace == NULL || thisEdge == NULL ) return;

    if( sheet->hasFace(thisFace) ) {
        sheet->ring =  1;
        return;
    }

    sheet->addFace(thisFace);
    if( thisCell == NULL  || !thisCell->isActive() ) return;

    Hexahedron *hex = Hexahedron::down_cast(thisCell);
    hex->get_topological_parallel_edges(thisEdge, par_edges);

    if( sheet->hasCell( thisCell ) ) {
        sheet->self_intersecting  =  1;
    } else {
        sheet->addElements(thisCell, thisEdge, par_edges[0], par_edges[1], par_edges[2] );
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

int HexDual :: getSheets(Face *f, vector<DualSheet*> &sheets)
{
    sheets.clear();
    if( !f->isActive() ) return 1;

    /*
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
    */
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int HexDual :: getSheets(Cell *s, vector<DualSheet*> &sheets)
{
    sheets.clear();
    vector<DualChord*> cellChords(3);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void HexDual :: expandSheet( Cell *thisCell, Edge *thisEdge, DualSheet* sheet)
{
    if( thisCell == NULL || thisEdge == NULL ) return;

    DualSheet::Key newKey(thisCell,thisEdge);

    if( sheet->hasKey(newKey) ) {
        sheet->ring =  1;
        return;
    }

    sheet->addKey(newKey);

    JEdgeSequence paredges;
    Hexahedron *hex = Hexahedron::down_cast(thisCell);
    hex->get_topological_parallel_edges(thisEdge, paredges);

    if( sheet->hasCell( thisCell ) )
        sheet->self_intersecting  =  1;

    sheet->addElements(thisCell, thisEdge, paredges[0], paredges[1], paredges[2] );

    JCellSequence cellneighs;
    for(int j = 0; j < 3; j++) {
        paredges[j]->getRelations(cellneighs);
        for(size_t i = 0; i < cellneighs.size(); i++) {
            Cell *nextCell = cellneighs[i];
            if( nextCell != thisCell)
                expandSheet(nextCell, paredges[j], sheet);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

DualSheet* HexDual :: getSheet(Edge *edge)
{
    if( pmesh == NULL ) return NULL;
    if( edge  == NULL ) return NULL;

    if( !edge->isActive() ) return NULL;

    if( pmesh->getAdjTable(1,3) == 0 )
        pmesh->buildRelations(1,3);

    JCellSequence cellneighs;
    edge->getRelations(cellneighs);

    if( cellneighs.empty() ) return NULL;

    DualSheet *sheet = new DualSheet;
    expandSheet( cellneighs[0], edge, sheet);

    return sheet;
}

////////////////////////////////////////////////////////////////////////////////

void HexDual :: getSheets(JEdgeSequence &eseq, vector<DualSheet*> &vsheets)
{
    vsheets.clear();
    for( size_t i = 0; i < eseq.size(); i++) {
        DualSheet *newsheet = getSheet( eseq[i] );
        if( newsheet ) vsheets.push_back(newsheet);
    }
}

////////////////////////////////////////////////////////////////////////////////

int HexDual :: getAllSheets( vector<DualSheet*> &sheets)
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

    vector<DualSheet*> sh;
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
            sh[j]->getFaces( chfaces);
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
                         vector<DualSheet*> &sheet1, vector<DualSheet*> &sheet2)
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
int HexDual :: removeColumn(DualChord *chord)
{
    /*
         if( !chord->isRemovable() ) return 1;

         JCellSequence cells;
         chord->getCells( cells );

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

         for( size_t i = 0; i < chord->collapseNodesPair.size()/2; i++) {
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
    */
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

    int nSize = chord1->getSize();

    JFaceSequence faceSeq1, faceSeq2;
    JCellSequence cellSeq1, cellSeq2;

    chord1->getFaces( faceSeq1 );
    chord2->getFaces( faceSeq2 );

    chord1->getCells( cellSeq1 );
    chord2->getCells( cellSeq2 );

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
    if(chord1->getSize() != chord2->getSize()  ) return 0;

    if( !chord1->isSimple() ) return 0;
    if( !chord2->isSimple() ) return 0;

    int nSize = chord1->getSize();

    JCellSequence cellSeq1, cellSeq2;
    chord1->getCells( cellSeq1 );
    chord2->getCells( cellSeq2 );

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
#endif



