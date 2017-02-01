#include "MeshDualViewer.hpp"
#include "MeshViewer.hpp"

using namespace Jaal;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
void
JMeshDualViewer::init()
{
    parallel_edges = 0;
    display_graph  = 0;
    display_duals  = 0;
    curr_sheet_dir = 0;
    currCounter    = 0;

    display_chord_faces = 0;
    display_chord_cells = 0;

    display_elements           = 1;
    display_touching_nodes     = 0;
    display_touching_edges     = 0;
    display_touching_faces     = 0;
    display_intersecting_faces = 0;
    display_intersecting_cells = 0;

    chordStyle   = 1;
    chordBeads   = 0;
    primalMesh   = nullptr;

    dualMeshViewer.reset(new JMeshViewer(viewManager));
//   primalMeshViewer = (JMeshViewer*)viewManager->getComponent("MeshViewer");
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualViewer :: setChord(const JDualChordPtr &ch)
{
    if( ch == nullptr ) return;
    vector<JDualChordPtr>  singleChord(1);
    singleChord[0] = ch;
    setChords( singleChord );
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDualViewer :: setChords( const vector<JDualChordPtr> &vch)
{
    /*
        chords = vch;
        JFaceSequence faces;

        JFaceRenderPtr fAttrib;
        if( primalMesh) {
            size_t numfaces = primalMesh->getSize(2);
            for( size_t i = 0; i < numfaces; i++) {
                JFacePtr f = primalMesh->getFaceAt(i);
                f->getAttribute("Render", fAttrib);
                fAttrib->display = 0;
            }
        }

        for( size_t i = 0; i < vch.size(); i++) {
            faces = vch[i]->getMeshFaces();
            for( size_t j = 0; j < faces.size(); j++) {
                faces[j]->getAttribute("Render", fAttrib);
                fAttrib->display = 1;
            }

        }

        chords = vch;
        if( chordStyle ) {
            for( size_t i = 0; i < chordTubes.size(); i++)
                chordTubes[i].clear();
            chordTubes.resize( vch.size() );
            for( size_t i = 0; i < vch.size(); i++) {
                JEdgeSequence segments = vch[i]->getSegments();
                chordTubes[i].setSegments( segments );
            }
        }
        refreshDisplay();
    */
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshDualViewer::draw_sheet(const JDualSheetPtr &sheet, int id)
{
#ifdef CSV
    sheetColor.setSheetID(id);

    JFaceSequence sheetfaces;
    sheet->getSheetFaces( sheetfaces );

    size_t numfaces = sheetfaces.size();
    if( numfaces == 0) return;

    Point3D xyz;
    glLineWidth(2.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor4f( 0.1, 0.1, 0.1, 0.5);
    glBegin(GL_QUADS);
    for (size_t i = 0; i < numfaces; i++) {
        JFacePtr f =  sheetfaces[i];
        xyz = f->getNodeAt(0)->getXYZCoords();
        glVertex3f( xyz[0], xyz[1], xyz[2] );

        xyz = f->getNodeAt(1)->getXYZCoords();
        glVertex3f( xyz[0], xyz[1], xyz[2] );

        xyz = f->getNodeAt(2)->getXYZCoords();
        glVertex3f( xyz[0], xyz[1], xyz[2] );

        xyz = f->getNodeAt(3)->getXYZCoords();
        glVertex3f( xyz[0], xyz[1], xyz[2] );
    }
    glEnd();
    glDisable( GL_BLEND );
    glLineWidth(1.0);

    Color color = sheetColor.getColor(id);
    glColor4f( color[0], color[1], color[2], 0.5 );

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glBegin(GL_QUADS);
    for (size_t i = 0; i < numfaces; i++) {
        JFacePtr f =  sheetfaces[i];
        xyz = f->getNodeAt(0)->getXYZCoords();
        glVertex3f( xyz[0], xyz[1], xyz[2] );

        xyz = f->getNodeAt(1)->getXYZCoords();
        glVertex3f( xyz[0], xyz[1], xyz[2] );

        xyz = f->getNodeAt(2)->getXYZCoords();
        glVertex3f( xyz[0], xyz[1], xyz[2] );

        xyz = f->getNodeAt(3)->getXYZCoords();
        glVertex3f( xyz[0], xyz[1], xyz[2] );
    }
    glEnd();


    if( display_elements ) {
        JCellSequence sheetcells;
        sheet->getCells(sheetcells);
        size_t numcells = sheetcells.size();

        for (size_t i = 0; i < numcells; i++) {
            JCellPtr c = sheetcells[i];
            activate(c);
        }
    }
#endif
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualViewer :: activate( const JFacePtr &face )
{
    JFaceRenderPtr fAttrib;
    face->getAttribute("Render", fAttrib);
    fAttrib->display = 1;
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDualViewer :: activate( const JCellPtr &cell )
{
    bool val = 1;
    cell->setAttribute("Display", val);
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshDualViewer::draw_chord( const JDualChordPtr &chord, int id)
{
    if( chord == nullptr ) return;

#ifdef CSV
    ChordNodeColor  chordNodeColor;
    ChordEdgeColor  chordEdgeColor;

    chordNodeColor.setChordID(id);
    chordEdgeColor.setChordID(id);

    size_t nSize;
    if( chordBeads ) {
        chord->getNodes( chordNodes );
        nSize = chordNodes.size();
        bool val = 1;
        for (size_t i = 0; i < nSize; i++) {
            chordNodeColor.assign(chordNodes[i]);
        }
    }
    /*
              if( display_duals ) {
                   DrawEdge *drawEdge = meshViewer->getDrawEdge();
                   if( drawEdge ) {
                        chordTubes[id].setRadius( 1.01*drawEdge->getCylinderRadius() );
                        chordTubes[id].setNumSides(drawEdge->getNumCylinderSides() );
                        chordTubes[id].draw();

                        if( parallel_edges) {
                             drawEdge->setGlyph(2);
                             chord->getEdges( parEdges );
                             drawEdge->setColorMethod( nullptr );
                             for( size_t i = 0; i < parEdges.size(); i++)
                                  drawEdge->draw(parEdges[i]);
                        }

                   }
              }
    */

    if( display_elements ) {
        glDisable( GL_LIGHTING );
        glDisable( GL_CULL_FACE );
        chordFaceColor.setChordID( id );
        chordFaces = chord->getMeshFaces();
        nSize = chordFaces.size();
        bool val = 1;
        for (size_t i = 0; i < nSize; i++) {
            JFacePtr f = chordFaces[i];
            f->setAttribute("ChordID", id );
            chordFaceColor.assign(f);
        }

        /*
                  chord->getCells(chordCells);
                  nSize = chordCells.size();
                  for (size_t i = 0; i < nSize; i++) {
                       Cell *c = chordCells[i];
                       activate(c);
                  }
        */
    }
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////

void JMeshDualViewer :: deleteAll()
{
    /*
         size_t nSize;

         nSize = chordTubes.size();
         for( size_t i = 0; i < nSize; i++)
              chordTubes[i].clear();
         chordTubes.clear();

         nSize = chords.size();
         for( size_t i = 0; i < nSize; i++)
              chords[i]->Delete();
         chords.clear();

         nSize = sheets.size();
         for( size_t i = 0; i < nSize; i++)
              sheets[i]->Delete();
         sheets.clear();

         mesh = meshViewer->getMesh(); // Since mesh object can change in MeshViewer.
         if( mesh ) {
              mesh->delete_edge_attribute("Steiner");
              mesh->delete_face_attribute("Steiner");
              mesh->delete_cell_attribute("Steiner");
         }
    */
}
///////////////////////////////////////////////////////////////////////////////////////////
void
JMeshDualViewer::draw()
{
//  if( dualMeshViewer == nullptr ) return;

//    dualMeshViewer->getDrawEdge()->setSteiner(1);

//    dualMeshViewer->draw();

    for( size_t i = 0; i < chords.size(); i++)
        draw_chord( chords[i], i);

    for( size_t i = 0; i < sheets.size(); i++)
        draw_sheet( sheets[i], i);

}


////////////////////////////////////////////////////////////////////////////////
