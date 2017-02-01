#include <assert.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <sstream>

#include "MeshQuadEditor.hpp"

using namespace Jaal;
using namespace std;

void
JMeshQuadEditorViewer::keyPressEvent(QKeyEvent *event)
{
#ifdef CSV
    if( event->modifiers() && Qt::ControlModifier) {
        display_sheet0 = 1;
        display_sheet1 = 1;
        display_sheet2 = 1;

        if (event->key() == Qt::Key_1) {
            display_sheet0 = 1;
            display_sheet1 = 0;
            display_sheet2 = 0;
            curr_sheet_dir = 0;
        }
        if (event->key() == Qt::Key_2) {
            display_sheet0 = 0;
            display_sheet1 = 1;
            display_sheet2 = 0;
            curr_sheet_dir = 1;
        }
        if (event->key() == Qt::Key_3) {
            display_sheet0 = 0;
            display_sheet1 = 0;
            display_sheet2 = 1;
            curr_sheet_dir = 2;
        }
        updateGL();
    }

    if (event->key() == Qt::Key_X) {
//        if( qdual ) qdual->dice( chords[0] );
//        if( hdual ) hdual->removeSheet( dsheets[curr_sheet_dir]  );
//        qdual->getChords(seedEdge, chords);
        /*
                  AllHexMeshGenerator::rotate_boundary_edge( mesh, seedEdge, parEdges);
                  if( hdual ) {
                       hdual->removeColumn( chords[0] );
                       JMeshViewer::init_mesh();
                       chords[0].clear();
                  }
                  seedFace = nullptr;
        */
        updateGL();
    }

    if (event->key() == Qt::Key_Z) {
        /*
         currCounter++;
        if( display_chords || display_sheets ) {
             dsheets[curr_sheet_dir].shrink();
        }
        */
//       chords[0].shrinkEdges();
        updateGL();
    }

    if (event->key() == Qt::Key_D) {
        currCounter = 0;
        if( display_chords) {
//              if( qdual ) qdual->getAllChords(chords);
            if( hdual ) hdual->getAllChords(chords);
            /*
                           chords.clear();
                           seedEdge = mesh->getRandomEdge();
                           chords.resize(1);
                           if( qdual ) qdual->getChord(seedEdge, chords[0] );
                           currChord = chords[0];
            */
            /*
                           chords.clear();
                           seedFace = mesh->getRandomFace();
                           chords.resize(1);
                           hdual->getChord(seedFace, chords[0] );
                           currChord = chords[0];
            */
        }

        /*
                 parEdges.clear();
                 seedEdge = mesh->getRandomEdge(1);
                 cout << "Debug : " << seedEdge->isBoundary() << endl;
                 cout <<  seedEdge->getNodeAt(0)->isBoundary() << " " << seedEdge->getNodeAt(1)->isBoundary() << endl;
        */
        /*
        static int ncounter = 0;
        if( ncounter == mesh->getSize(2)  ) ncounter = 0;
        for( int i = ncounter; i < mesh->getSize(2); i++) {
             seedFace = mesh->getFaceAt(i);
             if( seedFace->isBoundary() ) {
                 ncounter = i+1;
                 break;
             } else
                 seedFace = nullptr;
        }
        */

        //   seedFace = mesh->getRandomFace();

        /*
             if( display_chords) {
                  seedFace = mesh->getRandomFace(1);
                  chords.resize(1);
                  if( hdual ) hdual->getChord(seedFace, chords[0] );
             }
        */

        if( display_sheets ) {
            seedFace = mesh->getRandomFace(1);
            dsheets.clear();
            if( seedFace ) {
                JCellSequence cneigh;
                seedFace->getRelations(cneigh);
                Point3D p3d;
                seedFace->getAvgPos(p3d);
//                  setCenter(p3d);
//                  if( hdual ) hdual->getChord(seedFace, chords[0] );
                if( hdual ) hdual->getSheets(seedFace, dsheets);
                //      if( hdual ) hdual->getAllChords(chords);
            }
        }

        /*
        if( picked_faces.size()  )  {
             int fid  = *picked_faces.begin();
             seedFace = mesh->getFaceAt(fid);
             dual->getCurves(seedFace, dcurves);
        }
        */

        /*
        if( display_sheets) {
        seedCell = mesh->getRandomCell();
        Point3D p3d;
        seedCell->getAvgPos(p3d);
        setCenter(p3d);
        hdual->getSheets(seedCell, dsheets);
        }
        */

        /*
        if( picked_cells.size()  )  {
             int fid  = *picked_cells.begin();
             seedFace = mesh->getCellAt(fid);
             dual->getCurves(seedFace, dcurves);
        }
        */

        updateGL();
    }
#endif

    JMeshViewer::keyPressEvent(event);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshQuadEditorViewer :: animate()
{
    /*
         currCounter += 16;
         JMeshViewer::animate();
    */
}
///////////////////////////////////////////////////////////////////////////////

void
JMeshQuadEditorViewer::draw()
{
    bgColor[0] = background;
    bgColor[1] = background;
    bgColor[2] = background;
    glClearColor(bgColor[0], bgColor[1], bgColor[2], 0.0);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glCullFace( GL_BACK ) ;

    JMeshViewer::draw();
}

////////////////////////////////////////////////////////////////////////////////
void
JMeshQuadEditorViewer::init()
{
    JMeshViewer::init();
    assert(mesh);

    mesh->buildRelations(0,2);
}

////////////////////////////////////////////////////////////////////////////////

int
main(int argc, char** argv)
{
    assert(argc >= 2);
    QApplication application(argc, argv);

    JMeshQuadEditorViewer viewer;

    int iopt;
    while( (iopt = getopt( argc, argv, "d:i:") ) != -1) {
        switch(iopt) {
        case 'd':
            viewer.setDimension(atoi(optarg)) ;
            break;
            break;
        case 'i':
            viewer.addMeshFile(optarg);
            break;
        }
    }

    viewer.setWindowTitle("Jaal-QuadMeshEditor");
    NodeColor *fc = new IrregularNodeColor;
    viewer.setNodeColorMethod(fc);
//     FaceColor *fc = new PenroseTileColor;
//     EdgeColor *ec = new EdgeInterfaceColor;
//     viewer.setFaceColorMethod(fc);
//     viewer.setEdgeColorMethod(ec);

    // Make the viewer window visible on screen.
    viewer.show();

    // Run main loop.
    return application.exec();
}
////////////////////////////////////////////////////////////////////////////////
