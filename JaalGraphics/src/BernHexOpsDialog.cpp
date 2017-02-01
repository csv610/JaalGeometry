#include "BernHexOpsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JBernHexOpsDialog :: JBernHexOpsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    meshViewer  = nullptr;
    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JBernHexOpsDialog :: ~JBernHexOpsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

bool JBernHexOpsDialog :: isHexmesh()
{
    if( mesh == nullptr ) return 0;

    int dim = mesh->getTopology()->getDimension();
    if( dim == 3 ) {
        int etype = mesh->getTopology()->getElementsType(3);
        if( etype != 8 ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Bern's Hex operations require all hex-mesh ");
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) return 0;
        }
    }
    return 1;
}
///////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

    /*
         mesh = meshViewer->getMesh();
         if( mesh == nullptr ) return ;

         if( !isHexmesh() ) {
              mesh = nullptr;
              QMessageBox msg;
              msg.setIcon(QMessageBox::Warning);
              msg.setText("At present singlet operations are only for all quad mesh ");
              msg.setStandardButtons( QMessageBox::Ok);
              int ret = msg.exec();
              if( ret == QMessageBox::Ok ) {
                   return;
              }
         }
         bernOp.setMesh(mesh);
    */
}

///////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: demoOp1(int type)
{

#ifdef LATER
    JMeshPtr mytest = JMesh::newObject();

    JNodeSequence nodes;
    JCellSequence cells;

    if( type == 1 ) {
        JBernHexOps::getCanonical17(nodes, cells);
        mytest->addObjects( nodes );
        mytest->addObjects( cells );
    } else {
        JBernHexOps::getCanonical17(nodes, cells);
        mytest->addObjects( nodes );
        // Apply operations and get new nodesa and cells.
        JHexahedronPtr hex = Hexahedron::down_cast( cells[0] );
        JBernHexOps::Op1_7( hex, nodes, cells);
        mytest->addObjects( nodes );
        mytest->addObjects( cells );
    }

    meshViewer->getFaceDraw()->setTransparencyMethod( JFaceDraw::SCREENDOOR);
    meshViewer->addObject(mytest);
#endif
}
///////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: demoOp2(int type)
{
    JMeshPtr mytest = JMesh::newObject();

    JNodeSequence nodes;
    JCellSequence cells;

    if( type == 1 ) {
        JBernHexOps::getCanonical26(nodes, cells);
        mytest->addObjects( nodes );
        mytest->addObjects( cells );
    } else {
        JBernHexOps::getCanonical26(nodes, cells);
        mytest->addObjects( nodes );
        // Apply operations and get new nodesa and cells.
        JHexahedronPtr hex1 = JHexahedron::down_cast( cells[0] );
        JHexahedronPtr hex2 = JHexahedron::down_cast( cells[1] );

//        JBernHexOps::Op2_6( hex1, hex2, nodes, cells);
        mytest->addObjects( nodes );
        mytest->addObjects( cells );
    }

//    meshViewer->getFaceDraw()->setTransparencyMethod( JFaceDraw::SCREENDOOR);
    meshViewer->addObject(mytest);
}

////////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: demoOp3( int type)
{
    JMeshPtr mytest = JMesh::newObject();

    JNodeSequence nodes;
    JCellSequence cells;
    // Create one input cell...
    JBernHexOps::getCanonical314_5(nodes, cells);
    mytest->addObjects( nodes );
    mytest->addObjects( cells );

    if( type == 2) {
        JHexahedronPtr hex1 = JHexahedron::down_cast( cells[0] );
        JHexahedronPtr hex2 = JHexahedron::down_cast( cells[1] );
        JHexahedronPtr hex3 = JHexahedron::down_cast( cells[2] );
//        JBernHexOps::Op314_5( hex1, hex2, hex3, nodes, cells);
        mytest->addObjects( nodes );
        mytest->addObjects( cells );
    }

    meshViewer->addObject(mytest);
}

////////////////////////////////////////////////////////////////////////////////
void JBernHexOpsDialog :: demoOp4(int type)
{
    JMeshPtr mytest = JMesh::newObject();

    JNodeSequence nodes;
    JCellSequence cells;
    // Create one input cell...
    JBernHexOps::getCanonical316_5(nodes, cells);
    if( nodes.size() != 16 && cells.size() != 3 ) {
        return;
    }

    if( type == 2 ) {
        JHexahedronPtr hex1 = JHexahedron::down_cast( cells[0] );
        JHexahedronPtr hex2 = JHexahedron::down_cast( cells[1] );
        JHexahedronPtr hex3 = JHexahedron::down_cast( cells[2] );
//        JBernHexOps::Op316_5( hex1, hex2, hex3, cells);
    }

    mytest->addObjects( nodes );
    mytest->addObjects( cells );

    meshViewer->addObject(mytest);
}

////////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: demoOp5(int type)
{
    JMeshPtr mytest = JMesh::newObject();

    JNodeSequence nodes;
    JCellSequence cells;
    // Create one input cell...
    JBernHexOps::getCanonical416_4(nodes, cells);
    mytest->addObjects( nodes );

    if( type == 2 ) {
        JHexahedronPtr hex1 = JHexahedron::down_cast( cells[0] );
        JHexahedronPtr hex2 = JHexahedron::down_cast( cells[1] );
        JHexahedronPtr hex3 = JHexahedron::down_cast( cells[2] );
        JHexahedronPtr hex4 = JHexahedron::down_cast( cells[3] );
//          JBernHexOps::Op416_4( hex1, hex2, hex3, hex4, cells);
    }

    mytest->addObjects( cells );
    meshViewer->addObject(mytest);
}

////////////////////////////////////////////////////////////////////////////////
void JBernHexOpsDialog :: demoOp6(int type)
{
    JMeshPtr mytest = JMesh::newObject();

    JNodeSequence nodes;
    JCellSequence cells;
    // Create one input cell...
    JBernHexOps::getCanonical415_4(nodes, cells);

    /*
       Hexahedron *hex1 = Hexahedron::down_cast( cells[0] );
       Hexahedron *hex2 = Hexahedron::down_cast( cells[1] );
       BernHexOps::Op2_6( hex1, hex2, nodes, cells);
    */

    mytest->addObjects( nodes );
    mytest->addObjects( cells );
    meshViewer->addObject(mytest);
}

////////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: applyOp1()
{
    JMeshPtr mytest = JMesh::newObject();

    /*
         JNodeSequence nodes;
         JCellSequence cells;
         // Create one input cell...
         BernHexOps::getCanonical17(nodes, cells);

         // Apply operations and get new nodesa and cells.
         Hexahedron *hex = Hexahedron::down_cast( cells[0] );
         BernHexOps::Op1_7( hex, nodes, cells);

         mytest->addObjects( nodes );
         mytest->addObjects( cells );

         meshViewer->setNewMesh(mytest);
    */
}
///////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: applyOp2()
{
    JMeshPtr mytest = JMesh::newObject();

    /*
         JNodeSequence nodes;
         JCellSequence cells;
         BernHexOps::getCanonical26(nodes, cells);

         mytest->addObjects( nodes );

         // Apply operations and get new nodesa and cells.
         Hexahedron *hex1 = Hexahedron::down_cast( cells[0] );
         Hexahedron *hex2 = Hexahedron::down_cast( cells[1] );

         BernHexOps::Op2_6( hex1, hex2, nodes, cells);

         mytest->addObjects( nodes );
         mytest->addObjects( cells );

         meshViewer->setNewMesh(mytest);
    */
}

////////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: applyOp3()
{
    JMeshPtr mytest = JMesh::newObject();
    /*

         JNodeSequence nodes;
         JCellSequence cells;
         // Create one input cell...
         BernHexOps::getCanonical314_5(nodes, cells);
         mytest->addObjects( nodes );
         mytest->addObjects( cells );

         Hexahedron *hex1 = Hexahedron::down_cast( cells[0] );
         Hexahedron *hex2 = Hexahedron::down_cast( cells[1] );
         Hexahedron *hex3 = Hexahedron::down_cast( cells[2] );
         BernHexOps::Op314_5( hex1, hex2, hex3, nodes, cells);
         mytest->addObjects( nodes );
         mytest->addObjects( cells );

         meshViewer->setNewMesh(mytest);
    */
}

////////////////////////////////////////////////////////////////////////////////
void JBernHexOpsDialog :: applyOp4()
{
    /*
         Mesh *mytest = Mesh::newObject();

         JNodeSequence nodes;
         JCellSequence cells;
         // Create one input cell...
         BernHexOps::getCanonical316_5(nodes, cells);
         if( nodes.size() != 16 && cells.size() != 3 ) {
              return;
         }

         Hexahedron *hex1 = Hexahedron::down_cast( cells[0] );
         Hexahedron *hex2 = Hexahedron::down_cast( cells[1] );
         Hexahedron *hex3 = Hexahedron::down_cast( cells[2] );
         BernHexOps::Op316_5( hex1, hex2, hex3, cells);

         mytest->addObjects( nodes );
         mytest->addObjects( cells );

         meshViewer->setNewMesh(mytest);
    */
}

////////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: applyOp5()
{
    /*
         Mesh *mytest = Mesh::newObject();

         JNodeSequence nodes;
         JCellSequence cells;
         // Create one input cell...
         BernHexOps::getCanonical416_4(nodes, cells);
         mytest->addObjects( nodes );

         Hexahedron *hex1 = Hexahedron::down_cast( cells[0] );
         Hexahedron *hex2 = Hexahedron::down_cast( cells[1] );
         Hexahedron *hex3 = Hexahedron::down_cast( cells[2] );
         Hexahedron *hex4 = Hexahedron::down_cast( cells[3] );

         BernHexOps::Op416_4( hex1, hex2, hex3, hex4, cells);
         mytest->addObjects( cells );

         meshViewer->setNewMesh(mytest);
    */
}

////////////////////////////////////////////////////////////////////////////////
void JBernHexOpsDialog :: applyOp6()
{
    /*
         Mesh *mytest = Mesh::newObject();

         JNodeSequence nodes;
         JCellSequence cells;
         // Create one input cell...
         BernHexOps::getCanonical415_4(nodes, cells);

         mytest->addObjects( nodes );
         mytest->addObjects( cells );

         meshViewer->setNewMesh(mytest);
    */
}

////////////////////////////////////////////////////////////////////////////////


void JBernHexOpsDialog :: applyOp()
{
    if( op1RadioButton->isChecked() ) applyOp1();
    if( op2RadioButton->isChecked() ) applyOp2();
    if( op3RadioButton->isChecked() ) applyOp3();
    if( op4RadioButton->isChecked() ) applyOp4();
    if( op5RadioButton->isChecked() ) applyOp5();
    if( op6RadioButton->isChecked() ) applyOp6();
}

////////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: cancelOp()
{
}

////////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: searchOp1()
{
    /*
         bernOp.searchPattern_7_1( localCells);

         if( localCells.empty() ) {
              QMessageBox msg;
              msg.setIcon(QMessageBox::Information);
              msg.setText("No suitable hex found for Bern's 7-1 operation");
              msg.setStandardButtons( QMessageBox::Ok);
              int ret = msg.exec();
              if( ret == QMessageBox::Ok ) return;
         }
    */
}

void JBernHexOpsDialog :: searchOp2()
{
}

void JBernHexOpsDialog :: searchOp3()
{
}
void JBernHexOpsDialog :: searchOp4()
{
}
void JBernHexOpsDialog :: searchOp5()
{
}
void JBernHexOpsDialog :: searchOp6()
{
}

////////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: demoOp(int type)
{
    if( op1RadioButton->isChecked() ) demoOp1(type);
    if( op2RadioButton->isChecked() ) demoOp2(type);
    if( op3RadioButton->isChecked() ) demoOp3(type);
    if( op4RadioButton->isChecked() ) demoOp4(type);
    if( op5RadioButton->isChecked() ) demoOp5(type);
    if( op6RadioButton->isChecked() ) demoOp6(type);
}
////////////////////////////////////////////////////////////////////////////////
void JBernHexOpsDialog :: demoElem1()
{
    demoOp(1);
}
void JBernHexOpsDialog :: demoElem2()
{
    demoOp(2);
}

void JBernHexOpsDialog :: patternSearch()
{
    if( op1RadioButton->isChecked() ) searchOp1();
    if( op2RadioButton->isChecked() ) searchOp2();
    if( op3RadioButton->isChecked() ) searchOp3();
    if( op4RadioButton->isChecked() ) searchOp4();
    if( op5RadioButton->isChecked() ) searchOp5();
    if( op6RadioButton->isChecked() ) searchOp6();
}
////////////////////////////////////////////////////////////////////////////////

void JBernHexOpsDialog :: makeConnections()
{
    connect( closePushButton,  SIGNAL( clicked() ), this, SLOT( close() ));
    connect( applyPushButton,  SIGNAL( clicked() ), this, SLOT( applyOp() ));
    connect( patternPushButton,  SIGNAL( clicked() ), this, SLOT(patternSearch() ));
    connect( demo1PushButton,  SIGNAL( clicked() ), this, SLOT( demoElem1() ));
    connect( demo2PushButton,  SIGNAL( clicked() ), this, SLOT( demoElem2() ));
    connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( cancelOp() ));
}

///////////////////////////////////////////////////////////////////////////////
