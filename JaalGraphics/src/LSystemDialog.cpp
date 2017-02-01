#include "LSystemDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

LSystemDialog :: LSystemDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    viewManager = nullptr;
    meshViewer  = nullptr;
    makeConnections();
    numIterationsLineEdit->setText( QString::number(1) );
}

///////////////////////////////////////////////////////////////////////////////

void LSystemDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void LSystemDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void LSystemDialog :: genShape()
{
    if( meshViewer ==  nullptr ) return;

    QString qstr = numIterationsLineEdit->text();
    int niter  = qstr.toInt();

    JMeshPtr mesh;

    if( sierpinskiCarpetRadioButton->isChecked() ) {
        mesh = AllQuadMeshGenerator::getSierpinski(niter);
    }

    if( sierpinskiTrianglesRadioButton->isChecked() )  {
        mesh = AllTriMeshGenerator::getSierpinski(niter);
    }

    if( sierpinskiPyramidRadioButton->isChecked() ) {
    }


    if( mengerSpongeRadioButton->isChecked() ) {
    }

    if( mesh ) meshViewer->addObject(mesh);

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void LSystemDialog :: applyDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void LSystemDialog :: cancelDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void LSystemDialog :: makeConnections()
{

    connect( applyPushButton, SIGNAL( clicked() ), this, SLOT( genShape() ));
    connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( cancelDialog() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
