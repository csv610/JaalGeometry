#include "KnotsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JKnotsDialog :: JKnotsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    oldmesh = nullptr;
    newmesh = nullptr;
    newCurve   = new JCurve;
    meshViewer = nullptr;
    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JKnotsDialog :: ~JKnotsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JKnotsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}
///////////////////////////////////////////////////////////////////////////////

void JKnotsDialog :: readFromFile()
{
    /*
         QString qstr = QFileDialog::getOpenFileName(this,
                        *new QString("Select Mesh File "),
                        lastSelectedDirectory,
                        *new QString( "Mesh Format (*cub *.dat *.ele *msh *.node *.off *.obj *.smf *.xml *.vtk)"));

         string filename = qstr.toUtf8().constData();
         bool val = 1;
         JNodeSequence nodes;
         JEdgeSequence edges;
         JCurve *knotcurve;

         if (!filename.empty()) {
              knotcurve = JCurve::readFromFile(filename);
              knotcurve->setType(JCurve::CLOSE_CURVE);
              knotcurve->finalize();
              knotcurve->getNodes(nodes);
              knotcurve->getEdges(edges);
              meshViewer->attach(newCurve);
              for( size_t i = 0; i < nodes.size(); i++)
                   nodes[i]->setAttribute("Display", val);
              for( size_t i = 0; i < edges.size(); i++)
                   edges[i]->setAttribute("Display", val);
              newCurve->copyOf(knotcurve);
         }
    */
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JKnotsDialog :: closeDialog()
{
    this->parentWidget()->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JKnotsDialog :: makeConnections()
{
    connect( readFilePushButton,  SIGNAL( clicked() ), this, SLOT( readFromFile() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}

///////////////////////////////////////////////////////////////////////////////
