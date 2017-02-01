#include "MainWindow.hpp"

JaalMainWindow :: JaalMainWindow( QWidget *parent) : QMainWindow(parent)
{
    setupUi(this);
    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

JaalMainWindow :: ~JaalMainWindow()
{
    clearup();
}
///////////////////////////////////////////////////////////////////////////////
void JaalMainWindow :: resizeEvent( QResizeEvent *e)
{
   
    int w = centralwidget->width();
    int h = centralwidget->height();
    viewer->resize(w,h);
}

void JaalMainWindow :: clearup()
{
}

void JaalMainWindow :: Quit()
{
//  meshViewer->clearAll();
    clearup();
    cout << " The Coredump is because of Log4Cxx library  " << endl;
    exit(0);
}
///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openNewDataDialog()
{
    if( loadNewDataDialog.get() == nullptr )
        loadNewDataDialog.reset(new JLoadNewDataDialog(this));

    loadNewDataDialog->setViewManager( viewer );
    loadNewDataDialog->show();
}
void JaalMainWindow :: openGlobalSettingsDialog()
{
    if( globalSettingsDialog.get() == nullptr )
        globalSettingsDialog.reset(new JGlobalSettingsDialog());

    globalSettingsDialog->setViewManager(viewer);
    globalSettingsDialog->show();
}

///////////////////////////////////////////////////////////////////////////////
void JaalMainWindow :: openMeshToolsDialog()
{
    if( meshToolsDialog.get() == nullptr )
        meshToolsDialog.reset(new JMeshToolsDialog(this));

    meshToolsDialog->setViewManager( viewer );
    meshToolsDialog->show();
}

void JaalMainWindow :: openObjectsListDialog()
{
    if( objectsListDialog.get() == nullptr )
        objectsListDialog.reset(new JObjectsListDialog(this));

    objectsListDialog->setViewManager( viewer);
    objectsListDialog->setType(0);
    objectsListDialog->show();
}


void JaalMainWindow :: makeConnections()
{
    connect( actionOpen,  SIGNAL(triggered() ), this, SLOT( openNewDataDialog() ) );
    connect( actionQuit,  SIGNAL(triggered()), this, SLOT( Quit() ));
    connect( actionMeshTools,  SIGNAL(triggered() ), this, SLOT( openMeshToolsDialog() ) );
    connect( actionGlobalSettings,   SIGNAL(triggered() ), this, SLOT( openGlobalSettingsDialog() ) );
    connect( actionObjectViewer, SIGNAL(triggered()),  this, SLOT( openObjectsListDialog() ));
}

///////////////////////////////////////////////////////////////////////////////

/*
///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openSaveAnimationDialog()
{
    if( saveAnimationDialog.get() == nullptr )
        saveAnimationDialog.reset(new JSaveAnimationDialog());

    saveAnimationDialog->setViewManager( viewer );
    saveAnimationDialog->show();
}
///////////////////////////////////////////////////////////////////////////////
}
void JaalMainWindow :: openHeatConductionDialog()
{
    if( heatConductionDialog == nullptr )
        heatConductionDialog.reset(new JHeatConductionDialog(this));

    heatConductionDialog->setViewManager( viewer);
    heatConductionDialog->show();
}

///////////////////////////////////////////////////////////////////////////////


void JaalMainWindow :: openMeshTangleDialog()
{
    if( meshTangleDialog.get() == nullptr )
        meshTangleDialog.reset(new JMeshTangleDialog(this));

    meshTangleDialog->setViewManager( viewer);
    meshTangleDialog->show();

}

///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openTangleFEMTestsDialog()
{
    if( tanglefemtestsDialog == nullptr )
        tanglefemtestsDialog.reset(new JTangleFEMTestsDialog(this));

    tanglefemtestsDialog->setViewManager( viewer);
    tanglefemtestsDialog->show();
}

*/
