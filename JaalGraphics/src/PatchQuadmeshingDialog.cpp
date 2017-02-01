#include "PatchQuadmeshingDialog.hpp"
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

JPatchQuadmeshingDialog :: JPatchQuadmeshingDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

JPatchQuadmeshingDialog :: ~JPatchQuadmeshingDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JPatchQuadmeshingDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponent *c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_cast<JMeshViewer*>(c);

    if( meshViewer == nullptr ) return;

//   mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return ;
    qClean.reset(new JQuadCleanUp( mesh ));

    /*
         meshviewer->displayAll( 2, 1);
         meshviewer->getDrawFace()->setColorMethod( patchColor );

         mesh->buildRelations(0,2);
         meshviewer->getDrawNode()->setColorMethod( irregularColor );
    */

}

///////////////////////////////////////////////////////////////////////////////

void JPatchQuadmeshingDialog :: setColor()
{

}

///////////////////////////////////////////////////////////////////////////////

void JPatchQuadmeshingDialog :: remeshAll()
{
    if( mesh == nullptr ) return ;

//   mesh->buildRelations(0,2);

//  qClean->remesh_defective_patches();
//   meshViewer->setNewMesh(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JPatchQuadmeshingDialog :: accept()
{
//   meshViewer->getDrawFace()->setColorMethod( nullptr );
//   meshViewer->getDrawNode()->setColorMethod( nullptr );
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JPatchQuadmeshingDialog :: openTemplateDialog()
{
    if( meshingTemplatesDialog == nullptr ) {
        meshingTemplatesDialog.reset(new JQuadmeshingTemplatesDialog(this));
    }
    meshingTemplatesDialog->setViewManager( viewManager );
    meshingTemplatesDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JPatchQuadmeshingDialog :: makeConnections()
{
    PushButton( remeshPushButton, [=] {remeshAll();});
    PushButton( availableTemplatesPushButton, [=] {openTemplateDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
