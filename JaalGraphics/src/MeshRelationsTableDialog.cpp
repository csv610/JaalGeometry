#include "MeshRelationsTableDialog.hpp"


///////////////////////////////////////////////////////////////////////////////

JMeshRelationsTableDialog :: JMeshRelationsTableDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    mesh = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshRelationsTableDialog :: ~JMeshRelationsTableDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRelationsTableDialog :: init()
{
    if( mesh == nullptr) return;

    bool val;

    val = mesh->getAdjTable(0,0);
    rel00CheckBox->setChecked( val );

    val = mesh->getAdjTable(0,1);
    rel01CheckBox->setChecked( val );

    val = mesh->getAdjTable(0,2);
    rel02CheckBox->setChecked( val );

    val = mesh->getAdjTable(0,3);
    rel03CheckBox->setChecked( val );

    val = mesh->getAdjTable(1,0);
    rel10CheckBox->setChecked( val );

    val = mesh->getAdjTable(1,1);
    rel11CheckBox->setChecked( val );

    val = mesh->getAdjTable(1,2);
    rel12CheckBox->setChecked( val );

    val = mesh->getAdjTable(1,3);
    rel13CheckBox->setChecked( val );

    val = mesh->getAdjTable(2,0);
    rel20CheckBox->setChecked( val );

    val = mesh->getAdjTable(2,1);
    rel21CheckBox->setChecked( val );

    val = mesh->getAdjTable(2,2);
    rel22CheckBox->setChecked( val );

    val = mesh->getAdjTable(2,3);
    rel23CheckBox->setChecked( val );

    val = mesh->getAdjTable(3,0);
    rel30CheckBox->setChecked( val );

    val = mesh->getAdjTable(3,1);
    rel31CheckBox->setChecked( val );

    val = mesh->getAdjTable(3,2);
    rel32CheckBox->setChecked( val );

    val = mesh->getAdjTable(3,3);
    rel33CheckBox->setChecked( val );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRelationsTableDialog :: makeConnections()
{
    PushButton( closePushButton, [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
