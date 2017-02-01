#include "MeshEntityAttribListDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshEntityAttribListDialog :: JMeshEntityAttribListDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    mesh = nullptr;
    entity = 0;

    model = new QStringListModel( this );
    listView->setModel( model );
}

///////////////////////////////////////////////////////////////////////////////

JMeshEntityAttribListDialog :: ~JMeshEntityAttribListDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEntityAttribListDialog :: init()
{
    if( mesh == nullptr ) return ;
    switch( entity )
    {
    case 0:
        nodeRadioButton->setChecked(true);
        break;
    case 1:
        edgeRadioButton->setChecked(true);
        break;
    case 2:
        faceRadioButton->setChecked(true);
        break;
    case 3:
        cellRadioButton->setChecked(true);
        break;
    }

    vector<string> attribnames;
    mesh->getAttributeNames( attribnames, entity);

    QStringList list;
    for( size_t i = 0; i < attribnames.size(); i++)
        list << QString::fromStdString(attribnames[i]);
    model->setStringList(list);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEntityAttribListDialog :: makeConnections()
{
//    PushButton( closePushButton,  [=]{ close(); });
}

///////////////////////////////////////////////////////////////////////////////
