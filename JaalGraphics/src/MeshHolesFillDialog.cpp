#include "MeshHolesFillDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshHolesFillDialog :: JMeshHolesFillDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
    model = new QStandardItemModel();
    tableView->setModel( model );
}

///////////////////////////////////////////////////////////////////////////////

JMeshHolesFillDialog :: ~JMeshHolesFillDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshHolesFillDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshHolesFillDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////

void JMeshHolesFillDialog :: fillTable()
{
    model->clear();

    if( contours.empty() ) return;

    std::sort( contours.begin(), contours.end(), []( const JEdgeSequence  &a, const JEdgeSequence &b)
    {
        return a.size() < b.size();
    });

    model->setHorizontalHeaderItem(0, new QStandardItem(QString("ID")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("#Edges")));

    for( size_t i = 0; i < contours.size(); i++) {
        QList<QStandardItem*> newRow;
        QStandardItem* item1 = new QStandardItem(QString("%0").arg(i));
        newRow.append(item1);

        int nsize = contours[i].size();
        QStandardItem* item2 = new QStandardItem(QString("%1").arg(nsize));
        newRow.append(item2);

    }
}
///////////////////////////////////////////////////////////////////////////////

void JMeshHolesFillDialog :: searchHoles()
{
    if( mesh == nullptr) return;
    mesh->getTopology()->getBoundary(contours);

    if( contours.empty() ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Information);
        msg.setText("There are no boundary curves in the model");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
    fillTable();

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 1;
    mrender->displayEntity[2] = 0;
    mrender->displayEntity[3] = 0;
    size_t numEdges = mesh->getSize(1);
    JEdgeRenderPtr eAttrib;
    JColor red = JEntityColor::getColor("Red");

    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        edge->getAttribute("Render", eAttrib);
        if( edge->getNumRelations(2) == 1) {
            eAttrib->display   = 1;
            eAttrib->lineWidth = 3;
            eAttrib->color     = red;
        } else
            eAttrib->display   = 0;
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshHolesFillDialog ::  fillOne()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshHolesFillDialog ::  fillAll()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshHolesFillDialog ::  refineHole()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshHolesFillDialog ::  smoothHole()
{
}
///////////////////////////////////////////////////////////////////////////////

void JMeshHolesFillDialog ::  closeDialog()
{
}
///////////////////////////////////////////////////////////////////////////////

void JMeshHolesFillDialog :: makeConnections()
{
    PushButton( searchPushButton,    [=] {searchHoles(); });
    PushButton( fillOnePushButton, [=] {fillOne();});
    PushButton( fillAllPushButton, [=] {fillAll();});
    PushButton( closePushButton,     [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
