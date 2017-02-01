#include "EditMeshInterfaceDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JEditMeshInterfaceDialog :: JEditMeshInterfaceDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;

    model = new QStandardItemModel();
    tableView->setModel( model );
}

///////////////////////////////////////////////////////////////////////////////

JEditMeshInterfaceDialog :: ~JEditMeshInterfaceDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JEditMeshInterfaceDialog :: init()
{
    /*
        if( viewManager == nullptr ) return;
        JViewComponent *c = viewManager->getComponent("MeshViewer");
        meshViewer = dynamic_cast<JMeshViewer*>(c);
        if( meshViewer == nullptr ) return;
        setMesh( meshViewer->getCurrentMesh() );
    */
}

///////////////////////////////////////////////////////////////////////////////

void JEditMeshInterfaceDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////

void JEditMeshInterfaceDialog :: filltable()
{
    model->clear();

    model->setHorizontalHeaderItem(0, new QStandardItem(QString("ID")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("#Segments")));
    model->setHorizontalHeaderItem(2, new QStandardItem(QString("Length")));
    model->setHorizontalHeaderItem(3, new QStandardItem(QString("Dilation")));

    int numInterfaces = mp.getNumInterfaces();
    if( numInterfaces == 0) return;

    JEdgeSequence edges;
    vector<Item> items( numInterfaces);

    /*
         for( int i = 0; i < numInterfaces; i++) {
              mp.getInterface(i, edges);
              items[i].id = i;
              items[i].numSegments = edges.size();
              items[i].length = JEdgeGeometry::getLength(edges);
         }


         for( it = degreeCount.begin(); it != degreeCount.end(); ++it) {
              int ndegree = it->first;
              if( ndegree ) {
                  QList<QStandardItem*> newRow;
                  QStandardItem* item1 = new QStandardItem(QString("%0").arg(ndegree));
                  newRow.append(item1);

                  int nsize = it->second;
                  QStandardItem* item2 = new QStandardItem(QString("%1").arg(nsize));
                  newRow.append(item2);

                  double perc = nsize*100/(double)numNodes;
                  QStandardItem* item3 = new QStandardItem(QString("%2").arg(perc));
                  newRow.append(item3);
                  model->appendRow(newRow);
              }
          }
    */
}

///////////////////////////////////////////////////////////////////////////////

void JEditMeshInterfaceDialog :: makeConnections()
{
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
