#include "EditPartitionRegionDialog.hpp"


JEditPartitionRegionDialog :: JEditPartitionRegionDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;

    model = new QStandardItemModel();
    tableView->setModel( model );
}

///////////////////////////////////////////////////////////////////////////////

JEditPartitionRegionDialog :: ~JEditPartitionRegionDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JEditPartitionRegionDialog :: init()
{
    /*
        if( viewManager == nullptr ) return;
        JViewComponentPtr c = viewManager->getComponent("MeshViewer");
        meshViewer = dynamic_cast<JMeshViewer*>(c);
        if( meshViewer == nullptr ) return;
        setMesh( meshViewer->getCurrentMesh() );
    */
}

///////////////////////////////////////////////////////////////////////////////

void JEditPartitionRegionDialog :: filltable()
{
    model->clear();

    model->setHorizontalHeaderItem(0, new QStandardItem(QString("ID")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("#Elements")));
    model->setHorizontalHeaderItem(2, new QStandardItem(QString("Area")));

    int numParts = mp.getNumPartitions();
    if( numParts == 0) return;

    /*
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
void JEditPartitionRegionDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str()));
    filltable();

    JMeshPartitioner mp;
    mp.setMesh(mesh);

}
///////////////////////////////////////////////////////////////////////////////

void JEditPartitionRegionDialog :: makeConnections()
{
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
