#pragma once

#include <QDialog>
#include "MeshViewer.hpp"
#include "PolyCubesDialog.hpp"
#include "LegoBuilderDialog.hpp"
#include "MeshStackDialog.hpp"
#include "MeshVoxelizerDialog.hpp"

#include "Ui_HexMesherDialog.hpp"
#include "BernHexOpsDialog.hpp"
#include "MeshOctree.hpp"
#include "SphereHexMesherDialog.hpp"

//#include "AllHexMeshGenerator.hpp"

class JHexMesherDialog : public QDialog, public Ui::HexMesherDialog {
    Q_OBJECT

    /*
        struct MyThread : public QThread {
            void run();
            string cmd;
            JMeshPtr mesh;
            JTetMesherDialog *dialog;
        };
        friend class MyThread;
        boost::scoped_ptr<MyThread>  thread;
    */

public:
    JHexMesherDialog( QWidget *parent = 0);
    ~JHexMesherDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

protected:
    virtual void  keyPressEvent( QKeyEvent *e);

private slots:
//  void genNewMesh();
    void closeDialog();
//    void openOptionDialog();
    void getGeodeTemplate();
    void openPolyCubesDialog();
    void openLegoMeshDialog();
    void openBernHexOpsDialog();
    void openMeshStackDialog();
    void allTet2HexConversion();
    void openStructuredMeshDialog();
    void getSchneiderTemplate();
    void getOctreeMesh();
//  void openMeshVoxelizerDialog();
    void openSphHexMesherDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;
    boost::scoped_ptr<JBernHexOpsDialog>     bernHexOpsDialog;
    boost::scoped_ptr<JLegoBuilderDialog>    legoMeshDialog;
    boost::scoped_ptr<JPolyCubesDialog>      polycubesDialog;
    boost::scoped_ptr<JMeshStackDialog>      meshStackDialog;
    boost::scoped_ptr<JMeshVoxelizerDialog>  meshVoxelizerDialog;
    boost::scoped_ptr<JStructuredMeshDialog> structMeshDialog;
    boost::scoped_ptr<JSphereHexMesherDialog> sphHexMesherDialog;

    void init();
    void makeConnections();
};
