#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_StructuredMeshDialog.hpp"
#include "MeshViewer.hpp"

class JStructuredMeshDialog : public QDialog, public Ui::StructuredMeshDialog {
    Q_OBJECT

public:
    JStructuredMeshDialog( QWidget *parent = 0);
    ~JStructuredMeshDialog();

    void setViewManager(JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setOrigin( double x, double y, double z);
    void setLength( double x, double y, double z);
    void setDimension( int x, int y, int z);
    void setUVCoords( bool v = 0) {
        texCoords = v;
    }

signals:
    void meshCreated();

private slots:
    void genMesh();
    void keyPressEvent( QKeyEvent *e);
    void closeDialog();
    /*
         void newCubeIJK();
         void newCubeNodes();
    */

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr oldmesh, newmesh;
    bool  texCoords;

    void init();
    void makeConnections();
};
