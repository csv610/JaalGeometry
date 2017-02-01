#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshNormalsDialog.hpp"
#include "MeshViewer.hpp"

class JMeshNormalsDialog : public QDialog, public Ui::MeshNormalsDialog {
    Q_OBJECT

public:
    JMeshNormalsDialog( QWidget *parent = 0);
    ~JMeshNormalsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m, int pos);

    void setEntity( int e ) {
        normalAt = e;
    }

private slots:
    void setColor();
    void setScale();
    void setLength();
    void checkDisplay();
    void reverseAll();
    void getConsistent();
    void keyPressEvent( QKeyEvent *e);
    void recalculate();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    int   normalAt = 0;

    void init();
    void makeConnections();
};
