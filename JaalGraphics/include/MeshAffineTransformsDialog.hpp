#pragma once

#include <QDialog>

#include "MeshViewer.hpp"
#include "MeshAffineTransforms.hpp"

#include "Ui_AffineTransformDialog.hpp"

class JMeshAffineTransformsDialog : public QDialog, public Ui::AffineTransformDialog {
    Q_OBJECT

public:
    JMeshAffineTransformsDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void applyTransform();
    void normalize();
    void toCenter();
    void changeNodeCoords();
    void keyPressEvent( QKeyEvent *e);
    void setNewNodeID();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;

    boost::scoped_ptr<JMeshAffineTransform> affine;

    void geomInfo();
    void init();
    void makeConnections();
};

