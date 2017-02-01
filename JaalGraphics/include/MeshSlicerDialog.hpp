#pragma once

#include <QDialog>

#include "Ui_MeshSlicerDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshSlicer.hpp"
#include "MeshAffineTransforms.hpp"
#include "EdgeAttributesDialog.hpp"

class JMeshSlicerDialog : public QDialog, public Ui::MeshSlicerDialog {
    Q_OBJECT

    JaalViewer *viewManager;
public:
    JMeshSlicerDialog( QWidget *parent = 0);
    ~JMeshSlicerDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

protected:
    void keyPressEvent( QKeyEvent *e);

private slots:

    void setPlane();
    void sliceMesh();
    void setDisplay();
    void clearAll();
    void openEdgeAttribDialog();
    void closeDialog();

private:
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;

    JMeshPtr plane;
    JMeshPtr boundBox;
    JMeshPtr contourMesh;
    boost::scoped_ptr<JEdgeAttributesDialog>  edgeAttribDialog;

    Point3D  meshCenter;
    JMeshSlicer slicer;

    void init();
    void makeConnections();
};
