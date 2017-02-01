#pragma once

#include <QDialog>

#include "Ui_MeshRefine2DDialog.hpp"

#include "MeshViewer.hpp"
#include "MeshRefine.hpp"
#include "AllTriMeshGenerator.hpp"
#include "QuadMesherDialog.hpp"
#include "MeshSubdivisionDialog.hpp"

class JMeshRefine2DDialog : public QDialog, public Ui::MeshRefine2DDialog {
    Q_OBJECT

public:
    JMeshRefine2DDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

protected:
    void mouseReleaseEvent(QMouseEvent *e);

private slots:
    void refinetri13();
    void refinetri14();
    void refinetri16();
    void tri2quads();

    void refineQuad14();
    void refineQuad15();

    void quadtri2();
    void quadtri4();

    void insertPillows();
    void refineEdge();
    void flipEdge();
    void openTri2QuadsDialog();

    void quadBlocks();
    void setRefineSet();
    void openMeshSubdivisionDialog();
    void closeDialog();

private:
    JaalViewer *viewManager;

    JMeshViewerPtr meshViewer;
    boost::scoped_ptr<JQuadMesherDialog> quadmesherDialog;
    boost::scoped_ptr<JMeshSubdivisionDialog> meshSubdivisionDialog;

    JMeshPtr mesh;
    JMeshEntityPickerPtr entityPicker;
    JFaceSet   selectedFaces;

    int  refine_tri_type;
    int  refine_quad_type;
    void refinetrimesh();

    void init();
    void projectOnCircle();
    void makeConnections();
};

