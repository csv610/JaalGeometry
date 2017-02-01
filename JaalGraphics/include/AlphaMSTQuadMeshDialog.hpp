#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_AlphaMSTQuadMeshDialog.hpp"
#include "MeshViewer.hpp"

#include "CircleCover.hpp"
#include "DelaunayMesh.hpp"
#include "AlphaMSTQuadMesh.hpp"
#include "MeshOptDialog.hpp"
#include "MedialAxis2DDialog.hpp"
#include "MSTQuadMesherDialog.hpp"
#include "BasicShapes.hpp"

class JAlphaMSTQuadMeshDialog : public QDialog, public Ui::AlphaMSTQuadMeshDialog {
    Q_OBJECT

public:
    JAlphaMSTQuadMeshDialog( QWidget *parent = 0);
    ~JAlphaMSTQuadMeshDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

private slots:
    void  mouseReleaseEvent( QMouseEvent *e);
    void  showEvent( QShowEvent *e);

    void enablePicking();
    void hidePatch();
    void remesh();
    void getNewCircle();
    void showAllCircles();
    void showSingularities();
    void showCorners();
    void clearAll();
    void resetColor();
    void sweepAll();

    void openMeshOptDialog();
    void openPolyMesherDialog();
    void openMedialAxisDialog();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshEntityPickerPtr entityPicker;
    vector<JColor>  faceColors;

    JMeshPtr mesh;
    JMeshPtr coverMesh;
    JNodeSequence  patchNodes;
    vector<JCircle> medialCircles;
    boost::scoped_ptr<JAlphaMSTQuadMesh>    alphaMST;
    boost::scoped_ptr<JMedialAxis2DDialog>  medialAxisDialog;
    boost::scoped_ptr<JMSTQuadMesher>       mstMesher;
    boost::scoped_ptr<JMeshOptDialog>       meshOptDialog;
    boost::scoped_ptr<JMSTQuadMesherDialog>  mstQuadMesherDialog;

    void init();
    void makeConnections();
    void getNewPatch( const Point3D &);
};
////////////////////////////////////////////////////////////////////////////////
