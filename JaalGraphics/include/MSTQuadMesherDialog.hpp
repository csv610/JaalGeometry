#pragma once

#include <QDialog>
#include <fstream>

#include "Ui_MSTQuadMesherDialog.hpp"
#include "MeshViewer.hpp"

#include "MSTQuadMesher.hpp"
#include "MeshSingularityGraphDialog.hpp"
#include "MeshUntangle.hpp"
#include "MeshAffineTransforms.hpp"
#include "MeshPartitionColors.hpp"

class JMSTQuadMesherDialog : public QDialog, public Ui::MSTQuadMesherDialog {

    Q_OBJECT
public:
    JMSTQuadMesherDialog( QWidget *parent = 0);

    ~JMSTQuadMesherDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);
    void setPatchCenters( const JNodeSequence &n);

private slots:
    void openSingularityGraphDialog();

    void unittest();
    void remeshPatch();
    void randomSegments();
    void savePattern();

    void simplifyBoundary();
    void smoothBoundary();
    void smoothMesh();
    void getMedialAxis();

    void detectCorners();
    int  getDefectivePatch();
    void getDefectAt();
    void untangleMesh();
    void repeatSearch();
    void setSeedSelect();
    void displaySingularNodes();
    void reorderSingularities();

    void checkPoint();
    void checkPointRecovery();

    void clearAll();
    void closeDialog();

protected:
     void showEvent( QShowEvent *e);
    void  mouseReleaseEvent( QMouseEvent *e);
    void  keyPressEvent( QKeyEvent *e);

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    void init();
    JMeshPtr mesh, testmesh, quadmesh;
    JMeshEntityPickerPtr nodePicker;

    string         pattern;
    std::ofstream  mylogfile;
    JNodeSet       singularNodes;
    JMeshPtr       singularGraph;
    JFaceSet       pickedFaces;
    JNodeSequence  pickedNodes;
    JMeshGeodesics djkPath;
    JNodeSequence  specifiedPatchCenters;

    QDefectivePatchPtr  defectivePatch;

    JMSTQuadMesher mstMesher;
    boost::scoped_ptr<JMeshSingularityGraphDialog> meshSingularityGraphDialog;

    void errorMessage();
    void displayPatch();
    void displayMesh();
    void addPatch();

    void makeConnections();
};
