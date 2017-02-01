#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_PolyCubesDialog.hpp"
#include "TetMesherDialog.hpp"
#include "MeshLaplaceSmoothingDialog.hpp"
#include "MeshPaintingDialog.hpp"

#include "MeshViewer.hpp"
#include "PolyCubes.hpp"
#include "MeshEntityPicker.hpp"
#include "MeshPartitioner.hpp"

class JPolyCubesDialog : public QDialog, public Ui::PolyCubesDialog {
    Q_OBJECT

    JaalViewer *viewManager;
public:
    JPolyCubesDialog( QWidget *parent = 0);
    ~JPolyCubesDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

    void displaySide( int side, bool val);

private slots:

    void initialSegments();
    void checkSide();
    void checkTopology();
    void optSurfmesh();
    void optVolmesh();
    void realign();
    void showPatch();
    void regionGrow();
    void alignAlongXYZ();
    void openTetmesherDialog();
    void openPaintingDialog();
    void smoothInterface();
    void generate();
    void setCubeSide();
//   void meanShift();
    void loadPolyCubes();
    void showMesh();
    void integerSnap();

    void saveAs();

private:
    JMeshViewerPtr meshViewer;
    JPolyCubes   polyCubes;
    JMeshPtr     modelMesh, polyMesh1, polyMesh2, cubicalMesh, hexMesh;
    JMeshEntityPickerPtr entityPicker;
    JCellSequence  cells;
    int   countfaces[6];

    void straighten(const JEdgeSequence &seq);
    void flatten(const JMeshPtr &patch);

    boost::scoped_ptr<JMeshLaplaceSmoothingDialog> laplaceDialog;
    boost::scoped_ptr<JTetMesherDialog> tetmesherDialog;
    boost::scoped_ptr<JMeshPaintingDialog> meshpaintingDialog;

    void setMesh( const string &s);
    void init();
    void makeConnections();

    void assignColor( const JFacePtr &f);
    void assignColor( const JEdgePtr &e);
    void assignColor( const JNodePtr &v);
    void assignColors( const JMeshPtr &m);
};
