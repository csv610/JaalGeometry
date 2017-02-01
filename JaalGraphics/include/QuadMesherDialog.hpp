#pragma once

#include <QDialog>

#include "Ui_QuadMesherDialog.hpp"

#include "MeshViewer.hpp"

#include "Gmsh2DDialog.hpp"
#include "InstantMeshDialog.hpp"
#include "PolyPartitioner.hpp"
#include "QuadCheckerBoard.hpp"
#include "MSTQuadMesherDialog.hpp"
#include "QuadMeshCleanupDialog.hpp"
#include "QuadDominant2PureQuadsDialog.hpp"
#include "QuadVerdictDialog.hpp"
#include "MeshMinSingularity.hpp"
#include "AlphaMSTQuadMeshDialog.hpp"
#include "MeshSingularityGraphDialog.hpp"
#include "StructuredMeshDialog.hpp"
#include "RingQuadsDialog.hpp"
#include "MeshExtrudeDialog.hpp"
#include "CrossField2DDialog.hpp"

class JQuadMesherDialog : public QDialog, public Ui::QuadMesherDialog {
    Q_OBJECT

    JaalViewer *viewManager;
public:
    static const int  SIMPLE_ALGORITHM           = 0;
    static const int  NSIDED_TEMPLATES           = 1;
    static const int  TREE_MATCHING_ALGORITHM    = 2;
    static const int  GRAPH_MATCHING_ALGORITHM   = 3;
    static const int  BLOSSUM_MATCHING_ALGORITHM = 4;

    JQuadMesherDialog( QWidget *parent = 0);
    ~JQuadMesherDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    JMeshPtr getQuadMesh(int method);

private slots:
    void openStructMeshDialog();
    void openMSTQuadMesherDialog();
    void openAlphaMSTDialog();
    void openCleanupDialog();
    void openPureQuadsDialog();
    void openGmsh2DDialog();
    void openQualityDialog();
    void openPolygonSimplifyDialog();
    void openCrossFieldDialog();
    void openSingularityGraphDialog();
    void openRingQuadsDialog();
    void openMeshExtrudeDialog();
    void openInstantMeshDialog();

    void genmesh();
    void getCheckerBoardPattern();
    void displaySingularNodes();
    void getCyclicQuads();
    void getBaseQuadMesh();
    void refineAll();

    void closeDialog();

protected:
    void showEvent( QShowEvent *event);

private:
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    void init();
    void makeConnections();
    JNodeSequence  singularNodes;

    boost::scoped_ptr<JGmsh2DDialog>  gmsh2DDialog;
    boost::scoped_ptr<JQuadVerdictDialog>  quadVerdictDialog;
    boost::scoped_ptr<JCrossField2DDialog>  crossFieldDialog;
    boost::scoped_ptr<JQuadMeshCleanupDialog>  cleanupDialog;
    boost::scoped_ptr<JInstantMeshDialog>    instantMeshDialog;
    boost::scoped_ptr<JQuadDominant2PureQuadsDialog>  pureQuadsDialog;
    boost::scoped_ptr<JMSTQuadMesherDialog>  mstQuadMesherDialog;
    boost::scoped_ptr<JAlphaMSTQuadMeshDialog>  alphaMSTDialog;
    boost::scoped_ptr<JMeshSingularityGraphDialog>  meshSingularityGraphDialog;
    boost::scoped_ptr<JStructuredMeshDialog>       structmeshDialog;
    boost::scoped_ptr<JMeshExtrudeDialog>    meshExtrudeDialog;
    boost::scoped_ptr<JRingQuadsDialog>       ringQuadsDialog;
/*
    boost::scoped_ptr<JPolygonSimplifyDialog>  polySimplifyDialog;
    boost::scoped_ptr<JMixedIntegerQuadsDialog>  mixedIntegerQuadsDialog;
    boost::scoped_ptr<JQuadEditingDialog>  quadEditingDialog;
    boost::scoped_ptr<JSketchQuadsDialog>  sketchQuadsDialog;
*/
    void setMesh( const JMeshPtr &m);
    void displayNonQuads();
    void simpleTri2Quad();
    void treeMatching();
    void getEdmondMatching();
    void blossumMatching();
    void displayStructure();
    void greedyMatching();
    void getHamiltonQuads();
    void refineOneBoundEdge(const JEdgePtr &e);
};
