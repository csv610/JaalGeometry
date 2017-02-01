#pragma once

#include <QDialog>
#include "Ui_MeshSkeletonDialog.hpp"

#include "MeshViewer.hpp"
#include "ShapeViewer.hpp"
#include "MeshAffineTransforms.hpp"
#include "MeshRenderDialog.hpp"
#include "MeshPartitioner.hpp"
#include "NodeAttributesDialog.hpp"
#include "EdgeAttributesDialog.hpp"
#include "FaceAttributesDialog.hpp"
#include "MeshPartitionDialog.hpp"

/*
#include "MeshSkeleton.hpp"
#include "SimpleLaplaceSkeleton.hpp"
#include "MeshSkeletonContours.hpp"
#include "MeshSkeletonEditingDialog.hpp"
#include "MeshSkeletonShapesDialog.hpp"
*/


class JMeshSkeletonDialog : public QDialog, public Ui::MeshSkeletonDialog {
    Q_OBJECT

public:
    JMeshSkeletonDialog( QWidget *parent = 0);
    ~JMeshSkeletonDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);
    void setParentAfterClose( bool b) { showParentAfterClose = b; }

private slots:
    void contractGeometry();
    void getSkeleton();

    void openMeshRenderDialog();

    // Various type of nodes on skeleton
    void openNodesDialog();
    void openLeafNodesDialog();
    void openJunctionNodesDialog();

    void openSkelEditingDialog();
    void openSkelShapesDialog();

    void displaySkeleton();
    void openEdgesDialog();
    void createSurfaceParts();
    void displaySurfaceParts();
    void displayDisks();
    void displaySpheres();

    void loadSkeleton();
    void saveSkeleton();
    void rejectSkeleton();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JShapeViewerPtr shapeViewer;
    bool showParentAfterClose = 0;

    bool meanCurvature = 1;
    JMeshPtr inMesh;      // Input mesh can be anymesh
    JMeshPtr surfMesh;    // Derived Surface mesh from the input mesh.
    JMeshPtr skelGraph, refineSkelGraph;  // A skeleton is a graph containging nodes and edges...
    JMeshPtr contourDisks;

    boost::scoped_ptr<JMeshRenderDialog> meshRenderDialog;
    boost::scoped_ptr<JNodeAttributesDialog> nodeAttribDialog;
    boost::scoped_ptr<JEdgeAttributesDialog> edgeAttribDialog;

/*
    JMeshSkeletonPtr meshSkel;
    boost::scoped_ptr<JSimpleLaplacianSkeleton> lapSkel;
    boost::scoped_ptr<JMeshSkeletonShapesDialog>  skelShapesDialog;
    boost::scoped_ptr<JMeshSkeletonEditingDialog> skelEditingDialog;
*/

    void init();
    void initMesh();
    void makeConnections();
    void getSurfTriMesh();
};
////////////////////////////////////////////////////////////////////////////////
