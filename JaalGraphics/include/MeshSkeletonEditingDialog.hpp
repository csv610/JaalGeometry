#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshSkeletonEditingDialog.hpp"

#include "MeshViewer.hpp"
#include "ShapeViewer.hpp"

#include "MeshSkeleton.hpp"
#include "MeshSlicer.hpp"
#include "MeshSkeletonContours.hpp"

class JMeshSkeletonEditingDialog : public QDialog, public Ui::MeshSkeletonEditingDialog {
    Q_OBJECT

public:
    JMeshSkeletonEditingDialog( QWidget *parent = 0);
    ~JMeshSkeletonEditingDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setSkeleton( const JMeshSkeletonPtr &s);

private slots:
    void selectBranch();
    void displayBranch();
    void getContours();
    void deleteContour();
    void removeBranch();
    void splitBranch();
    void reparameterize();
    void enumNodes();
    void getCap();
    void getPoles();
    void closeDialog();
    
private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JShapeViewerPtr shapesViewer;

    JMeshPtr surfMesh;
    JMeshPtr skelGraph;
    JMeshSkeletonPtr meshSkel;
    JMeshPtr selectedBranchMesh;
    JMeshPtr otherBranchesMesh;

    JSlicePtr currSlice;
    boost::scoped_ptr<JMeshSkeletonContours> skelContourPtr;

    void setMesh( const JMeshPtr &m);
    int  getClosestContour(const vector<JEdgeSequence> &c, const Point3D &p);

    void displayContours( const JMeshPtr &m);
    void getAllBranchContours( int branchID);
    void getOneBranchContour( int branchID);
    void getTwoBranchContours( int branchID);

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
