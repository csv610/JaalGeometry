#pragma once

#include <QDialog>

#include "Ui_MeshOptBoundaryLayerDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshOptBoundaryLayer.hpp"


class JContourNormalsViewer : public JViewComponent
{
public:
    void setMesh( const JMeshPtr &m);
    void draw();

    float scale = 0.01;
    bool  displayEdgeNormals = 0;
    bool  displayNodeNormals = 1;
private:
    JMeshPtr mesh;
    vector<float> eNormalsHead, eNormalsTail;
    vector<float> vNormalsHead, vNormalsTail;
};


class JMeshOptBoundaryLayerDialog : public QDialog, public Ui::MeshOptBoundaryLayerDialog {
    Q_OBJECT

public:
    JMeshOptBoundaryLayerDialog( QWidget *parent = 0);
    ~JMeshOptBoundaryLayerDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void displayNormals();
    void optLayer();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;
    boost::shared_ptr<JContourNormalsViewer>  contourNormalsViewer;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
