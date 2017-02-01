#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_SurfaceParameterizationDialog.hpp"

#include "MeshViewer.hpp"
#include "MeshPartitioner.hpp"
#include "SurfaceParameterization.hpp"
#include "MeshAffineTransforms.hpp"
#include "MeshPartitionDialog.hpp"
#include "MeshOptDialog.hpp"

////////////////////////////////////////////////////////////////////////////////

class JSurfaceMap2DViewer : public JViewComponent
{
public:
    JSurfaceMap2DViewer();
    ~JSurfaceMap2DViewer() {}

    void draw();

    void setPatch( const JMeshPtr &m);

    void setShowChart( bool b) {
        showChart = b;
    }

private:
    JMeshPtr patch;
    int  projCenter[2];
    int  projWidth, projHeight;
    std::map<int,JMeshPtr> mapPatches;

    bool showChart = 0;
    int  numRowPatches = 10;

    GLuint  texID;

    void init();
    void drawPatch();
    void getImage();
    void drawImage();
    void drawChart();
};

////////////////////////////////////////////////////////////////////////////////

class JSurfaceParameterizationDialog : public QDialog, public Ui::SurfaceParameterizationDialog {
    Q_OBJECT

public:
    JSurfaceParameterizationDialog( QWidget *parent = 0);
    ~JSurfaceParameterizationDialog() {}

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

private slots:
    void getUVCharts();
    void openMeshPartitionDialog();
    void smoothUVBoundary();
    void untangleUVMesh();
    void displayCharts();
    void checkDisplay();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    boost::shared_ptr<JSurfaceMap2DViewer >  surfMapViewer;

private:
    JMeshPtr mesh;

    vector<JMeshPtr> uvMeshes;

    boost::scoped_ptr<JMeshPartitioner>  meshPart;
    boost::scoped_ptr<JSurfaceParameterization> surfParam;
    boost::scoped_ptr<JMeshPartitionDialog>  meshPartDialog;

    void setParameters();

    void init();
    int  initMesh();
    void makeConnections();
};
