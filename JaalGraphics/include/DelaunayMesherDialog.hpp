#pragma once

#include <QDialog>

#include "Ui_DelaunayMesherDialog.hpp"
#include "MeshViewer.hpp"
#include "DelaunayMesh.hpp"
#include "MeshGeometricQualityDialog.hpp"
#include "MeshTopologyQualityDialog.hpp"

/////////////////////////////////////////////////////////////////////////////////
class JDelaunayEdgeColor : public JEdgeColor {
public:
    string getName() const {
        return "DelaunayEdgeColor";
    }
    bool isPerEdge() const {
        return 1;
    }

    JDelaunayEdgeColor()
    {
        highlightColor[0] = 1.0;
        highlightColor[1] = 0.0;
        highlightColor[2] = 0.0;
        highlightColor[3] = 1.0;

        defaultColor[0]   = 0.6;
        defaultColor[1]   = 1.0;
        defaultColor[2]   = 0.6;
        defaultColor[3]   = 1.0;
    }

    int  assign( const JEdgePtr &e);
private:
    JColor  defaultColor;
};

/////////////////////////////////////////////////////////////////////////////////

class JDelaunayFaceColor : public JFaceColor {
public:
    string getName() const {
        return "DelaunayFaceColor";
    }

    bool isPerFace() const {
        return 1;
    }

    JDelaunayFaceColor();

    int  assign(const JFacePtr &);

    int  operator() (const JFacePtr &f) {
        return assign(f);
    }
};

/////////////////////////////////////////////////////////////////////////////////

class JDelaunayMesherDialog : public QDialog, public Ui::DelaunayMesherDialog {
    Q_OBJECT

    JaalViewer *viewManager;
public:
    JDelaunayMesherDialog( QWidget *parent = 0);
    ~JDelaunayMesherDialog();

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void generate();
    void openGeomQualityDialog();
    void openTopoQualityDialog();
    void closeDialog();

private:
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    JMeshPtr newMesh;

    boost::scoped_ptr<JDelaunayMesh2D> mesher;
    boost::scoped_ptr<JMeshTopologyQualityDialog>  topoQualityDialog;
    boost::scoped_ptr<JMeshGeometricQualityDialog> geomQualityDialog;

    void init();
    void makeConnections();
    void generateNewMesh();
};

/////////////////////////////////////////////////////////////////////////////////

