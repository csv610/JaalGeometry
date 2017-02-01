#pragma once

#include <QDialog>
#include "Ui_MeshTangleDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshTangle.hpp"

class JMeshTangleViewer : public JViewComponent
{
public:

    void clear() {
        intersectPoints.clear();
    }
    void setIntersectPoints( vector<Point2D> &p) {
        intersectPoints = p;
    }
    void setNegativeFaces( JFaceSequence &f) {
        negativeFaces = f;
    }

    void draw();
private:
    vector<Point2D> intersectPoints;
    JFaceSequence    negativeFaces;
};

class JMeshTangleDialog : public QDialog, public Ui::MeshTangleDialog
{
    Q_OBJECT
    JaalViewer *viewManager;
public:
    JMeshTangleDialog( QWidget *parent = 0);
    ~JMeshTangleDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void getOrgMesh();
    void getTangle();
    void setColor();
    void moveNodes();

private:
    JMeshPtr mesh;
    JMeshViewerPtr meshViewer;
    boost::shared_ptr<JMeshTangleViewer> tangleViewer;

    JMeshTangle meshtangle;
    vector<double> orgCoords;
    vector<size_t> l2g;

    JFaceSequence tangledFaces;
    JFaceSequence negativeFaces;
    JEdgeSequence intersectEdges;

    void init();
    void makeConnections();
    void clearColors();
    void getNegativeFaces();
    void setNegativeColor();

    void setOverlapColor();
    void getOverlappedFaces();

    void setEdgesColor();
};
