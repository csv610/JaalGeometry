#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_CurveGenDialog.hpp"
#include "Curve.hpp"
#include "JaalViewer.hpp"
#include "MeshViewer.hpp"
#include "KnotsDialog.hpp"

class JCurveGenDialog : public QDialog, public Ui::CurveGenDialog
{
    Q_OBJECT

public:
    JCurveGenDialog( QWidget *parent = 0);
    ~JCurveGenDialog();

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }

private:
    void mouseMoveEvent( QMouseEvent *e);

private slots:

    void closeCurve();
    void deleteCurve();
    void deletePoint();
    void openKnotsDialog();
    void getKnots();
    void refineSegments();
    void startSketching();
    void createCanvas();
    void setCanvasColor();
    void closeDialog();

private:
    JMeshPtr mesh;
    JMeshPtr canvas;
    vector<JMeshPtr> meshdb;

    int    nCount;
    JCurve  *newCurve;

    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JNodePtr lastVertex;
    boost::scoped_ptr<JKnotsDialog> knotsDialog;

    void init();
    void makeConnections();
    void startCurve();
};
