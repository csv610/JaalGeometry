#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshStackDialog.hpp"
#include "MeshViewer.hpp"
#include "CurveGenDialog.hpp"
#include "MeshRefine2DDialog.hpp"


class JMeshStackDialog : public QDialog, public Ui::MeshStackDialog
{
    Q_OBJECT

public:
    JMeshStackDialog( QWidget *parent = 0);
    ~JMeshStackDialog();

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }

private slots:

    void openCurveDialog();
    void refineQuads();
    void gen3DMesh();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr oldmesh, baseQuadmesh, hexmesh;
    vector<JMeshPtr> meshPlanes;
    double   rotAngle;

    JCurve      *trajCurve;
    boost::scoped_ptr<JCurveGenDialog> curveDialog;
    boost::scoped_ptr<JMeshRefine2DDialog> refineDialog;

    void init();
    void makeConnections();
};

