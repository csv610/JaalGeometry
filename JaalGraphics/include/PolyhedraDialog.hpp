#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_PolyhedraDialog.hpp"
#include "MeshViewer.hpp"
#include "BoySurface.hpp"

class JPolyhedraDialog : public QDialog, public Ui::PolyhedraDialog
{
    Q_OBJECT

public:
    JPolyhedraDialog( QWidget *parent = 0);
    ~JPolyhedraDialog();

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }

private slots:
    void genRegularPolytopes();
    void genBoySurface();

private:
    JMeshPtr mesh, newmesh;
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    void init();
    void platonicSolid();
    void archimedSolid();
    void makeConnections();
};

