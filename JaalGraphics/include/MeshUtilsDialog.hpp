#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshUtilsDialog.hpp"
#include "MeshViewer.hpp"
#include "QuadParamLines.hpp"


class JMeshUtilsDialog : public QDialog, public Ui::MeshUtilsDialog {
    Q_OBJECT

public:
    JMeshUtilsDialog( QWidget *parent = 0);
    ~JMeshUtilsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void init();
    void getQuadParamLines();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    boost::shared_ptr<JEdgeColor> paramLinesColor;

    void makeConnections();
};
