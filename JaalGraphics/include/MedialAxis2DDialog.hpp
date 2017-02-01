#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MedialAxis2DDialog.hpp"
#include "MeshViewer.hpp"
#include "DelaunayMesh.hpp"
#include "NodeAttributesDialog.hpp"
#include "EdgeAttributesDialog.hpp"

class JMedialAxis2DDialog : public QDialog, public Ui::MedialAxis2DDialog {
    Q_OBJECT

public:
    JMedialAxis2DDialog( QWidget *parent = 0);
    ~JMedialAxis2DDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void getAxis();
    void closeDialog();
    void openNodeAttribsDialog();
    void openEdgeAttribsDialog();
    void displayDelaunay();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    JMeshPtr medialAxis;
    JMeshPtr delmesh;

    boost::scoped_ptr<JNodeAttributesDialog> nodeAttribDialog;
    boost::scoped_ptr<JEdgeAttributesDialog> edgeAttribDialog;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
