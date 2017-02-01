#pragma once

#include "Ui_MeshRenderDialog.hpp"

#include <QDialog>
#include <QColorDialog>
#include <QMessageBox>

#include "MeshViewer.hpp"
#include "MeshNodesDialog.hpp"
#include "MeshEdgesDialog.hpp"
#include "MeshFacesDialog.hpp"
#include "MeshCellsDialog.hpp"
#include "ObjectsListDialog.hpp"
#include "SuggestiveContoursDialog.hpp"

class JMeshRenderDialog : public QDialog, public Ui::MeshRenderDialog {
    Q_OBJECT

public:
    JMeshRenderDialog( QWidget *parent = 0);
    ~JMeshRenderDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

protected:
    void  showEvent(QShowEvent *e);
    void  keyPressEvent( QKeyEvent *e);

private slots:
    void  openMaterialDialog();
    void  openMeshlistDialog();

    void  openBoxColorDialog();

    void  setEnclosure();

    void  resetCamera();
    void  fitBoundSphere();
    void  fitBoundBox();

    void  openNodesDialog();
    void  openEdgesDialog();
    void  openFacesDialog();
    void  openCellsDialog();
    void  setPOVRayScene();

    void  closeDialog();

private:
    JaalViewer *viewManager;

    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;

    boost::scoped_ptr<JObjectsListDialog> meshlistDialog;
    boost::scoped_ptr<JMeshNodesDialog> meshNodesDialog;
    boost::scoped_ptr<JMeshEdgesDialog> meshEdgesDialog;
    boost::scoped_ptr<JMeshFacesDialog> meshFacesDialog;
    boost::scoped_ptr<JMeshCellsDialog> meshCellsDialog;

    void init();
    void makeConnections();
    void warnMessage();
};

