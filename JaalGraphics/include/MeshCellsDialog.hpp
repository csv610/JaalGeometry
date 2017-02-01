#pragma once

#include <QDialog>

#include "MeshViewer.hpp"
#include "Ui_MeshCellsDialog.hpp"
#include "MeshWavefrontsDialog.hpp"
#include "MeshSlicerDialog.hpp"

class MeshEntityAttributesDialog;

class JMeshCellsDialog : public QDialog, public Ui::MeshCellsDialog {
    Q_OBJECT

public:
    JMeshCellsDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh(const JMeshPtr &m);

private slots:
    void closeDialog();
    void checkState();
    /*
         void explode();
         void showOne();
         void nextCell();
         void getCellID();
         void setDefaultColor();
    */

    void setNumVisible();
    void checkDisplay();
    void saveSurfmesh();
    void reverseAll();

    void openMeshSlicerDialog();
    void openAttribListDialog();
    void openWavefrontDialog();
    void deleteAll();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh, exploded_mesh;

    boost::scoped_ptr<JMeshWavefrontsDialog> wavefrontDialog;
    boost::scoped_ptr<JMeshSlicerDialog> meshSlicerDialog;
    boost::scoped_ptr<JMeshEntityAttribListDialog> attribListDialog;

    size_t currCellID;
    void init();
    void activate( JCellPtr cell );
    void deactivate( JCellPtr cell );

    void makeConnections();
};

