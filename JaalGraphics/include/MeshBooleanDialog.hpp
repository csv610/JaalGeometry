#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshBooleanDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshBoolean.hpp"
#include "ObjectsListDialog.hpp"

class JMeshBooleanDialog : public QDialog, public Ui::MeshBooleanDialog {
    Q_OBJECT

public:
    JMeshBooleanDialog( QWidget *parent = 0);
    ~JMeshBooleanDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &mA, const JMeshPtr &mB);

private slots:
    void showEvent(QShowEvent *e);

    void applyOp();
    void loadMeshA();
    void loadMeshB();

    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    bool  selectMesh[2];

    JMeshPtr meshA, meshB, meshC;
    boost::shared_ptr<JMeshBoolean> meshBoolean;
    boost::scoped_ptr<JObjectsListDialog> objectsListDialog;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
