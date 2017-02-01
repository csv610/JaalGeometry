#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshComponentsDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshTopology.hpp"
#include "MeshBoolean.hpp"

////////////////////////////////////////////////////////////////////////////////

class JMeshComponentsDialog : public QDialog, public Ui::MeshComponentsDialog {
    Q_OBJECT

public:
    JMeshComponentsDialog( QWidget *parent = 0);
    ~JMeshComponentsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void searchComponents();
    void removeComponent();
    void mergeComponents();
    void select2Components();
    void mergeAll();
    void displayComponent();

    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;
    vector<JColor> colors;

    boost::scoped_ptr<JMeshComponents> jc;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
