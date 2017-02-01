#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_InstantMeshDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshRenderDialog.hpp"

class JInstantMeshDialog : public QDialog, public Ui::InstantMeshDialog {
    Q_OBJECT

public:
    JInstantMeshDialog( QWidget *parent = 0);
    ~JInstantMeshDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMeshType( int t)  { meshtype = t; }

protected:
    void showEvent( QShowEvent *event);

private slots:
    void genMesh();
    void rejectMesh();
    void openRenderMeshDialog();
    void showOrgBoundary();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    int meshtype = 4;
    int inEulerCharacteristic;

    JMeshPtr mesh;
    JMeshPtr newMesh;
    boost::scoped_ptr<JMeshRenderDialog> meshRenderDialog;

    void setMesh( const JMeshPtr &m);
    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
