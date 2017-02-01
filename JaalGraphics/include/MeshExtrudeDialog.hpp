#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshExtrudeDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshExtrude.hpp"

class JMeshExtrudeDialog : public QDialog, public Ui::MeshExtrudeDialog {
    Q_OBJECT

public:
    JMeshExtrudeDialog( QWidget *parent = 0);
    ~JMeshExtrudeDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void genMesh();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    boost::shared_ptr<JMeshExtrude> meshExtrude;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
