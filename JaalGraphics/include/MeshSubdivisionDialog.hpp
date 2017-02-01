#pragma once

#include <QDialog>

#include "Ui_MeshSubdivisionDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshSubdivision.hpp"

class JMeshSubdivisionDialog : public QDialog, public Ui::MeshSubdivisionDialog {
    Q_OBJECT

public:
    JMeshSubdivisionDialog( QWidget *parent = 0);
    ~JMeshSubdivisionDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh(const JMeshPtr &m);

private slots:
    void getRefinedMesh();
    void closeDialog();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
//  boost::scoped_ptr<JMeshSubdivision> meshSubdiv;

    void init();
    void makeConnections();
};
