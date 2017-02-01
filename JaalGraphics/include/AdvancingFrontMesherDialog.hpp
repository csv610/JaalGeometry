#pragma once

#include <QDialog>

#include "Ui_AdvancingFrontMesherDialog.hpp"
#include "MeshViewer.hpp"

class JAdvancingFrontMesherDialog : public QDialog, public Ui::AdvancingFrontMesherDialog {
    Q_OBJECT

public:
    JAdvancingFrontMesherDialog( QWidget *parent = 0);
    ~JAdvancingFrontMesherDialog();

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }

private slots:

    void removeAll();
    void accept();
    void reject();

private:
    JaalViewer *viewManager;
    JMeshViewer *meshViewer;
    Mesh *mesh;

    void init();
    void makeConnections();
    void generateNewMesh();
};
