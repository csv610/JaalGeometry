#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshHolesFillDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshHolesFiller.hpp"

class JMeshHolesFillDialog : public QDialog, public Ui::MeshHolesFillDialog {
    Q_OBJECT

public:
    JMeshHolesFillDialog( QWidget *parent = 0);
    ~JMeshHolesFillDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void searchHoles();
    void fillOne();
    void fillAll();
    void refineHole();
    void smoothHole();

    void closeDialog();
private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    QStandardItemModel *model;

    JMeshPtr mesh;
    vector<JEdgeSequence> contours;
    JMeshHolesFiller  holeFiller;

    void init();
    void fillTable();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
