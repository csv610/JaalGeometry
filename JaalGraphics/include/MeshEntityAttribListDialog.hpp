#pragma once

#include <QDialog>
#include <QStringListModel>

#include "Ui_MeshEntityAttribListDialog.hpp"
#include "Mesh.hpp"

using namespace Jaal;

class JMeshEntityAttribListDialog : public QDialog, public Ui::MeshEntityAttribListDialog {
    Q_OBJECT

public:
    JMeshEntityAttribListDialog( QWidget *parent = 0);
    ~JMeshEntityAttribListDialog();

    void setMesh( const JMeshPtr &m, int e) {
        mesh = m;
        entity = e;
        init();
    }

private slots:
    void init();

private:
    int  entity;
    JMeshPtr mesh;
    QStringListModel *model;

    void makeConnections();
};


