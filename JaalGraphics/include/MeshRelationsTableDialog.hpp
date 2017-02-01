#pragma once

#include "Ui_MeshRelationsTableDialog.hpp"

#include <QDialog>
#include "MeshViewer.hpp"

class JMeshRelationsTableDialog : public QDialog, public Ui::MeshRelationsTableDialog
{
    Q_OBJECT

public:
    JMeshRelationsTableDialog( QWidget *parent = 0);
    ~JMeshRelationsTableDialog();

    void setMesh( const JMeshPtr &v) {
        mesh = v;
        init();
    }

private:
    JMeshPtr mesh;
    void  init();
    void makeConnections();
};

