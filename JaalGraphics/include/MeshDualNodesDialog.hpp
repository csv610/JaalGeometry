#pragma once

#include <QDialog>

#include "Ui_MeshDualNodesDialog.hpp"

#include "MeshDualViewer.hpp"
#include "NodeAttributesDialog.hpp"

class JMeshEntityAttribListDialog;

class JMeshDualNodesDialog : public QDialog, public Ui::MeshDualNodesDialog {
    Q_OBJECT

public:
    JMeshDualNodesDialog( QWidget *parent = 0);

    void setMeshViewer( JMeshDualViewer::shared_ptr v) {
        dualViewer = v;
        init();
    }

private slots:
    void keyPressEvent( QKeyEvent *e);

    void openAttribDialog();
    void checkNodes();
    void deleteNodes();

private:
    JMeshPtr dualGraph;

    JMeshDualViewer::shared_ptr dualViewer;
    boost::scoped_ptr<JNodeAttributesDialog> attribDialog;

    void makeConnections();
    void init();
};
