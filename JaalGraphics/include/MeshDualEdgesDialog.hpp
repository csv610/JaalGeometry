#pragma once

#include <QDialog>

#include "Ui_MeshDualEdgesDialog.hpp"

#include "MeshDualViewer.hpp"
#include "EdgeAttributesDialog.hpp"

class JMeshEntityAttribListDialog;

class JMeshDualEdgesDialog : public QDialog, public Ui::MeshDualEdgesDialog {
    Q_OBJECT

public:
    JMeshDualEdgesDialog( QWidget *parent = 0);

    void setMeshViewer( boost::shared_ptr<JMeshDualViewer> &v) {
        dualViewer = v;
        init();
    }

private slots:
    void keyPressEvent( QKeyEvent *e);
    void openAttribDialog();
    void checkEdges();
    void deleteEdges();

private:
    JMeshPtr dualGraph;
    boost::shared_ptr<JMeshDualViewer> dualViewer;
    boost::scoped_ptr<JEdgeAttributesDialog> attribDialog;

    void makeConnections();
    void init();
};
