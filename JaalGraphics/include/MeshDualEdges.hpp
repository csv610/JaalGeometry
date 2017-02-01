#pragma once

#ifndef JMESHDUALNODE_H
#define JMESHDUALNODE_H

#include <QDialog>

#include "Ui_MeshDualNodesDialog.hpp"

#include "MeshDualViewer.hpp"
#include "NodeAttributesDialog.hpp"

class JMeshEntityAttribListDialog;

class JMeshDualNodesDialog : public QDialog, public Ui::MeshDualNodesDialog {
    Q_OBJECT

public:
    JMeshDualNodesDialog( QWidget *parent = 0);

    void setMeshViewer( JMeshDualViewer *v) {
        dualViewer = v;
    }

private slots:
    void keyPressEvent( QKeyEvent *e);
    void openAttribDialog();
    void checkNodes();

private:
    JMeshDualViewer *dualViewer;
    Mesh *dualGraph;

    DrawNode  *drawNode;
    MeshNodeColor *nodeColor;
    JNodeAttributesDialog *attribDialog;

    void makeConnections();
};

#endif
