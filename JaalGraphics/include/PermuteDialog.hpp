#pragma once

#include <QDialog>

#include "Ui_PermuteDialog.hpp"

#include "MeshViewer.hpp"

#include "NodeAttributesDialog.hpp"

class JPermuteDialog : public QDialog, public Ui::PermuteDialog {
    Q_OBJECT

public:
    JPermuteDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void permute();
//   void  closeDialog();
    void  keyPressEvent( QKeyEvent *e);
private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;

    void init();
    void makeConnections();
};
