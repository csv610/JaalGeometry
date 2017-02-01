#pragma once

#include <QDialog>

#include "Ui_FaceAttributesDialog.hpp"
#include "MeshViewer.hpp"

using namespace Jaal;

class JFaceAttributesDialog : public QDialog, public Ui::FaceAttributesDialog
{
    Q_OBJECT

public:
    JFaceAttributesDialog( QWidget *parent = 0);

    void setViewManager(JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh(const JMeshPtr &m);

    void updateDefaults( bool b) {
        update_default_values = b;
    }
    void setFaces( JFaceSequence &es);

private slots:
    void setColor();
    void setMaterial();
    void checkDisplay();
//   void setOffset();

    void keyPressEvent( QKeyEvent *e);

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr     mesh;

    JFaceSequence faces;
    bool update_default_values;

    void init();
    void makeConnections();
};
