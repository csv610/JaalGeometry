#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MixedIntegerQuadsDialog.hpp"
#include "MeshViewer.hpp"

#include "MixedIntegerQuads.hpp"

class JMixedIntegerQuadsDialog : public QDialog, public Ui::MixedIntegerQuadsDialog {
    Q_OBJECT


public:
    JMixedIntegerQuadsDialog( QWidget *parent = 0);
    ~JMixedIntegerQuadsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void generateQuads();
    void checkDisplay();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    boost::scoped_ptr<JMixedIntegerQuads>  mixedIntegerQuads;

    JMeshPtr mesh;
    JMeshPtr uvMesh;

    JMeshPtr crossField, bisectorField, bisectorCombinedField;
    JFaceSequence singularFaces;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
