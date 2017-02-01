#pragma once

#include <QDialog>
#include "Ui_MeshSamplesDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshSamples.hpp"

class JMeshSamplesDialog : public QDialog, public Ui::MeshSamplesDialog {
    Q_OBJECT

public:
    JMeshSamplesDialog( QWidget *parent = 0);
    ~JMeshSamplesDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void getSamples();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;
    boost::scoped_ptr<JMeshSamples> meshSamples;

    int  algo;
    int  numSamples;

    void init();
    void makeConnections();
    void getNodeSamples();
    void getEdgeSamples();
    void getFaceSamples();
    void getCellSamples();

};
////////////////////////////////////////////////////////////////////////////////
