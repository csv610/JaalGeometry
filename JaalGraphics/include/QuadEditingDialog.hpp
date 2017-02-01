#pragma once

#include <QDialog>

#include "Ui_QuadEditingDialog.hpp"
#include "MeshViewer.hpp"
#include "MSTQuadMesher.hpp"
#include "MeshPartitionDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

class JQuadEditingViewer  : public JViewComponent {
public:
    JQuadEditingViewer();
    ~JQuadEditingViewer() { }

    void setViewManager(JaalViewer *v) {
        viewManager = v;
    }
    JaalViewer* getViewManager() const {
        return viewManager;
    }

    void setMesh(const JMeshPtr &m) {
        mesh = m;
    }

    void setSingularNodes( const JNodeSequence &v) {
        singularNodes = v;
    }

    void  draw();
private:
    JColor   colors[5];
    JMeshPtr  mesh;
    JaalViewer *viewManager;
    JNodeSequence singularNodes;
};

///////////////////////////////////////////////////////////////////////////////

class JQuadEditingDialog : public QDialog, public Ui::QuadEditingDialog {
    Q_OBJECT

public:
    JQuadEditingDialog( QWidget *parent = 0);
    ~JQuadEditingDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

protected:

    void  mouseReleaseEvent( QMouseEvent *e);
//    void  keyPressEvent( QKeyEvent *e);

private slots:
    void collapse();
    void closeDialog();
    void clearAll();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    boost::shared_ptr<JQuadEditingViewer>  quadEditViewer;

private:
    JMeshPtr mesh;
    JMeshEntityPickerPtr picker;
    JNodeSequence  pickedNodes;
    JMSTQuadMesher mstMesher;
//    JMotorcycleGraph   motorGraph;
    JFaceColorPtr faceColor;
    JEdgeColorPtr edgeColor;
    bool stop_picking;

    void init();
    int  isQuadMesh();
    void makeConnections();
    void getParitions();
    void displaySingularNodes();
    void assignPartitionColors();
};
////////////////////////////////////////////////////////////////////////////////

