#pragma once

#include <QDialog>

#include "Ui_MeshGeodesicsDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshGeodesics.hpp"
#include "NodeAttributesDialog.hpp"
#include "EdgeAttributesDialog.hpp"
#include "DDG_MeshGeodesics.hpp"
#include <igl/jet.h>

class JMeshGeodesicsDialog : public QDialog, public Ui::MeshGeodesicDialog {
    Q_OBJECT

public:
    JMeshGeodesicsDialog( QWidget *parent = 0);
    ~JMeshGeodesicsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

protected:
    void  mouseReleaseEvent( QMouseEvent *e);
    void  keyPressEvent( QKeyEvent *e);

private slots:
    void getDijkstraPath();
    void checkState();
    void closeDialog();
    void openNodeAttribDialog();
    void openEdgeAttribDialog();
    void startNewPath();
    void deleteLastSegment();
    void deleteAll();
    void heatFlow();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshEntityPickerPtr picker;

    JMeshPtr mesh;
    boost::scoped_ptr<JMeshGeodesics>  jgeodesic;

    boost::scoped_ptr<JNodeAttributesDialog> nodeAttribDialog;
    boost::scoped_ptr<JEdgeAttributesDialog> edgeAttribDialog;

    bool oneSrcDst;
    vector<double> dist;
    JNodeSet srcSet, srcNodes;
    JNodeSet dstSet, dstNodes;
    vector<JEdgeSequence> pathEdges;
    JNodeSequence  pickedNodes;

    void donePicking();
    void setSrcDstColor( const JNodePtr &v, int src_or_dst);
    void setPathColor( const JEdgePtr &e);
    void getPath( const JNodePtr &src, const JNodePtr &dst);

    void init();
    void makeConnections();
    void setMesh( const string &n);
    void setNodesDistanceColor();
};

