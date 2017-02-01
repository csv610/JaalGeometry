#pragma once

#include <QDialog>
#include <QColorDialog>
#include <igl/jet.h>
#include <igl/parula.h>

#include "Ui_MeshUntangleDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshUntangle.hpp"
#include "LoadNewDataDialog.hpp"
#include "MeshMapQualityDialog.hpp"

////////////////////////////////////////////////////////////////////////////////

class JMeshUntangleDialog : public QDialog, public Ui::MeshUntangleDialog {
    Q_OBJECT

public:
    static const int  JAAL_METHOD  = 0;
    static const int  SACHT_METHOD = 1;

    JMeshUntangleDialog( QWidget *parent = 0);
    ~JMeshUntangleDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);
    void setMethod( int m) {
        method = m;
    }

private slots:
    void execute();
    void loadSrcMesh();
    void setOriginal();
    void openMapQualityDialog();
    void backProjection();
    void setColor();
    void setLaplaceSmoothing();
    void closeDialog();
private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    int method = JAAL_METHOD;
    JMeshPtr srcMesh, inflatedMesh, deformedMesh;
    JMeshUntangle untangle;
    boost::scoped_ptr<JMeshMapQualityDialog> mapQualityDialog;
    std::map<JNodePtr, JColor>   prevNodeColor;
    std::map<JEdgePtr, JColor>   prevEdgeColor;
    std::map<JFacePtr, JColor>   prevFaceColor;
    JColor highlightColor;
    bool initialized;

    vector<double> orgCoords;
    vector<size_t> l2g;

    void init();
    void initMesh();
    void makeConnections();
    void displayOffset();
    void useSachtMethod();
    void useJaalMethod();
    void displayInverted( const JMeshPtr &m);
};
////////////////////////////////////////////////////////////////////////////////
