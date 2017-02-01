#pragma once

#include <QDialog>
#include <QColorDialog>

#include <boost/filesystem.hpp>
#include <ostream>

#include "Ui_CurveShorteningFlowDialog.hpp"
#include "MeshViewer.hpp"
#include "CurveShorteningFlow.hpp"

class JCurveShorteningFlowDialog : public QDialog, public Ui::CurveShorteningFlowDialog {
    Q_OBJECT

public:
    JCurveShorteningFlowDialog( QWidget *parent = 0);
    ~JCurveShorteningFlowDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void startflow();
    void getOriginal();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;
    int      currStep;

    vector<double> orgCoords;
    vector<size_t> l2g;

    boost::scoped_ptr<JCurveShorteningFlow>  csflow;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
