#pragma once

#include <QDialog>
#include <QColorDialog>

#include "QwtHistogram.hpp"

#include "Ui_MeshTopologyQualityDialog.hpp"
#include "MeshViewer.hpp"

class JMeshTopologyQualityDialog : public QDialog, public Ui::MeshTopologyQualityDialog {
    Q_OBJECT
public:
    JMeshTopologyQualityDialog( QWidget *parent = 0);
    ~JMeshTopologyQualityDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void closeDialog();
    void setNodesDomain();
    void displayDegrees();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;

    QStandardItemModel *model;
    boost::scoped_ptr<QwtPlotGrid>  grid;
    boost::scoped_ptr<JHistogram>   qhistogram;
    boost::scoped_ptr<QwtPlotCurve> plotcurve;
    boost::scoped_ptr<JNodeDegreeColor>  nodeColor;
    boost::scoped_ptr<QwtPlotMagnifier>  magnifier;
    boost::scoped_ptr<QwtPlotPanner>     panner;

    int nodepos;
    int maxDegree;

    map<int,size_t> degreeCount;
    vector<double> xData, yData;

    void init();
    void getMeasure();
    void makeConnections();
};
