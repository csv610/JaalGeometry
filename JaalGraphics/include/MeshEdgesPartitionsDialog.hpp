#pragma once

#include <QDialog>
#include <QColorDialog>

#include "MeshViewer.hpp"
#include "Ui_MeshEdgesPartitionsDialog.hpp"
#include "EdgeAttributesDialog.hpp"

#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_magnifier.h>
#include <qwt/qwt_plot_panner.h>

class JMeshEdgesPartitionsDialog : public QDialog, public Ui::MeshEdgesPartitionsDialog {
    Q_OBJECT

public:
    JMeshEdgesPartitionsDialog( QWidget *parent = 0);
    ~JMeshEdgesPartitionsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void viewPartition();
    void openEdgeAttribsDialog();
    void displayEdges();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;
    JEdgeSequence highlightEdges;

    boost::scoped_ptr<JMeshPartitioner> meshpart;
    boost::scoped_ptr<JEdgeAttributesDialog> attribsDialog;

    boost::scoped_ptr<QwtPlotCurve> plotcurve;
    boost::scoped_ptr<QwtPlotGrid>  grid;
    boost::scoped_ptr<QwtPlotMagnifier> magnifier;
    boost::scoped_ptr<QwtPlotPanner>    panner;

    void init();
    void initMesh();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
