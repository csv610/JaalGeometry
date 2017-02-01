#pragma once

#include <QDialog>
#include <QColorDialog>
#include <numeric>

#include "Ui_MeshGeometricQualityDialog.hpp"

#include "MeshViewer.hpp"
#include "MeshQuality.hpp"
#include "StatisticsDialog.hpp"

#include "QwtHistogram.hpp"

#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_magnifier.h>
#include <qwt/qwt_plot_panner.h>

class JEdgeQualityColor : public JEdgeColor {
public:

    string getName() const {
        return "EdgeQualityColor";
    }

    JEdgeQualityColor() {
        attribname = "Quality";
    }

    bool isPerEdge() const {
        return 1;
    }

    void setLowerCutoff( double l) {
        lowerVal = l;
    }

    void setUpperCutoff( double l) {
        upperVal = l;
    }

    void setAttribute( const string &s) {
        attribname = s;
    }

    int assign(const JEdgePtr &edge);

    int  operator() (const JEdgePtr &edge) {
        return assign(edge);
    }

private:
    double val, lowerVal, upperVal;
    string attribname;
};

///////////////////////////////////////////////////////////////////////////////

class JFaceQualityColor : public JFaceColor {
public:

    string getName() const {
        return "FaceQualityColor";
    }

    JFaceQualityColor() {
        attribname = "Quality";
    }

    bool isPerFace() const {
        return 1;
    }

    void setLowerCutoff( double l) {
        lowerVal = l;
    }

    void setUpperCutoff( double l) {
        upperVal = l;
    }

    void setAttribute( const string &s) {
        attribname = s;
    }

    int assign( const JFacePtr &face);

    int  operator() (const JFacePtr &f) {
        return assign(f);
    }

private:
    double val, lowerVal, upperVal;
    string attribname;
};

class JCellQualityColor : public JCellColor {
public:
    string getName() const {
        return "CellQualityColor";
    }

    bool isPerCell() const {
        return 1;
    }

    JCellQualityColor() {
        attribname = "Quality";
    }

    void setLowerCutoff( double l) {
        lowerVal = l;
    }

    void setUpperCutoff( double l) {
        upperVal = l;
    }

    void setAttribute( const string &s) {
        attribname = s;
    }

    int assign(const JCellPtr &cell);

private:
    double val, lowerVal, upperVal;
    string attribname;
};

/////////////////////////////////////////////////////////////////////////////////////

class JMeshGeometricQualityDialog : public QDialog, public Ui::MeshGeometricQualityDialog
{
    Q_OBJECT

public:
    static const int HISTOGRAM = 0;
    static const int XYPLOT    = 1;
    static const int LOGPLOT   = 2;

    JMeshGeometricQualityDialog( QWidget *parent = 0);
    ~JMeshGeometricQualityDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh(const JMeshPtr &m);

    void setHistogramBins( int n) {
        numBins = n;
    }
    void setValues();

private slots:
    void setPlotType();
    void newQuality();
    void newHistogram();
    void saveAs();
    void openStatisticsDialog();
    void closeDialog();

protected:
    void showEvent( QShowEvent *e);
    void keyPressEvent( QKeyEvent *e);

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    boost::scoped_ptr<QwtPlotGrid>  grid;
    boost::scoped_ptr<QwtPlotCurve> plotcurve;
    boost::scoped_ptr<JHistogram>   qhistogram;
    boost::scoped_ptr<JMeshQuality> meshQual;
    boost::scoped_ptr<JStatisticsDialog> statisticsDialog;
    boost::scoped_ptr<QwtPlotMagnifier>  magnifier;
    boost::scoped_ptr<QwtPlotPanner>     panner;

    map<string,int> name2id;
    vector<double>  xData, yData, qData;

    int   topDim;
    int   plotType;
    int   numBins;
    bool  savedata;
    double maxVal,minVal;
    string filename;

    void  init();
    void  addItems();
    void  makeConnections();
    void  getMeasure();
};
