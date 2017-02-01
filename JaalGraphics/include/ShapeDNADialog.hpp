#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_ShapeDNADialog.hpp"
#include "MeshViewer.hpp"

#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_plot_grid.h>

class JShapeDNADialog : public QDialog, public Ui::ShapeDNADialog {
    Q_OBJECT

public:
    JShapeDNADialog( QWidget *parent = 0);
    ~JShapeDNADialog();

    void setMesh( JMeshPtr m) {
        mesh = m;
    }

private slots:
    void getEigenSpectrum();

private:
    JMeshPtr mesh;
    EigenColor  *eigenColor;
    QwtPlotGrid  *grid;
    QwtPlotCurve *plotcurve;
    vector<double> xdata, eigenvalues, eigenvector;

    void makeConnections();
    int  isTriMesh();
    int  openEigFile();

};
