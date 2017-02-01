#pragma once

#include <vector>
#include <QPalette>
#include <qwt/qwt_plot_histogram.h>
#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_magnifier.h>
#include <qwt/qwt_plot_panner.h>


using namespace std;

class JHistogram: public QwtPlotHistogram
{
public:
    JHistogram( const QString &, const QColor & );
    ~JHistogram() {
        binData.clear();
    }

    void setColor( const QColor & );
    void setValues( const vector<double> &v, int &numBins);

private:
    vector<double> binData;
};

