#include "GraphPlot.hpp"


JGraphPlot :: JGraphPlot( const QwtPlot *p)
{
    qwtPlot = p;
    qwtPlot->setCanvasBackground( Qt::white);
    grid.reset(new QwtPlotGrid());
    grid->setPen( QPen( Qt::black, 0, Qt::DotLine ) );
    grid->attach( qwtPlot );

    qhistogram.reset(new JHistogram( "Quality", Qt::red ));
    plotcurve.reset( new QwtPlotCurve());
    magnifier.reset(new QwtPlotMagnifier(qwtPlot->canvas()));
    magnifier->setMouseButton(Qt::MidButton);

    panner.reset(new QwtPlotPanner(qwtPlot->canvas()));
    panner->setMouseButton(Qt::LeftButton);
}

