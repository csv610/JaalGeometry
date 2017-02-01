#include "QwtHistogram.hpp"

///////////////////////////////////////////////////////////////////////////////
JHistogram::JHistogram( const QString &title, const QColor &symbolColor ):
    QwtPlotHistogram( title )
{
    setStyle( QwtPlotHistogram::Columns );

    setColor( symbolColor );
}

///////////////////////////////////////////////////////////////////////////////

void JHistogram::setColor( const QColor &symbolColor )
{
    QColor color = symbolColor;
    color.setAlpha( 180 );

    setPen( QPen( Qt::black ) );
    setBrush( QBrush( color ) );

    QwtColumnSymbol *symbol = new QwtColumnSymbol( QwtColumnSymbol::Box );
    symbol->setFrameStyle( QwtColumnSymbol::Raised );
    symbol->setLineWidth( 2 );
    symbol->setPalette( QPalette( color ) );
    setSymbol( symbol );
}

///////////////////////////////////////////////////////////////////////////////

void JHistogram::setValues( const vector<double> &values, int &numBins )
{
    size_t numValues = values.size();

    double minVal = *min_element( values.begin(), values.end());
    double maxVal = *max_element( values.begin(), values.end());

    if( maxVal != minVal ) {
        double dx = (maxVal - minVal)/(double)numBins;
        binData.resize(numBins);
        for( int i = 0; i < numBins; i++)
            binData[i] = 0;

        for( size_t i = 0; i < numValues; i++) {
            int bid =  (values[i] - minVal)/dx;
            binData[bid]++;
        }

    } else {
        numBins = 1;
        binData.resize(1);
        binData[0] = numValues;
    }

    QVector<QwtIntervalSample> samples( numBins );
    for ( int i = 0; i < numBins; i++ ) {
        QwtInterval interval( double( i ), i + 1.0 );
        interval.setBorderFlags( QwtInterval::ExcludeMaximum );

        samples[i] = QwtIntervalSample( binData[i], interval );
    }
    setData( new QwtIntervalSeriesData( samples ) );
}
///////////////////////////////////////////////////////////////////////////////
