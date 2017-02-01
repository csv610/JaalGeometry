#include "ShapeDNADialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JShapeDNADialog :: JShapeDNADialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    mesh = nullptr;
//   eigenColor = new EigenColor;

    numEigenValuesLineEdit->setText( QString::number(2) );
    displayEigenVecLineEdit->setText( QString::number(2) );

    qwtPlot->setCanvasBackground( Qt::white);

    grid = new QwtPlotGrid();
//     grid->setMajPen( QPen( Qt::black, 0, Qt::DotLine ) );
    grid->attach( qwtPlot );

    qwtPlot->setAxisTitle( QwtPlot::yLeft, "Eigenvalues");
    qwtPlot->setAxisTitle( QwtPlot::xBottom, "Number");
    plotcurve  = new QwtPlotCurve();

}

///////////////////////////////////////////////////////////////////////////////

JShapeDNADialog :: ~JShapeDNADialog()
{
    delete eigenColor;
}

///////////////////////////////////////////////////////////////////////////////

int JShapeDNADialog :: isTriMesh()
{
    if( mesh == nullptr ) return 0;

    int dim = mesh->getTopology()->getDimension();
    if( dim == 2 ) {
        int etype = mesh->getTopology()->getElementsType(2);
        if( etype != 3 ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Shape DNA requires triangle mesh ");
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) {
                return 0;
            }
        }
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////
void JShapeDNADialog :: getEigenSpectrum()
{
    QString qstr;

    qstr = numEigenValuesLineEdit->text();
    int numEigs = qstr.toInt();

    qstr = displayEigenVecLineEdit->text();
    int eigVec = qstr.toInt();

    if( !isTriMesh() ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

//   mesh->saveAs("tmp.off");
    exit(0);

    ostringstream oss;
    oss << "shapeDNA --mesh tmp.off --outfile eigen.dat --degree 1 --num " << numEigs;

    string cmd = oss.str();

    int stat = system( cmd.c_str() );

    if( stat < 0) return;

    ifstream infile("eigen.dat", ios::in);
    if( infile.fail() ) return;

    eigenvalues.resize(numEigs);

    string str;
    while(!infile.eof() ) {
        infile >> str;
        if( str == "Eigenvalues:") {
            infile >> str;
            assert( str == "{");
            for( int i = 0; i < numEigs; i++) {
                infile >> eigenvalues[i];
                infile >> str;
            }
        }
    }

    plotcurve->attach( nullptr );
    xdata.resize(numEigs);
    double energy = 0.0;
    for( int i = 0; i < numEigs; i++) {
        xdata[i] = i;
        energy += fabs(eigenvalues[i] );
    }

    double minVal  = *boost::min_element( eigenvalues );
    double maxVal  = *boost::max_element( eigenvalues );

    plotcurve->setSymbol( new QwtSymbol( QwtSymbol::Ellipse, QColor(Qt::blue),
                                         QPen( Qt::red ), QSize( 5, 5 ) ) );
    plotcurve->setPen( QColor( Qt::black ) );
    plotcurve->setStyle( QwtPlotCurve::Dots );

    plotcurve->setRawSamples( &xdata[0], &eigenvalues[0], xdata.size() );
    qwtPlot->setAxisScale( QwtPlot::xBottom, 0, numEigs-1);
    qwtPlot->setAxisScale( QwtPlot::yLeft, minVal, maxVal);
    plotcurve->attach( qwtPlot );

    qwtPlot->replot();

    graphEnergyLineEdit->setText( QString::number(energy) );
}
////////////////////////////////////////////////////////////////////////


void JShapeDNADialog :: makeConnections()
{
    PushButton( applyPushButton,  [=] {getEigenSpectrum();});
    PushButton( closePushButton,  [=] {close();});

}

///////////////////////////////////////////////////////////////////////////////
