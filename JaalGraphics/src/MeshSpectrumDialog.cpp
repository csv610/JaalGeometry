#include "MeshSpectrumDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshSpectrumDialog :: JMeshSpectrumDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    meshSpectrum.reset( new JMeshSpectrum);
}

///////////////////////////////////////////////////////////////////////////////

JMeshSpectrumDialog :: ~JMeshSpectrumDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSpectrumDialog :: init()
{
    qwtPlot->setCanvasBackground( Qt::white);
    grid.reset(new QwtPlotGrid());
    grid->setPen( QPen( Qt::black, 0, Qt::DotLine ) );
    grid->attach( qwtPlot );

    plotcurve.reset( new QwtPlotCurve());
    magnifier.reset(new QwtPlotMagnifier(qwtPlot->canvas()));
    magnifier->setMouseButton(Qt::MidButton);
    panner.reset(new QwtPlotPanner(qwtPlot->canvas()));
    panner->setMouseButton(Qt::LeftButton);
    plotcurve->setSymbol( new QwtSymbol( QwtSymbol::Ellipse, QColor(Qt::blue),
                                         QPen( Qt::red ), QSize( 5, 5 ) ) );
    plotcurve->setPen( QColor( Qt::black ) );
    plotcurve->setStyle( QwtPlotCurve::Dots );

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSpectrumDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    meshSpectrum->setMesh(mesh);

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->setSurfaceShade(JRender::SMOOTH_SHADE);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSpectrumDialog :: genVectors()
{
    JWaitCursor waitCursor;
    waitCursor.start();

    int n = genEigenVecSpinBox->value();
    viewEigenVecSpinBox->setMaximum(n-1);
    meshSpectrum->genEigenVectors(n);

}
///////////////////////////////////////////////////////////////////////////////
void JMeshSpectrumDialog :: viewVector()
{
    int id = viewEigenVecSpinBox->value();
    vector<double> ev = meshSpectrum->getEigenVector(id);
    if( ev.empty() ) return;

    size_t numnodes = mesh->getSize(0);

    Eigen::VectorXd K;
    K.resize(numnodes);
    for( size_t i = 0; i < numnodes; i++)
        K[i] = ev[i];

    Eigen::MatrixXd colors;
    igl::jet(K, true, colors);

    JNodeRenderPtr vAttrib;
    JColor  color;
    int     index = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) {
            v->getAttribute("Render", vAttrib);
            color[0] = colors.coeff(index, 0);
            color[1] = colors.coeff(index, 1);
            color[2] = colors.coeff(index, 2);
            color[3] = 1.0;
            vAttrib->color = color;
            index++;
        }
    }
    meshViewer->updateBuffers(mesh);

    int nsize = ev.size();
    vector<double> xData;
    xData.resize( nsize);
    for( size_t i = 0; i < nsize; i++)
        xData[i] = 1.0*i;

    sort( ev.begin(), ev.end() );
    double minVal  = *boost::min_element( ev );
    double maxVal  = *boost::max_element( ev );
    plotcurve->attach( nullptr );

    qwtPlot->updateCanvasMargins();
    plotcurve->setRawSamples( &xData[0], &ev[0], nsize );
    qwtPlot->setAxisScale( QwtPlot::xBottom, 0, nsize);
    qwtPlot->setAxisScale( QwtPlot::yLeft, minVal, maxVal);
    plotcurve->attach( qwtPlot );
    qwtPlot->replot();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSpectrumDialog :: makeConnections()
{
    PushButton( genEigenVectorsPushButton, [=] {genVectors();});
    PushButton( viewPushButton, [=] {viewVector();});
    SpinBoxi(  viewEigenVecSpinBox, [=] {viewVector();});
    PushButton( closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
