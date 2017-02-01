#include "IsotropicTriMeshDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JIsotropicTriMeshDialog :: JIsotropicTriMeshDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    qwtPlot->setCanvasBackground( Qt::white);
    grid.reset(new QwtPlotGrid());
    grid->setPen( QPen( Qt::black, 0, Qt::DotLine ) );
    grid->attach( qwtPlot );

    plotcurve.reset( new QwtPlotCurve());
    magnifier.reset(new QwtPlotMagnifier(qwtPlot->canvas()));
    magnifier->setMouseButton(Qt::MidButton);

    panner.reset(new QwtPlotPanner(qwtPlot->canvas()));
    panner->setMouseButton(Qt::LeftButton);
}

///////////////////////////////////////////////////////////////////////////////

JIsotropicTriMeshDialog :: ~JIsotropicTriMeshDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JIsotropicTriMeshDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JIsotropicTriMeshDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////

void JIsotropicTriMeshDialog :: getEdgesLength()
{
    if( mesh == nullptr) return;

    boost::scoped_ptr<JMeshQuality> meshQual;
    meshQual.reset( new JMeshQuality);
    meshQual->setMesh(mesh);

    int pos = JMeshEntity::ANY_ENTITY;
    qData = meshQual->getEdgesQuality( JMeshQuality::EDGE_LENGTH, pos, 1);

    size_t nsize = qData.size();
    if( nsize == 0) return;

    xData.clear();
    yData.clear();

    if( nsize ) {
        xData.resize(nsize);
        yData.resize(nsize);
    }

    for( size_t i = 0; i < nsize; i++) {
        xData[i] = 1.0*i;
        yData[i] = qData[i];
    }

    double minVal  = *boost::min_element( yData );
    double maxVal  = *boost::max_element( yData );
    plotcurve->attach( nullptr );
    double avg     = 0.5*(maxVal + minVal );
    edgeLengthLineEdit->setText( QString::number(avg));

    qwtPlot->updateCanvasMargins();
    qwtPlot->setAxisTitle( QwtPlot::yLeft, "EdgeLength");
    plotcurve->setSymbol( new QwtSymbol( QwtSymbol::Ellipse, QColor(Qt::blue),
                                         QPen( Qt::red ), QSize( 5, 5 ) ) );
    plotcurve->setPen( QColor( Qt::black ) );
    plotcurve->setStyle( QwtPlotCurve::Dots );
    plotcurve->setRawSamples( &xData[0], &yData[0], nsize );

    plotcurve->attach( qwtPlot );
    qwtPlot->replot();

}

///////////////////////////////////////////////////////////////////////////////

void JIsotropicTriMeshDialog :: genMesh()
{
    if( mesh == nullptr) return;
    JIsotropicMesh isomesh;
    isomesh.setMesh(mesh);
    QString qstr = edgeLengthLineEdit->text();
    double elen = qstr.toDouble();
    if( elen < 1.0E-10) return;

    JWaitCursor wCursor;
    wCursor.start();
    isomesh.setEdgeLength(elen);
    JMeshPtr newmesh = isomesh.getMesh();
    meshViewer->addObject(newmesh);
}

///////////////////////////////////////////////////////////////////////////////
void JIsotropicTriMeshDialog :: setSharpEdges()
{
}
void JIsotropicTriMeshDialog :: closeDialog()
{
    this->close();
    parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JIsotropicTriMeshDialog :: makeConnections()
{
    PushButton( sharpEdgesPushButton, [=] {setSharpEdges();});
    PushButton( getEdgesLengthPushButton, [=] {getEdgesLength();});
    PushButton( genMeshPushButton, [=] {genMesh();});
    PushButton( closePushButton,   [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
