#include "MeshEdgesPartitionsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshEdgesPartitionsDialog :: JMeshEdgesPartitionsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshEdgesPartitionsDialog :: ~JMeshEdgesPartitionsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesPartitionsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    qwtPlot->setCanvasBackground( Qt::white);
    plotcurve.reset( new QwtPlotCurve());
    qwtPlot->setAxisTitle( QwtPlot::xBottom, "InterfaceID ");
    qwtPlot->setAxisTitle( QwtPlot::yLeft, "#Segments");
    plotcurve->setSymbol( new QwtSymbol( QwtSymbol::Ellipse, QColor(Qt::blue),
                                         QPen( Qt::red ), QSize( 5, 5 ) ) );
    plotcurve->setPen( QColor( Qt::black ) );
    plotcurve->setStyle( QwtPlotCurve::Dots );

    magnifier.reset(new QwtPlotMagnifier(qwtPlot->canvas()));
    magnifier->setMouseButton(Qt::MidButton);

    panner.reset(new QwtPlotPanner(qwtPlot->canvas()));
    panner->setMouseButton(Qt::LeftButton);

    grid.reset(new QwtPlotGrid());
    grid->setPen( QPen( Qt::black, 0, Qt::DotLine ) );
    grid->attach( qwtPlot );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesPartitionsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    initMesh();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesPartitionsDialog :: initMesh()
{
    meshpart.reset( new JMeshPartitioner);
    meshpart->setMesh(mesh);

    int err = meshpart->searchInterfaces();
    if( err ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning: Elements have not been partitioned: interface search failed ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    int numInterfaces = meshpart->getNumInterfaces();
    numInterfacesLineEdit->setText( QString::number(numInterfaces) );

    vector<double> xData(numInterfaces);
    vector<double> yData(numInterfaces);
    JEdgeSequence edges;
    for( int i = 0; i < numInterfaces; i++) {
        meshpart->getInterface(i, edges);
        xData[i] = i;
        yData[i] = edges.size();
    }

    boost::sort(yData);
    double maxVal = *boost::max_element(yData);
    plotcurve->attach( nullptr );
    qwtPlot->updateCanvasMargins();
    qwtPlot->setAxisScale( QwtPlot::xBottom, 0, numInterfaces);
    qwtPlot->setAxisScale( QwtPlot::yLeft, 0, maxVal);
    plotcurve->setSamples( &xData[0], &yData[0], numInterfaces );
    plotcurve->attach( qwtPlot );
    qwtPlot->replot();
    viewPartitionSpinBox->setMaximum(numInterfaces-1);

    displayEdges();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesPartitionsDialog :: displayEdges()
{
    JEdgeRenderPtr   eAttrib;
    JEdgeSequence edges;

    bool val = !displayInterfaceOnlyCheckBox->isChecked();

    JColor white = JEntityColor::getColor("White");
    JColor black = JEntityColor::getColor("Black");

    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        int err = edge->getAttribute("Render", eAttrib);
        eAttrib->display = 0;
//        eAttrib->color   = black;
        if( edge->hasAttribute("Interface") )  {
             eAttrib->display = 1;
//           eAttrib->color   = white;
        }
    }

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesPartitionsDialog :: viewPartition()
{
    if( meshpart == nullptr) return;
    JEdgeRenderPtr eAttrib;

    for( const JEdgePtr &edge : highlightEdges) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->lineWidth = 1.0;
    }

    int id = viewPartitionSpinBox->value();
    meshpart->getInterface(id, highlightEdges);
    if( highlightEdges.empty() ) return;

    numSegmentsLineEdit->setText( QString::number(highlightEdges.size()));

    for( const JEdgePtr &edge : highlightEdges) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->lineWidth = 5.0;
    }
    const Point3D p3d = highlightEdges[0]->getNodeAt(0)->getXYZCoords();
    qglviewer::Vec vec;
    vec[0] = p3d[0];
    vec[1] = p3d[1];
    vec[2] = p3d[2];
    viewManager->setSceneCenter(vec);
    meshViewer->updateBuffers(mesh);

}
///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesPartitionsDialog :: openEdgeAttribsDialog()
{
    if( attribsDialog.get() == nullptr )
        attribsDialog.reset(new JEdgeAttributesDialog( this ));

    attribsDialog->setViewManager(viewManager);
    attribsDialog->setMesh(mesh);

    JEdgeSequence eseq;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &e = mesh->getEdgeAt(i);
        if( e->hasAttribute("Interface") )  eseq.push_back(e);
    }
    attribsDialog->setEdges(eseq);
    attribsDialog->show();
    this->hide();

}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesPartitionsDialog :: closeDialog()
{
    parentWidget()->show();
    close();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesPartitionsDialog :: makeConnections()
{
    CheckBox( displayInterfaceOnlyCheckBox, [=] {displayEdges();});
    SpinBoxi( viewPartitionSpinBox, [=] {viewPartition();});
    PushButton( attribsPushButton, [=] {openEdgeAttribsDialog();});
    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
