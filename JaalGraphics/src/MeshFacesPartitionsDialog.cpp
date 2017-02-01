#include "MeshFacesPartitionsDialog.hpp"
#include "MeshPartitionDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshFacesPartitionsDialog :: JMeshFacesPartitionsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshFacesPartitionsDialog :: ~JMeshFacesPartitionsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFacesPartitionsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    qwtPlot->setCanvasBackground( Qt::white);
    plotcurve.reset( new QwtPlotCurve());
    qwtPlot->setAxisTitle( QwtPlot::xBottom, "PartitionID ");
    qwtPlot->setAxisTitle( QwtPlot::yLeft, "#Faces");
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

void JMeshFacesPartitionsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    initMesh();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshFacesPartitionsDialog :: initMesh()
{
    meshpart.reset( new JMeshPartitioner);
    meshpart->setMesh(mesh);

    int numParts = meshpart->getNumPartitions();

    numPartitionsLineEdit->setText( QString::number(numParts) );

    if( numParts == 0) return;
    vector<double> xData(numParts);
    vector<double> yData(numParts);
    JEdgeSequence edges;
    for( int i = 0; i < numParts; i++) {
        JMeshPtr patch = meshpart->getSubMesh(i);
        assert( patch ) ;
        xData[i] = i;
        yData[i] = patch->getSize(2);
    }

    boost::sort(yData);
    double maxVal = *boost::max_element(yData);
    plotcurve->attach( nullptr );
    qwtPlot->updateCanvasMargins();
    qwtPlot->setAxisScale( QwtPlot::xBottom, 0, numParts);
    qwtPlot->setAxisScale( QwtPlot::yLeft, 0, maxVal);
    plotcurve->setSamples( &xData[0], &yData[0], numParts );
    plotcurve->attach( qwtPlot );
    qwtPlot->replot();
    displayPatchSpinBox->setMaximum(numParts-1);

    JFacePartitionColor faceColor;
    faceColor.setMesh(mesh);
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshFacesPartitionsDialog :: checkDisk()
{
    int numParts = meshpart->getNumPartitions();
    if( numParts == 0) return;
    QMessageBox msg;
    int ret;

    int nCount = 0;
    for(int i = 0; i < numParts; i++) {
        JMeshPtr patch = meshpart->getSubMesh(i);
        if( !patch->getTopology()->isDisk() ) {
            nCount++;
            msg.setIcon(QMessageBox::Warning);
            QString qstr = QString("Patch ") + QString::number(i) + QString(" is not a topological disk");
            msg.setText(qstr);
            msg.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel);
            ret = msg.exec();
            if( ret == QMessageBox::Cancel ) return;
        }
    }

    if( nCount == 0) {
        msg.setIcon(QMessageBox::Information);
        QString qstr = QString("All patches have disk topology");
        msg.setText(qstr);
        msg.setStandardButtons( QMessageBox::Ok );
        ret = msg.exec();
        if( ret == QMessageBox::Ok) return;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshFacesPartitionsDialog :: displayPatch()
{
    if( mesh == nullptr) return;

    JFaceRenderPtr   fAttrib;
    int val = 0;
    size_t numFaces = mesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        int err = face->getAttribute("Render", fAttrib);
        if( !err )  fAttrib->display = val;
    }

    int id = displayPatchSpinBox->value();
    JMeshPtr patch = meshpart->getSubMesh(id);
    numFaces = patch->getSize(2);
    val = 1;
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = patch->getFaceAt(i);
        int err = face->getAttribute("Render", fAttrib);
        if( !err )  fAttrib->display = val;
    }
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFacesPartitionsDialog :: viewAll()
{
    int val = 1;
    size_t numFaces = mesh->getSize(2);
    JFaceRenderPtr   fAttrib;
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        int err = face->getAttribute("Render", fAttrib);
        if( !err )  fAttrib->display = val;
    }
    meshViewer->updateBuffers( mesh );

}
///////////////////////////////////////////////////////////////////////////////
void JMeshFacesPartitionsDialog :: randomColors()
{
    if( mesh == nullptr) return;
    JFacePartitionColor faceColor;
    faceColor.setMesh(mesh);
    meshViewer->updateBuffers( mesh );
}
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesPartitionsDialog :: closeDialog()
{
    parentWidget()->show();
    close();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesPartitionsDialog :: makeConnections()
{
    SpinBoxi( displayPatchSpinBox, [=] {displayPatch(); });
    PushButton( viewAllPushButton, [=] {viewAll();});
    PushButton( applyRandomColorsPushButton, [=] {randomColors(); });
    PushButton( checkDiskPushButton, [=] {checkDisk(); });
    PushButton( closePushButton,     [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
