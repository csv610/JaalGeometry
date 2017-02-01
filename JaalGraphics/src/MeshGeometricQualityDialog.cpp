#include "MeshGeometricQualityDialog.hpp"

using namespace std;
using namespace boost;

///////////////////////////////////////////////////////////////////////////////
int JEdgeQualityColor :: assign(const JEdgePtr &edge )
{
    edge->getAttribute(attribname, val);

    JEdgeRenderPtr eAttrib;
    edge->getAttribute("Render", eAttrib);

    if( val < lowerVal ) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        color[3] = 1.0;
        eAttrib->color = color;
        return 0;
    }

    if( val > upperVal ) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 1.0;
        color[3] = 1.0;
        eAttrib->color = color;
        return 0;
    }

    color[0] = 0.0;
    color[1] = 1.0;
    color[2] = 0.0;
    color[3] = alpha;
    eAttrib->color = color;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JFaceQualityColor :: assign(const JFacePtr &face )
{
    face->getAttribute(attribname, val);

    JFaceRenderPtr fAttrib;
    face->getAttribute("Render", fAttrib);

    if( val < lowerVal ) {
        color[0] = 0.5;
        color[1] = 0.0;
        color[2] = 0.0;
        color[3] = 1.0;
        fAttrib->color = color;
        return 0;
    }

    if( val > upperVal ) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 0.5;
        color[3] = 1.0;
        fAttrib->color = color;
        return 0;
    }

    color[0] = 0.0;
    color[1] = 0.5;
    color[2] = 0.0;
    color[3] = 1.0;
    fAttrib->color = color;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JCellQualityColor :: assign( const JCellPtr &cell)
{
    cell->getAttribute(attribname, val);

    color[3] = 1.0;
    if( val < lowerVal ) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        cell->setAttribute("Color", color);
        return 0;
    }

    if( val > upperVal ) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 1.0;
        cell->setAttribute("Color", color);
        return 0;
    }
    color[0] = 0.0;
    color[1] = 1.0;
    color[2] = 0.0;
    cell->setAttribute("Color", color);
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

JMeshGeometricQualityDialog :: JMeshGeometricQualityDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    plotType = XYPLOT;

    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;

    mesh = nullptr;
    savedata = 0;
    filename = "quality.dat";

    numBinsLineEdit->setText( QString::number(100) );

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
///////////////////////////////////////////////////////////////////////////////
JMeshGeometricQualityDialog :: ~JMeshGeometricQualityDialog()
{ }
///////////////////////////////////////////////////////////////////////////////

void JMeshGeometricQualityDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshGeometricQualityDialog :: showEvent( QShowEvent *e)
{
   if( meshViewer ) setMesh( meshViewer->getCurrentMesh());
   QDialog::showEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshGeometricQualityDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    objectNameLineEdit->setText( QString(mesh->getName().c_str() ));

    topDim = mesh->getTopology()->getDimension();
    addItems();
    meshQual.reset(new JMeshQuality(mesh));
    meshQual->addAttribute("Quality");
}
///////////////////////////////////////////////////////////////////////////////

void JMeshGeometricQualityDialog :: addItems()
{
    name2id["Area"]  = JMeshQuality::AREA;
    name2id["AspectBeta"] = JMeshQuality::ASPECT_BETA;
    name2id["AspectGamma"] = JMeshQuality::ASPECT_GAMMA;
    name2id["AspectRatio"] = JMeshQuality::ASPECT_RATIO;
    name2id["ConditionNumber"] = JMeshQuality::CONDITION_NUMBER;
    name2id["DiagonalRatio"] = JMeshQuality::DIAGONAL_RATIO;
    name2id["Distortion"] = JMeshQuality::DISTORTION;
    name2id["EdgeLength"] = JMeshQuality::EDGE_LENGTH;
    name2id["Jacobian"] = JMeshQuality::JACOBIAN;
    name2id["MaxAngle"] = JMeshQuality::MAX_ANGLE;
    name2id["MinAngle"] = JMeshQuality::MIN_ANGLE;
    name2id["Oddy"] = JMeshQuality::ODDY;
    name2id["Relative Size Squared"] = JMeshQuality::RELATIVE_SIZE_SQUARED;
    name2id["Scaled Jacobian"] = JMeshQuality::SCALED_JACOBIAN;
    name2id["Shape"] = JMeshQuality::SHAPE;
    name2id["Shape and Size"] = JMeshQuality::SHAPE_AND_SIZE;
    name2id["Shear"] = JMeshQuality::SHEAR;
    name2id["Shear and Size"] = JMeshQuality::SHEAR_AND_SIZE;
    name2id["Skew"]    = JMeshQuality::SKEW;
    name2id["Stretch"] = JMeshQuality::STRETCH;
    name2id["Taper"]   = JMeshQuality::TAPER;
    name2id["Warpage"] = JMeshQuality::WARPAGE;
    name2id["Volume"]  = JMeshQuality::VOLUME;

    int nSize = comboBox->count();
    for( int i = 0; i < nSize; i++)
        comboBox->removeItem( comboBox->currentIndex() );

    int entityType = mesh->getTopology()->getElementsType(topDim);

    comboBox->addItem( QString::fromStdString("EdgeLength") );

    if( topDim == 2 ) {
        comboBox->addItem( QString::fromStdString("Area") );
        comboBox->addItem( QString::fromStdString("AspectRatio") );
        comboBox->addItem( QString::fromStdString("ConditionNumber") );
        comboBox->addItem( QString::fromStdString("Distortion") );
        comboBox->addItem( QString::fromStdString("MinAngle") );
        comboBox->addItem( QString::fromStdString("MaxAngle") );
        comboBox->addItem( QString::fromStdString("Shape") );
        comboBox->addItem( QString::fromStdString("Shape and Size") );
        comboBox->addItem( QString::fromStdString("Scaled Jacobian") );
        comboBox->addItem( QString::fromStdString("Relative Size Squared") );

        if( entityType == JFace::QUADRILATERAL) {
            comboBox->addItem( QString::fromStdString("Jacobian"));
            comboBox->addItem( QString::fromStdString("Oddy"));
            comboBox->addItem( QString::fromStdString("Shape and Size"));
            comboBox->addItem( QString::fromStdString("Shear"));
            comboBox->addItem( QString::fromStdString("Skew"));
            comboBox->addItem( QString::fromStdString("Stretch"));
            comboBox->addItem( QString::fromStdString("Taper"));
            comboBox->addItem( QString::fromStdString("Warpage"));
        }
    }

    if( topDim == 3 ) {
        comboBox->addItem( QString::fromStdString("ConditionNumber"));
        comboBox->addItem( QString::fromStdString("Distortion"));
        comboBox->addItem( QString::fromStdString("Jacobian"));
        comboBox->addItem( QString::fromStdString("Shape"));
        comboBox->addItem( QString::fromStdString("Shape and Size"));
        comboBox->addItem( QString::fromStdString("Relative Size Squared"));
        comboBox->addItem( QString::fromStdString("Scaled Jacobian"));
        comboBox->addItem( QString::fromStdString("Volume"));

        if( entityType == JCell::TETRAHEDRON ) {
            comboBox->addItem( QString::fromStdString("AspectBeta"));
            comboBox->addItem( QString::fromStdString("AspectGamma"));
        }

        if( entityType == JCell::HEXAHEDRON) {
            comboBox->addItem( QString::fromStdString("AspectRatio"));
            comboBox->addItem( QString::fromStdString("DiagonalRatio"));
            comboBox->addItem( QString::fromStdString("Oddy"));
            comboBox->addItem( QString::fromStdString("Shear"));
            comboBox->addItem( QString::fromStdString("Shear and Size"));
            comboBox->addItem( QString::fromStdString("Stretch"));
            comboBox->addItem( QString::fromStdString("Taper"));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeometricQualityDialog :: closeDialog()
{
    parentWidget()->show();
    this->close();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshGeometricQualityDialog :: setValues()
{
    if( mesh == nullptr ) return ;

    QString qstr = comboBox->currentText();
    string str = qstr.toUtf8().constData();

    if( name2id.find(str) == name2id.end() ) {
        cout << "Warning: No such quality measure for the entity" << str << endl;
        return;
    }

    int item = name2id[str];

    JWaitCursor waitCursor;
    waitCursor.start();

    int err;
    int measured_entity = -1;
    switch( topDim ) {
    case 2:
        measured_entity = 2;
        qData = meshQual->getFacesQuality( item, 0, 1);
        break;
    case 3:
        measured_entity = 3;
        qData = meshQual->getCellsQuality( item, 0, 1);
        break;
    }

    if( str == "EdgeLength" || str == "EdgeWeight") {
        qData = meshQual->getEdgesQuality( item, 0, 1);
        measured_entity = 1;
    }

    size_t nsize = qData.size();

    numSamplesLineEdit->setText( QString::number( nsize ) );

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

    minVal  = *boost::min_element( yData );
    maxVal  = *boost::max_element( yData );
    plotcurve->attach( nullptr );
    qhistogram->attach( nullptr );

    if( plotType == XYPLOT ) {
        qwtPlot->updateCanvasMargins();
//      qwtPlot->setAxisTitle( QwtPlot::xBottom, "Element");
        qwtPlot->setAxisTitle( QwtPlot::yLeft, str.c_str());
        plotcurve->setSymbol( new QwtSymbol( QwtSymbol::Ellipse, QColor(Qt::blue),
                                             QPen( Qt::red ), QSize( 5, 5 ) ) );
        plotcurve->setPen( QColor( Qt::black ) );
        plotcurve->setStyle( QwtPlotCurve::Dots );
        plotcurve->setRawSamples( &xData[0], &yData[0], nsize );
        qwtPlot->setAxisScale( QwtPlot::xBottom, 0, nsize);
        qwtPlot->setAxisScale( QwtPlot::yLeft, minVal, 1.1*maxVal);
        plotcurve->attach( qwtPlot );
    }

    qstr = numBinsLineEdit->text() ;
    int  numBins  = qstr.toInt();
    numBins = min( size_t(numBins), nsize );

    numBinsLineEdit->setText( QString::number(numBins) );

    if( plotType == HISTOGRAM ) {
        qwtPlot->setAxisTitle( QwtPlot::xBottom, "Bin number");
        qwtPlot->setAxisTitle( QwtPlot::yLeft, "#Elements");
        qhistogram->setValues( qData, numBins );
        qwtPlot->setAxisScale( QwtPlot::xBottom, 0, numBins);
        qwtPlot->setAxisAutoScale( QwtPlot::yLeft, 1);
        qhistogram->attach( qwtPlot );
    }
    numBinsLineEdit->setText( QString::number(numBins) );

    if( savedata ) {
        if( !filename.empty() ) {
            ofstream outfile( filename.c_str(), ios::out);
            for( size_t i = 0; i < qData.size(); i++)
                outfile << qData[i] << endl;
            outfile.close();
        }
    }

    if( statisticsDialog ) {
        statisticsDialog->setMesh(mesh);
        statisticsDialog->setEntityDimension(measured_entity);
        statisticsDialog->setQualityName(str);
        statisticsDialog->setValues();
    }
    qwtPlot->replot();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshGeometricQualityDialog :: setPlotType()
{
    if( histogramRadioButton->isChecked() ) plotType = HISTOGRAM;
    if( xyplotRadioButton->isChecked() )    plotType = XYPLOT;

    setValues();

}
///////////////////////////////////////////////////////////////////////////////
void JMeshGeometricQualityDialog :: saveAs()
{
    /*
         savedata = saveAsCheckBox->isChecked();
         QString qstr = filenameLineEdit->text();
         filename  = qstr.toUtf8().constData();
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshGeometricQualityDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshGeometricQualityDialog :: newQuality()
{
    numBinsLineEdit->setText( QString::number(100) );
    setValues();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshGeometricQualityDialog :: newHistogram()
{
    setValues();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshGeometricQualityDialog :: openStatisticsDialog()
{
    if( statisticsDialog == nullptr)
        statisticsDialog.reset(new JStatisticsDialog( this ));
    statisticsDialog->setViewManager( viewManager );
    statisticsDialog->show();
    setValues();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeometricQualityDialog :: makeConnections()
{
    PushButton( statisticsPushButton, [=] {openStatisticsDialog(); });
    PushButton( saveAsPushButton, [=] {saveAs();});
    ComboBox( comboBox, [=] {newQuality();});
    RadioButton( histogramRadioButton, [=] {setPlotType();});
    RadioButton( xyplotRadioButton, [=] {setPlotType();});
    LineEdit( numBinsLineEdit,  [=] {newHistogram();});
    PushButton( closePushButton, [=] {closeDialog();});

}
///////////////////////////////////////////////////////////////////////////////
