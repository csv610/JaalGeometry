#include "MeshTopologyQualityDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshTopologyQualityDialog :: JMeshTopologyQualityDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    model = new QStandardItemModel();
    tableView->setModel( model );

    makeConnections();
    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    lowDegreeSpinBox->setMinimum( 0 );
    equalDegreeSpinBox->setMinimum( 0 );
    highDegreeSpinBox->setMinimum( 0);

    qwtPlot->setCanvasBackground( Qt::white);
    grid.reset(new QwtPlotGrid());
    grid->setPen( QPen( Qt::black, 0, Qt::DotLine ) );
    grid->attach( qwtPlot );
    qhistogram.reset(new JHistogram( " Degree Distribution ", Qt::red ));
    plotcurve.reset( new QwtPlotCurve());

    magnifier.reset(new QwtPlotMagnifier(qwtPlot->canvas()));
    magnifier->setMouseButton(Qt::MidButton);

    panner.reset(new QwtPlotPanner(qwtPlot->canvas()));
    panner->setMouseButton(Qt::LeftButton);

    nodeColor.reset( new JNodeDegreeColor());
}

///////////////////////////////////////////////////////////////////////////////

JMeshTopologyQualityDialog :: ~JMeshTopologyQualityDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQualityDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQualityDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    mesh->buildRelations(0,0);

    setNodesDomain();

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str() ));
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQualityDialog :: closeDialog()
{
    this->close();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQualityDialog :: setNodesDomain()
{
    if( allNodesRadioButton->isChecked() ) nodepos = JMeshEntity::ANY_ENTITY;
    if( boundNodesRadioButton->isChecked() ) nodepos = JMeshEntity::BOUNDARY_ENTITY;
    if( internalNodesRadioButton->isChecked() ) nodepos = JMeshEntity::INTERNAL_ENTITY;
    getMeasure();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQualityDialog :: getMeasure()
{
    if( mesh == nullptr ) return;

    plotcurve->attach( nullptr );
    plotcurve->setSymbol( new QwtSymbol( QwtSymbol::Ellipse, QColor(Qt::yellow),
                                         QPen( Qt::red,2 ), QSize( 15, 15 ) ) );
    plotcurve->setPen( QPen( Qt::blue,5) );
    plotcurve->setStyle( QwtPlotCurve::Sticks );

    degreeCount = mesh->getTopology()->getNodesDegreeDistribution( nodepos );

    size_t nsize = degreeCount.size();
    if( nsize  == 0)  return;

    map<int,size_t>::const_iterator it;
    int xmax = 0;
    int numSamples = 0;
    for( it = degreeCount.begin(); it != degreeCount.end(); ++it) {
        xmax = max( xmax, it->first);
        numSamples += it->second;
    }
    maxDegree = xmax;


    lowDegreeSpinBox->setMaximum( maxDegree );
    equalDegreeSpinBox->setMaximum( maxDegree );
    highDegreeSpinBox->setMaximum( maxDegree );

    numSamplesLineEdit->setText( QString::number(numSamples) );

    xData.resize( xmax + 1 );
    yData.resize( xmax + 1 );

    for( int i = 0; i < xmax+1; i++) {
        xData[i] = i;
        yData[i] = 0;
    }

    for( it = degreeCount.begin(); it != degreeCount.end(); ++it) {
        yData[it->first] = it->second;
    }

    plotcurve->setRawSamples( &xData[0], &yData[0], xData.size() );

    size_t numNodes = mesh->getSize(0);

    model->clear();

    model->setHorizontalHeaderItem(0, new QStandardItem(QString("Degree")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("#Elements")));
    model->setHorizontalHeaderItem(2, new QStandardItem(QString("Percentage")));
    for( it = degreeCount.begin(); it != degreeCount.end(); ++it) {
        int ndegree = it->first;
        if( ndegree ) {
            QList<QStandardItem*> newRow;
            QStandardItem* item1 = new QStandardItem(QString("%0").arg(ndegree));
            newRow.append(item1);

            int nsize = it->second;
            QStandardItem* item2 = new QStandardItem(QString("%1").arg(nsize));
            newRow.append(item2);

            double perc = nsize*100/(double)numNodes;
            QStandardItem* item3 = new QStandardItem(QString("%2").arg(perc));
            newRow.append(item3);
            model->appendRow(newRow);
        }
    }

    plotcurve->attach( qwtPlot );
    qwtPlot->replot();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQualityDialog :: displayDegrees()
{
    if( mesh == nullptr ) return;

    QString str;
    int nd;

    if( lowDegreeCheckBox->isChecked() ) {
        nd  = lowDegreeSpinBox->value();
        nodeColor->setLowDegree(nd);
    } else
        nodeColor->setLowDegree(0);

    if( equalDegreeCheckBox->isChecked() ) {
        nd  = equalDegreeSpinBox->value();
        nodeColor->setEqualDegree(nd);
    } else
        nodeColor->setEqualDegree(0);

    if( highDegreeCheckBox->isChecked() ) {
        nd  = highDegreeSpinBox->value();
    } else
        nodeColor->setHighDegree(maxDegree);

    nodeColor->setMesh(mesh);

    meshViewer->updateBuffers(mesh);
}

/////////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQualityDialog :: makeConnections()
{
    RadioButton( allNodesRadioButton, [=] {setNodesDomain();});
    RadioButton( boundNodesRadioButton, [=] {setNodesDomain(); });
    RadioButton( internalNodesRadioButton, [=] {setNodesDomain(); });

    CheckBox( lowDegreeCheckBox, [=] {displayDegrees();});
    CheckBox( equalDegreeCheckBox, [=] {displayDegrees(); });
    CheckBox( highDegreeCheckBox, [=] {displayDegrees();});

    PushButton( closePushButton, [=] {closeDialog();});
}
///////////////////////////////////////////////////////////////////////////////
