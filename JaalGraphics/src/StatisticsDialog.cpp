#include "StatisticsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog:: display( const JEdgePtr &edge)
{
    if( !edge->isActive() ) return;

    JEdgeRenderPtr eAttrib;
    edge->getAttribute("Render", eAttrib);

    double val;
    edge->getAttribute("Quality", val);

    if( displayLower && val < lowerCutoff) {
        eAttrib->display = 1;
    }

    if( displayUpper && val > upperCutoff) {
        eAttrib->display = 1;
    }

    if( displayMiddle && val >= lowerCutoff && val <= upperCutoff) {
        eAttrib->display = 1;
    }
}


////////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog:: display_edges_quality()
{
    if( meshViewer == nullptr ) return;
    if( mesh == nullptr ) return;

    size_t numedges = mesh->getSize(1);
    #pragma omp parallel for
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge= mesh->getEdgeAt(i);
        display(edge);
    }
    meshViewer->updateBuffers(mesh);

}

///////////////////////////////////////////////////////////////////////////////
void JStatisticsDialog:: display( const JFacePtr &face)
{
    if( !face->isActive() ) return;
    double val;

    JFaceRenderPtr fAttrib;
    JEdgeRenderPtr eAttrib;

    int err = face->getAttribute("Quality",val);
    if( err ) {
        cout << "Warning: Face does not have quality attribute" << endl;
        return;
    }

    face->getAttribute("Render", fAttrib);
    fAttrib->display = 0;

    if( displayLower && val < lowerCutoff) {
        fAttrib->display = 1;
        fAttrib->color   = redColor;
        for( int j = 0; j < face->getSize(1); j++) {
            const JEdgePtr &edge = face->getEdgeAt(j);
            edge->getAttribute("Render", eAttrib);
            eAttrib->display = 1;
        }
    }

    if( displayUpper && val > upperCutoff) {
        fAttrib->display = 1;
        fAttrib->color   = blueColor;
        for( int j = 0; j < face->getSize(1); j++) {
            const JEdgePtr &edge = face->getEdgeAt(j);
            edge->getAttribute("Render", eAttrib);
            eAttrib->display = 1;
        }
    }

    if( displayMiddle && val >= lowerCutoff && val <= upperCutoff) {
        fAttrib->display = 1;
        fAttrib->color   = greenColor;
        for( int j = 0; j < face->getSize(1); j++) {
            const JEdgePtr &edge = face->getEdgeAt(j);
            edge->getAttribute("Render", eAttrib);
            eAttrib->display = 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog:: display_faces_quality()
{
    if( mesh == nullptr ) return;

    size_t numfaces = mesh->getSize(2);

//  #pragma omp parallel for
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face =  mesh->getFaceAt(i);
        display(face);
    }

    meshViewer->updateBuffers(mesh);

}

///////////////////////////////////////////////////////////////////////////////
void JStatisticsDialog::display( const JCellPtr &cell)
{
    if( !cell->isActive() ) return;

    double val;
    int err = cell->getAttribute("Quality",val);
    if( err ) return;

    JCellRenderPtr cAttrib;
    cell->getAttribute("Render", cAttrib);
    cAttrib->display = 0;

    JFaceRenderPtr fAttrib;
    bool display_edges = 0;
    if( displayLower && val < lowerCutoff) {
        cAttrib->display = 1;
        for( int j = 0; j < cell->getSize(2); j++) {
            const JFacePtr &face = cell->getFaceAt(j);
            face->getAttribute("Render", fAttrib);
            fAttrib->display = 1;
        }
        display_edges = 1;
    }

    if( displayUpper && val > upperCutoff) {
        cAttrib->display = 1;
        for( int j = 0; j < cell->getSize(2); j++) {
            const JFacePtr &face = cell->getFaceAt(j);
            face->getAttribute("Render", fAttrib);
            fAttrib->display = 1;
        }
        display_edges = 1;
    }

    if( displayMiddle && val >= lowerCutoff && val <= upperCutoff) {
        cAttrib->display = 1;
        for( int j = 0; j < cell->getSize(2); j++) {
            const JFacePtr &face = cell->getFaceAt(j);
            face->getAttribute("Render", fAttrib);
            fAttrib->display = 1;
        }
        display_edges = 1;
    }

    JEdgeRenderPtr eAttrib;
    if( display_edges ) {
        for( int j = 0; j < cell->getSize(1); j++) {
            const JEdgePtr &edge = cell->getEdgeAt(j);
            edge->getAttribute("Render", eAttrib);
            eAttrib->display = 1;
            eAttrib->scale   = 2;
        }
    }

}

/////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog::display_cells_quality()
{
    if( mesh == nullptr ) return;

    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);
    #pragma omp parallel for private(fAttrib)
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        f->getAttribute("Render", fAttrib);
        fAttrib->display = 0;
    }

    JEdgeRenderPtr eAttrib;
    size_t numedges = mesh->getSize(1);
    #pragma omp parallel for private(eAttrib)
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &e = mesh->getEdgeAt(i);
        e->getAttribute("Render", eAttrib);
        eAttrib->display = 0;
        if( e->isBoundary()) eAttrib->display = 1;
    }

    size_t numcells = mesh->getSize(3);

    #pragma omp parallel for
    for( size_t i = 0; i < numcells; i++) {
        const JCellPtr &cell =  mesh->getCellAt(i);
        display(cell);
    }

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

JStatisticsDialog :: JStatisticsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    redColor   = JEntityColor::getColor("Pink");
    greenColor = JEntityColor::getColor("LightGreen");
    blueColor  = JEntityColor::getColor("LightBlue");
}

///////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog :: getCount()
{
    size_t nrange = 0;
    size_t nabove = 0;
    size_t nbelow = 0;

    QString qstr;
    qstr = lowerCutoffLineEdit->text();
    lowerCutoff = qstr.toDouble();

    qstr = upperCutoffLineEdit->text();
    upperCutoff = qstr.toDouble();

    for( size_t i = 0; i < data.size(); i++) {
        if( data[i] < lowerCutoff ) nbelow ++;
        if( data[i] > upperCutoff ) nabove ++;
    }

    nrange = data.size() - nabove - nbelow;

    numSamplesLineEdit->setText( QString::number( data.size() ) );
    numBelowCutoffLineEdit->setText( QString::number( nbelow ) );
    numAboveCutoffLineEdit->setText( QString::number( nabove ) );
    numInRangeLineEdit->setText( QString::number( nrange ) );

    displayLower  = displayLowerCheckBox->isChecked();
    displayUpper  = displayUpperCheckBox->isChecked();
    displayMiddle = displayMiddleCheckBox->isChecked();

    switch(entityDim)
    {
    case 2:
        display_faces_quality();
        break;
    case 3:
        display_cells_quality();
        break;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    data.clear();
}

///////////////////////////////////////////////////////////////////////////////
void JStatisticsDialog :: setQualityName( const string &s)
{
    measureLineEdit->setText( QString(s.c_str()));
}
///////////////////////////////////////////////////////////////////////////////
void JStatisticsDialog :: setValues()
{
    data.clear();
    if( mesh == nullptr) return;

    double q;
    size_t  numedges = mesh->getSize(1);
    size_t  numfaces = mesh->getSize(2);
    size_t  numcells = mesh->getSize(3);

    switch( entityDim ) {
    case 1:
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                int err = edge->getAttribute("Quality", q);
                if( !err) data.push_back(q);
            }
        }
        break;
    case 2:
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                int err = face->getAttribute("Quality", q);
                if( !err) data.push_back(q);
            }
        }
        break;
    case 3:
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                int err = cell->getAttribute("Quality", q);
                if( !err) data.push_back(q);
            }
        }
        break;
    }

    if( data.empty() ) return;

    minVal  = *boost::min_element( data );
    maxVal  = *boost::max_element( data );

    double sumval  =  std::accumulate(data.begin(), data.end(), 0.0);
    double avgval  =  JMath::average_value(data);
    double meanval =  JMath::mean_value(data);
    double stddev  =  JMath::standard_deviation(data);

    lowerCutoff = minVal;
    upperCutoff = maxVal;

    numSamplesLineEdit->setText( QString::number(data.size()) );
    sumLineEdit->setText( QString::number(sumval) );
    avgQualityLineEdit->setText( QString::number(avgval) );
    meanQualityLineEdit->setText( QString::number(meanval) );
    stdevQualityLineEdit->setText( QString::number(stddev) );

    minQualityLineEdit->setText( QString::number(minVal) );
    maxQualityLineEdit->setText( QString::number(maxVal) );
    relativeMinLineEdit->setText( QString::number(minVal/meanval) );
    relativeMaxLineEdit->setText( QString::number(maxVal/meanval) );


    lowerCutoffSlider->setSliderPosition( 0 );
    upperCutoffSlider->setSliderPosition( upperCutoffSlider->maximum() );

    lowerCutoffLineEdit->setText( QString::number( minVal ) );
    upperCutoffLineEdit->setText( QString::number( maxVal ) );

    getCount();
}

///////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog :: resetLowerCutoff()
{
    lowerCutoffSlider->setSliderPosition( 0 );
    lowerCutoff = minVal;
    getCount();
}
///////////////////////////////////////////////////////////////////////////////


void JStatisticsDialog :: resetUpperCutoff()
{
    upperCutoffSlider->setSliderPosition( upperCutoffSlider->maximum() );
    upperCutoff = maxVal;
    getCount();
}

///////////////////////////////////////////////////////////////////////////////
/*
void JStatisticsDialog :: getLowerEntities()
{
    int intval = lowerCutoffSlider->value();
    lowerCutoff = minVal + 0.010*(intval)*(maxVal-minVal);
    lowerCutoffLineEdit->setText( QString::number( lowerCutoff  ) );
    getCount();
}

///////////////////////////////////////////////////////////////////////////////
void JStatisticsDialog :: getUpperEntities()
{
    int intval = upperCutoffSlider->value();
    upperCutoff = minVal + 0.01*(intval)*(maxVal-minVal);
    upperCutoffLineEdit->setText( QString::number( upperCutoff  ) );
    getCount();
}
*/
///////////////////////////////////////////////////////////////////////////////
void JStatisticsDialog :: setSliderValue()
{
    int intval;
    intval = lowerCutoffSlider->value();
    lowerCutoff = minVal + 0.0010*(intval)*(maxVal-minVal);
    lowerCutoffLineEdit->setText( QString::number( lowerCutoff  ) );

    intval = upperCutoffSlider->value();
    upperCutoff = minVal + 0.001*(intval)*(maxVal-minVal);
    upperCutoffLineEdit->setText( QString::number( upperCutoff  ) );
}
///////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog :: editLowerCutoff()
{
    QString str = lowerCutoffLineEdit->text() ;
    double val  = str.toDouble();
    if( val < minVal ) val = minVal;
    if( val > maxVal ) val = maxVal;
    int pos = (val-minVal)*1000.0/(maxVal-minVal);
    lowerCutoffSlider->setValue(pos);

    lowerCutoff = val;

    getCount();
}

///////////////////////////////////////////////////////////////////////////////
void JStatisticsDialog :: editUpperCutoff()
{
    QString str = upperCutoffLineEdit->text() ;
    double val  = str.toDouble();
    if( val < minVal ) val = minVal;
    if( val > maxVal ) val = maxVal;
    int pos = (val-minVal)*1000.0/(maxVal-minVal);
    upperCutoffSlider->setValue(pos);

    upperCutoff = val;

    getCount();
}
///////////////////////////////////////////////////////////////////////////////

void JStatisticsDialog :: makeConnections()
{
    PushButton( closePushButton,  [=] {close(); });
    PushButton( lowerCutoffPushButton, [=] {resetLowerCutoff();});
    PushButton( upperCutoffPushButton, [=] {resetUpperCutoff();});

    connect( lowerCutoffSlider, SIGNAL( sliderReleased() ),  this, SLOT( getCount() ));
    connect( upperCutoffSlider, SIGNAL( sliderReleased() ),  this, SLOT( getCount() ));
    connect( lowerCutoffSlider, SIGNAL( valueChanged(int) ),  this, SLOT( setSliderValue() ));
    connect( upperCutoffSlider, SIGNAL( valueChanged(int) ),  this, SLOT( setSliderValue() ));

    CheckBox( displayLowerCheckBox,  [=] { getCount();});
    CheckBox( displayUpperCheckBox,  [=] { getCount();});
    CheckBox( displayMiddleCheckBox, [=] { getCount();});

    LineEdit( lowerCutoffLineEdit,  [=] {editLowerCutoff();});
    LineEdit( upperCutoffLineEdit,  [=] {editUpperCutoff();});
}

///////////////////////////////////////////////////////////////////////////////
