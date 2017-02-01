#include "QuadVerdictDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JQuadVerdictDialog :: JQuadVerdictDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

}

///////////////////////////////////////////////////////////////////////////////

JQuadVerdictDialog :: ~JQuadVerdictDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JQuadVerdictDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JQuadVerdictDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;
    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    size_t numfaces = mesh->getSize(2);
    size_t nCount = 0;
    for( size_t i = 0; i <  numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        if( f->getType() == JFace::QUADRILATERAL) nCount++;
    }

    numQuadsLineEdit->setText(QString::number( nCount) );
}

///////////////////////////////////////////////////////////////////////////////
void JQuadVerdictDialog :: setData( const vector<double> &quality, double minAccept, QLineEdit *minLineEdit,
                                    double maxAccept,  QLineEdit *maxLineEdit)
{
    if( quality.empty() ) return;

    QPalette acceptable, unacceptable;

    acceptable.setColor(QPalette::Base, QColor(0, 255, 0, 10));
    unacceptable.setColor(QPalette::Base, QColor(255, 0, 0, 10));

    double minVal = *boost::min_element(quality);

    if( minVal < minAccept)
        minLineEdit->setPalette(unacceptable);
    else
        minLineEdit->setPalette(acceptable);

    minLineEdit->setText( QString::number(minVal));


    double maxVal = *boost::max_element(quality);
    if( maxVal > maxAccept)
        maxLineEdit->setPalette( unacceptable);
    else
        maxLineEdit->setPalette( acceptable);

    maxLineEdit->setText( QString::number(maxVal));
}

///////////////////////////////////////////////////////////////////////////////

double JQuadVerdictDialog :: getAcceptable( const vector<double> &quality, double minval, double maxval)
{
    if( quality.empty() ) return 0.0;
    size_t nCount = 0;
    for( double val : quality )
        if( val >= minval && val <= maxval ) nCount++;
    return 100.0*nCount/(double)quality.size();
}
///////////////////////////////////////////////////////////////////////////////

void JQuadVerdictDialog :: getData()
{
    if( mesh == nullptr) return;
    JMeshQuality mq;
    mq.setMesh(mesh);

    double neginf = -0.99*std::numeric_limits<double>::infinity();
    double posinf =  0.99*std::numeric_limits<double>::infinity();
    int    pos    = 0;
    QString qstr = sampleComboBox->currentText();
    string  str  = StdString(qstr);
    if( str == "All") pos = 0;
    if( str == "Boundary") pos = 1;
    if( str == "Internal") pos = 2;

    vector<double>  quality;
    double goodval;

    //
    quality = mq.getFacesQuality( JMeshQuality::AREA, pos, 1);
    setData( quality, 0.0, minAreaLineEdit, posinf, maxAreaLineEdit);
    goodval = getAcceptable( quality, 0, posinf);
    acceptedAreaLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::ASPECT_RATIO, pos, 1);
    setData( quality, 1.0, minAspectRatioLineEdit, 1.3, maxAspectRatioLineEdit);
    goodval = getAcceptable( quality, 1.0, 1.3);
    acceptedAspectRatioLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::CONDITION_NUMBER, pos, 1);
    setData( quality, 1.0, minConditionLineEdit, 4.0, maxConditionLineEdit);
    goodval = getAcceptable( quality, 1.0, 4.0);
    acceptedConditionNumberLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::DISTORTION, pos, 1);
    setData( quality, 0.5, minDistortionLineEdit, 1.0, maxDistortionLineEdit);
    goodval = getAcceptable( quality, 0.5, 1.0);
    acceptedDistortionLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::JACOBIAN, pos, 1);
    setData( quality, 0.0, minJacobianLineEdit, posinf, maxJacobianLineEdit);
    goodval = getAcceptable( quality, 0.0, posinf);
    acceptedJacobianLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::MIN_ANGLE, pos, 1);
    setData( quality, 45.0, minminAngleLineEdit, 90.0, minmaxAngleLineEdit);
    goodval = getAcceptable( quality, 45.0, 90.0);
    acceptedMinAngleLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::MAX_ANGLE, pos, 1);
    setData( quality, 90.0, maxminAngleLineEdit, 135.0, maxmaxAngleLineEdit);
    goodval = getAcceptable( quality, 90.0, 135.0);
    acceptedMaxAngleLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::RELATIVE_SIZE_SQUARED, pos, 1);
    setData( quality, 0.3, minRelativeSizeLineEdit, 1.0, maxRelativeSizeLineEdit);
    goodval = getAcceptable( quality, 0.3, 1.0);
    acceptedRelativeSizeLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::SCALED_JACOBIAN, pos, 1);
    setData( quality, 0.5, minScaledJacobianLineEdit, 1.0, maxScaledJacobianLineEdit);
    goodval = getAcceptable( quality, 0.5, 1.0);
    acceptedScaledJacobianLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::SHAPE, pos, 1);
    setData( quality, 0.3, minShapeLineEdit, 1.0, maxShapeLineEdit);
    goodval = getAcceptable( quality, 0.3, 1.0);
    acceptedShapeLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::SHEAR, pos, 1);
    setData( quality, 0.3, minShearLineEdit, 1.0, maxShearLineEdit);
    goodval = getAcceptable( quality, 0.3, 1.0);
    acceptedShearLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::SHEAR_AND_SIZE, pos, 1);
    setData( quality, 0.2, minShearSizeLineEdit, 1.0, maxShearSizeLineEdit);
    goodval = getAcceptable( quality, 0.2, 1.0);
    acceptedShearSizeLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::SHAPE_AND_SIZE, pos, 1);
    setData( quality, 0.2, minShapeSizeLineEdit, 1.0, maxShapeSizeLineEdit);
    goodval = getAcceptable( quality, 0.2, 1.0);
    acceptedShapeSizeLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::SKEW, pos, 1);
    setData( quality, 0.2, minSkewLineEdit, 1.0, maxSkewLineEdit);
    goodval = getAcceptable( quality, 0.2, 1.0);
    acceptedSkewLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::STRETCH, pos, 1);
    setData( quality, 0.25, minStretchLineEdit, 1.0, maxStretchLineEdit);
    goodval = getAcceptable( quality, 0.25, 1.0);
    acceptedStretchLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::TAPER, pos, 1);
    setData( quality, 0.25, minTaperLineEdit, 1.0, maxTaperLineEdit);
    goodval = getAcceptable( quality, 0.25, 1.0);
    acceptedTaperLineEdit->setText( QString::number(goodval) );

    quality = mq.getFacesQuality( JMeshQuality::WARPAGE, pos, 1);
    setData( quality, 0.90, minWarpageLineEdit, 1.0, maxWarpageLineEdit);
    goodval = getAcceptable( quality, 0.90, 1.0);
    acceptedWarpageLineEdit->setText( QString::number(goodval) );

    int nCount = quality.size();
    numQuadsLineEdit->setText(QString::number(nCount) );

}
///////////////////////////////////////////////////////////////////////////////
void JQuadVerdictDialog :: makeConnections()
{
    PushButton( applyPushButton, [=] {getData(); });
    PushButton( closePushButton,  [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
