#include "MeshMapQualityDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshMapQualityDialog :: JMeshMapQualityDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshMapQualityDialog :: ~JMeshMapQualityDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshMapQualityDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");

    qwtPlot->setCanvasBackground( Qt::white);
    grid.reset(new QwtPlotGrid());
    grid->setPen( QPen( Qt::black, 0, Qt::DotLine ) );
    grid->attach( qwtPlot );

    magnifier.reset(new QwtPlotMagnifier(qwtPlot->canvas()));
    magnifier->setMouseButton(Qt::MidButton);

    panner.reset(new QwtPlotPanner(qwtPlot->canvas()));
    panner->setMouseButton(Qt::LeftButton);

    curve1.reset( new QwtPlotCurve());
    curve2.reset( new QwtPlotCurve());
}

///////////////////////////////////////////////////////////////////////////////

void JMeshMapQualityDialog :: setSourceMesh( const JMeshPtr &m)
{
    sourceMesh = m;
    if( sourceMesh == nullptr ) return ;

    string name = sourceMesh->getName();
    sourceNameLineEdit->setText(QString(name.c_str()));

    if( targetMesh)  checkQuality();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshMapQualityDialog :: setTargetMesh( const JMeshPtr &m)
{
    targetMesh = m;
    if( targetMesh == nullptr ) return ;

    string name = targetMesh->getName();
    targetNameLineEdit->setText(QString(name.c_str()));

    if( sourceMesh ) checkQuality();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshMapQualityDialog :: setValues( const vector<double> &aval, const vector<double> &bval)
{
    size_t nsize = aval.size();

    xData.resize(nsize);
    for( size_t i = 0; i < nsize; i++)
        xData[i] = 1.0*i;

    double minVal;
    minVal  = *boost::min_element(aval);
    minVal  = min(minVal, *boost::min_element( bval));

    double maxVal;
    maxVal  = *boost::max_element( aval );
    maxVal  = max(maxVal, *boost::max_element( bval ));

    qwtPlot->setAxisScale( QwtPlot::xBottom, 0, nsize);
    qwtPlot->setAxisScale( QwtPlot::yLeft, minVal, 1.1*maxVal);
    qwtPlot->updateCanvasMargins();

    curve1->attach( nullptr );
    curve1->setSymbol( new QwtSymbol( QwtSymbol::Triangle, QColor(Qt::blue),
                                      QPen( Qt::red ), QSize( 5, 5 ) ) );
    curve1->setPen( QColor( Qt::red ) );
    curve1->setStyle( QwtPlotCurve::Lines );
    curve1->setRawSamples( &xData[0], &aval[0], nsize );
    curve1->attach( qwtPlot );

    curve2->attach( nullptr );
    curve2->setSymbol( new QwtSymbol( QwtSymbol::Rect, QColor(Qt::blue),
                                      QPen( Qt::red ), QSize( 5, 5 ) ) );
    curve2->setPen( QColor( Qt::blue ) );
    curve2->setStyle( QwtPlotCurve::Lines );
    curve2->setRawSamples( &xData[0], &bval[0], nsize );
    curve2->attach( qwtPlot );

    qwtPlot->replot();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshMapQualityDialog :: checkQuality()
{
    QMessageBox msg;
    if( sourceMesh == nullptr || targetMesh == nullptr) {
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning: A source and target mesh must be specified");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    if( sourceMesh->getSize(0) != targetMesh->getSize(0)) {
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning: A source and target mesh must have same number of nodes");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    if( sourceMesh->getSize(2) != targetMesh->getSize(2)) {
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning: A source and target mesh must have same number of faces");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    QString qstr = qualityComboBox->currentText();

    JEdgeSequence edges;
    JFaceSequence faces;
    if( qstr == QString("Length") ) {
        sourceMesh->enumerate(1);
        targetMesh->enumerate(1);
        size_t numedges = sourceMesh->getSize(1);
        edges.resize(numedges);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = sourceMesh->getEdgeAt(i);
            double l1 = JEdgeGeometry::getLength(edge);
            edge->setAttribute("Length", l1);
            edges[i] = edge;
        }

        std::sort( edges.begin(), edges.end(), []( const JEdgePtr &edge1, const JEdgePtr &edge2)
        {
            double l1, l2;
            edge1->getAttribute("Length", l1);
            edge2->getAttribute("Length", l2);
            return l1 < l2;
        });
        aData.resize(numedges);
        bData.resize(numedges);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &srcEdge = edges[i];
            aData[i] = JEdgeGeometry::getLength(srcEdge);
            const JEdgePtr &dstEdge = targetMesh->getEdgeAt(srcEdge->getID());
            bData[i] = JEdgeGeometry::getLength(dstEdge);
        }
        setValues(aData, bData);
    }

    if( qstr == QString("MinAngle") ) {
        sourceMesh->enumerate(2);
        targetMesh->enumerate(2);
        double angle;
        vector<double> angles;
        size_t numfaces = sourceMesh->getSize(2);
        faces.resize(numfaces);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = sourceMesh->getFaceAt(i);
            JFaceGeometry::getAngles(face, angles, ANGLE_IN_DEGREES);
            angle = *boost::min_element( angles);
            face->setAttribute("Angle", angle);
            faces[i] = face;
        }

        std::sort( faces.begin(), faces.end(), []( const JFacePtr &face1, const JFacePtr &face2)
        {
            double val1, val2;
            face1->getAttribute("Angle", val1);
            face2->getAttribute("Angle", val2);
            return val1 < val2;
        });

        aData.resize(numfaces);
        bData.resize(numfaces);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &srcFace = faces[i];
            srcFace->getAttribute("Angle", angle);
            aData[i] = angle;

            const JFacePtr &dstFace = targetMesh->getFaceAt(srcFace->getID());
            JFaceGeometry::getAngles(dstFace, angles, ANGLE_IN_DEGREES);
            bData[i] = *boost::min_element( angles );
        }
        setValues(aData, bData);
    }

    if( qstr == QString("MaxAngle") ) {
        sourceMesh->enumerate(2);
        targetMesh->enumerate(2);
        double angle;
        vector<double> angles;
        size_t numfaces = sourceMesh->getSize(2);
        faces.resize(numfaces);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = sourceMesh->getFaceAt(i);
            JFaceGeometry::getAngles(face, angles, ANGLE_IN_DEGREES);
            angle = *boost::max_element( angles );
            face->setAttribute("Angle", angle);
            faces[i] = face;
        }

        std::sort( faces.begin(), faces.end(), []( const JFacePtr &face1, const JFacePtr &face2)
        {
            double val1, val2;
            face1->getAttribute("Angle", val1);
            face2->getAttribute("Angle", val2);
            return val1 < val2;
        });
        aData.resize(numfaces);
        bData.resize(numfaces);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &srcFace = faces[i];
            srcFace->getAttribute("Angle", angle);
            aData[i] = angle;

            const JFacePtr &dstFace = targetMesh->getFaceAt(srcFace->getID());
            JFaceGeometry::getAngles(dstFace, angles, ANGLE_IN_DEGREES);
            bData[i] = *boost::max_element( angles );
        }
        setValues(aData, bData);
    }

    if( qstr == QString("Area") ) {
        sourceMesh->enumerate(2);
        targetMesh->enumerate(2);
        double area;
        size_t numfaces = sourceMesh->getSize(2);
        faces.resize(numfaces);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = sourceMesh->getFaceAt(i);
            area = JFaceGeometry::getArea(face);
            face->setAttribute("Area", area);
            faces[i] = face;
        }

        std::sort( faces.begin(), faces.end(), []( const JFacePtr &face1, const JFacePtr &face2)
        {
            double val1, val2;
            face1->getAttribute("Area", val1);
            face2->getAttribute("Area", val2);
            return val1 < val2;
        });
        aData.resize(numfaces);
        bData.resize(numfaces);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &srcFace = faces[i];
            srcFace->getAttribute("Area", area);
            aData[i] = area;

            const JFacePtr &dstFace = targetMesh->getFaceAt(srcFace->getID());
            bData[i] = JFaceGeometry::getArea(dstFace);
        }
        setValues(aData, bData);
    }
}
///////////////////////////////////////////////////////////////////////////////

void JMeshMapQualityDialog :: makeConnections()
{
    ComboBox( qualityComboBox, [=] { checkQuality(); });
    PushButton( closePushButton, [=] { close(); });
}

///////////////////////////////////////////////////////////////////////////////
