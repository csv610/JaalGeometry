#include "MeshFeaturesDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

JMeshFeaturesDialog :: JMeshFeaturesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;

    vecScaleLineEdit->setText( QString::number(1.0) );
}

///////////////////////////////////////////////////////////////////////////////

JMeshFeaturesDialog :: ~JMeshFeaturesDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: init()
{

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );

}
///////////////////////////////////////////////////////////////////////////////

void JMeshFeaturesDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;

    if( mesh == nullptr) return;
    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: angleDefect()
{
    if( mesh == nullptr ) return;

    double  angle  = defectAngleSlider->value();
    mesh->getGeometry()->setAngleDefects(angle);
//   Vertex::featureSet->addAttribute( "AngleDefect");
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: sharpEdges()
{
    if( mesh == nullptr ) return;

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));

    double  angle  = creaseAngleSlider->value();
    mesh->getGeometry()->setSharpEdges(angle);
//   meshViewer->displayAll(1,1);

    SharpEdgeColor edgeColor;
//   EdgeColor::assign( mesh, &edgeColor);

    size_t numedges = mesh->getSize(1);

    JEdgeRenderPtr eAttrib;
    size_t nCount = 0;
    for( size_t i = 0; i < numedges; i++) {
        JEdgePtr edge = mesh->getEdgeAt(i);
        if( edge->hasAttribute("CreaseAngle") ) {
            nCount++;
            edge->getAttribute("Render", eAttrib);
            eAttrib->display = 1;
        }
    }

    sharpEdgesLineEdit->setText( QString::number(nCount) );
    QApplication::restoreOverrideCursor();

//   Edge::featureSet->addAttribute( "SharpEdge");
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: setNodesColor( int id)
{
    size_t nsize =   mesh->getActiveSize(0);
    size_t numnodes = mesh->getSize(0);

    Eigen::VectorXd K;
    K.resize(nsize);
    int  index = 0;
    Array3D val;

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Curvature", val);
            if( id == 3)
                K(index++) = 0.5*(val[1]+ val[2]);
            else
                K(index++) = val[id];
        }
    }

    Eigen::MatrixXd Color;
    igl::jet(K, true, Color);

    index = 0;
    JNodeRenderPtr nAttrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", nAttrib);
            nAttrib->color[0] = Color(index,0);
            nAttrib->color[1] = Color(index,1);
            nAttrib->color[2] = Color(index,2);
            index++;
        }
    }
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->setSurfaceShade(JRender::SMOOTH_SHADE);
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: setGaussianCurvature()
{
    if( mesh == nullptr) return;

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));

    JMeshCurvature mcurve;
    mcurve.setMesh(mesh);
    int err = mcurve.setGaussianCurvature();
    if( !err) setNodesColor( 0 );

    QApplication::restoreOverrideCursor();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: setMeanCurvature()
{
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->setSurfaceShade(JRender::SMOOTH_SHADE);

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));

    JMeshCurvature mcurve;
    mcurve.setMesh(mesh);
    mcurve.setMeanCurvature();
    setNodesColor( 3 );

    QApplication::restoreOverrideCursor();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: getCurvatureDirections()
{
    if( mesh == nullptr) return;

    JMeshPtr minK, maxK;
    mesh->getAttribute("MinCurvatureDir", minK);
    if( minK ) meshViewer->removeObject(minK);

    mesh->getAttribute("MaxCurvatureDir", maxK);
    if( maxK ) meshViewer->removeObject(maxK);

    if( mesh->getSize(2) < 1) return;

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));

    JMeshCurvature mcurve;
    mcurve.setMesh(mesh);
    QString qstr = vecScaleLineEdit->text();
    mcurve.setScale( qstr.toDouble() );

    std::pair<JMeshPtr,JMeshPtr> minmax = mcurve.getCurvatureDirections();

    minK = minmax.first;
    maxK = minmax.second;
    meshViewer->addObject( minK );
    meshViewer->addObject( maxK );

    JColor redColor;
    redColor[0] = 0.9;
    redColor[1] = 0.1;
    redColor[2] = 0.1;
    redColor[3] = 1.0;
    JEdgeRenderPtr  eAttrib;

    if( minK ) {
        size_t numedges = minK->getSize(1);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = minK->getEdgeAt(i);
            if( edge->isActive() ) {
                edge->getAttribute("Render", eAttrib);
                eAttrib->color =  redColor;
            }
        }
        meshViewer->updateBuffers(minK);
    }

    JColor blueColor;
    blueColor[0] = 0.1;
    blueColor[1] = 0.1;
    blueColor[2] = 0.9;
    blueColor[3] = 1.0;

    if( maxK ) {
        size_t numedges = maxK->getSize(1);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = maxK->getEdgeAt(i);
            if( edge->isActive() ) {
                edge->getAttribute("Render", eAttrib);
                eAttrib->color =  blueColor;
            }
        }
        meshViewer->updateBuffers(maxK);
    }
    QApplication::restoreOverrideCursor();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: setEdgeAttributes( const JMeshPtr &kmesh)
{

    if( edgeAttribsDialog.get() == nullptr )
        edgeAttribsDialog.reset(new JEdgeAttributesDialog( this ));

    edgeAttribsDialog->setViewManager(viewManager);
    edgeAttribsDialog->setMesh(kmesh);

    JEdgeSequence eseq;
    size_t numedges = kmesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &e = kmesh->getEdgeAt(i);
        if( e->isActive() && !e->isBoundary()  ) eseq.push_back(e);
    }
    edgeAttribsDialog->setEdges(eseq);
    edgeAttribsDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshFeaturesDialog :: openMinKAttribDialog()
{
    if( mesh == nullptr) return;

    JMeshPtr minK;
    mesh->getAttribute("MinCurvatureVecField", minK);
    if( minK == nullptr) return;
    setEdgeAttributes( minK);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: openMaxKAttribDialog()
{
    if( mesh == nullptr) return;

    JMeshPtr maxK;
    mesh->getAttribute("MaxCurvatureVecField", maxK);
    if( maxK == nullptr) return;
    setEdgeAttributes( maxK );
}
///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: setNodesColor()
{
    QString qstr = nodesColorComboBox->currentText();
    string str  = qstr.toUtf8().constData();

    if( str == "Gaussian Curvature")
        setNodesColor(0);

    if( str == "Min Curvature")
        setNodesColor(1);

    if( str == "Max Curvature")
        setNodesColor(2);

    if( str == "Mean Curvature")
        setNodesColor(3);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshFeaturesDialog :: makeConnections()
{
    connect( defectAngleSlider, SIGNAL( valueChanged(int) ), this, SLOT( angleDefect() ));
    connect( creaseAngleSlider, SIGNAL( valueChanged(int) ), this, SLOT( sharpEdges() ));

    ComboBox( nodesColorComboBox, [=] {setNodesColor();});
    SpinBoxi( displayEigenVectorSpinBox, [=] {displayEigenVector(); });

    PushButton( minCurvaturePushButton,   [=] {openMinKAttribDialog();});
    PushButton( maxCurvaturePushButton,   [=] {openMaxKAttribDialog();});
    PushButton( gaussCurvaturePushButton,  [=] {setGaussianCurvature();});
    PushButton( meanCurvaturePushButton,   [=] {setMeanCurvature();});
    PushButton( curvatureDirectionsPushButton,  [=] {getCurvatureDirections();});

    PushButton( calculateEigenVectorsPushButton,  [=] {getEigenVectors();});

    PushButton( closePushButton,  [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////

void JMeshFeaturesDialog :: getEigenVectors()
{
    /*
        if( mesh == nullptr) return;

        if( meshSpectrum == nullptr) meshSpectrum.reset( new JMeshSpectrum);

        int  neigen = numEigenVectorsSpinBox->value();
        displayEigenVectorSpinBox->setMaximum( neigen-1);

        QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));
        meshSpectrum->setMesh(mesh);
        meshSpectrum->getEigenVectors(neigen);
        displayEigenVector();
        QApplication::restoreOverrideCursor();
    */

}
///////////////////////////////////////////////////////////////////////////////
void JMeshFeaturesDialog :: displayEigenVector()
{
    if( mesh == nullptr) return;
    /*
        int  neig = displayEigenVectorSpinBox->value();
        vector<double> eVals = meshSpectrum->getAbsEigenVector(neig);

        if( eVals.empty() ) return;

        vector<JColor> nodesColor;
        JColorMap::jet(eVals, nodesColor);
        JNodeRenderPtr nAttrib;
        size_t numnodes = mesh->getSize(0);
        int index = 0;
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->isActive() ) {
                vtx->getAttribute("Render", nAttrib);
                nAttrib->color = nodesColor[index++];
            }
        }
        meshViewer->getFaceDraw()->setRenderMode(JRenderMode::SMOOTH_SHADE);
    */

    meshViewer->updateBuffers(mesh);
}

