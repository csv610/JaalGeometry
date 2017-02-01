#include "MeshUntangleDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshUntangleDialog :: JMeshUntangleDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;
    highlightColor = JEntityColor::getColor("Red");
}

///////////////////////////////////////////////////////////////////////////////

JMeshUntangleDialog :: ~JMeshUntangleDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshUntangleDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshUntangleDialog :: setMesh( const JMeshPtr &m)
{
    srcMesh = m;
    if( srcMesh == nullptr ) return ;
    initialized = 0;
    inflatedMesh.reset();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshUntangleDialog :: initMesh()
{
    if( initialized ) return;

    srcMesh->getGeometry()->getCoordsArray(orgCoords,l2g);

    string name = srcMesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    JNodeRenderPtr vAttrib;
    size_t numnodes =  srcMesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = srcMesh->getNodeAt(i);
        if( v->isActive() ) {
            v->getAttribute("Render", vAttrib);
            prevNodeColor[v] = vAttrib->color;
        }
    }

    JEdgeRenderPtr eAttrib;
    size_t numedges =  srcMesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &e = srcMesh->getEdgeAt(i);
        if( e->isActive() ) {
            e->getAttribute("Render", eAttrib);
            prevEdgeColor[e] = eAttrib->color;
        }
    }

    JFaceRenderPtr fAttrib;
    size_t numfaces =  srcMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = srcMesh->getFaceAt(i);
        if( f->isActive() ) {
            f->getAttribute("Render", fAttrib);
            prevFaceColor[f] = fAttrib->color;
        }
    }

    untangle.setMesh(srcMesh);
    int numInverted = untangle.countInverted(srcMesh);
    if( numInverted  == 0) return;

    inflatedMesh = untangle.getInflatedMesh();
    meshViewer->addObject(inflatedMesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshUntangleDialog :: setColor()
{
    JMeshRenderPtr mrender;
    if( srcMeshRadioButton) {
        displayInverted(srcMesh);
        inflatedMesh->getAttribute("Render", mrender);
        mrender->display = 0;
    } else {
        displayInverted(inflatedMesh);
        srcMesh->getAttribute("Render", mrender);
        mrender->display = 0;
    }
}
///////////////////////////////////////////////////////////////////////////////
/*
void JMeshUntangleDialog :: setColor()
{
    QColor color = QColorDialog::getColor();
    if( meshViewer == nullptr ) return;

    float rgb[3];
    highlightColor[0] = color.red()/255.0;
    highlightColor[1] = color.green()/255.0;
    highlightColor[2] = color.blue()/255.0;
}
*/
///////////////////////////////////////////////////////////////////////////////

void JMeshUntangleDialog :: setLaplaceSmoothing()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshUntangleDialog :: backProjection()
{
    int iter = numBackProjectionsSpinBox->value();
    untangle.setMaxBackProjectionSteps(iter);
    untangle.startBackProjection();

    double  achanged =  untangle.getAreaChanged();
    areaVolChangeLineEdit->setText(QString::number(achanged));
    displayOffset();

    meshViewer->updateBuffers(srcMesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshUntangleDialog :: displayInverted( const JMeshPtr &mesh)
{
    if( meshViewer == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->display = 1;

    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str()));

    meshViewer->getFaceDraw()->setFrontFace(GL_CCW);
    meshViewer->getFaceDraw()->setCulling( GL_BACK );

    /*
        if( displayInvertedCheckBox->isChecked() ) {
            meshViewer->getFaceDraw()->setCulling( GL_FRONT );
            meshViewer->getFaceDraw()->setFrontFace(GL_CW);
        }
    */

    int numInv = mesh->getGeometry()->getNumOfInvertedElements();

    JColor red  = JEntityColor::getColor("Red");
    JColor black = JEntityColor::getColor("Black");

    int entityDim = mesh->getTopology()->getDimension();

    if( entityDim == 2) {
        meshViewer->getFaceDraw()->setLights(0);

        size_t numnodes = mesh->getSize(0);
        JNodeRenderPtr vAttrib;
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() ) {
                v->getAttribute("Render", vAttrib);
                vAttrib->color   = black;
                vAttrib->glyph   = 0;
                vAttrib->display = 1;
            }
        }

        size_t numedges = mesh->getSize(1);
        JEdgeRenderPtr eAttrib;
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &e = mesh->getEdgeAt(i);
            if( e->isActive() && !e->isBoundary() )  {
                e->getAttribute("Render", eAttrib);
                eAttrib->display = 0;
            }
        }

        size_t numfaces = mesh->getSize(2);
        JFaceRenderPtr fAttrib;
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() ) {
                f->getAttribute("Render", fAttrib);
                fAttrib->color = highlightColor;
            }
        }

        if( mesh->getTopology()->getElementsType(2) == JFace::TRIANGLE)  {
            numInvertedTrisLineEdit->setText(QString::number(numInv));
        }

        if( mesh->getTopology()->getElementsType(2) == JFace::QUADRILATERAL) {
            numInvertedQuadsLineEdit->setText(QString::number(numInv));
        }
        meshViewer->updateBuffers(mesh);
    }

    if( entityDim == 3) {
        size_t numnodes = mesh->getSize(0);
        JNodeRenderPtr vAttrib;
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() ) {
                v->getAttribute("Render", vAttrib);
                vAttrib->color   = black;
                vAttrib->display = 1;
            }
        }

        size_t numedges = mesh->getSize(1);
        JEdgeRenderPtr eAttrib;
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &e = mesh->getEdgeAt(i);
            if( e->isActive() )  {
                e->getAttribute("Render", eAttrib);
                eAttrib->display = 0;
            }
        }

        size_t numfaces = mesh->getSize(2);
        JFaceRenderPtr fAttrib;
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() ) {
                f->getAttribute("Render", fAttrib);
                fAttrib->display = 0;
            }
        }

        size_t numcells = mesh->getSize(3);
        JCellRenderPtr cAttrib;
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &c = mesh->getCellAt(i);
            if( c->isActive() ) {
                c->getAttribute("Render", cAttrib);
                cAttrib->color   = red;
                cAttrib->display = 0;
            }
        }

        JCellSequence inverted = untangle.getInvertedCells();
        for( const JCellPtr &cell : inverted ) {
            cell->getAttribute("Render", cAttrib);
            cAttrib->color   = highlightColor;
            cAttrib->display = 1;
        }

        if( mesh->getTopology()->getElementsType(3) == JCell::TETRAHEDRON)
            numInvertedTetsLineEdit->setText(QString::number(numInv));

        if( mesh->getTopology()->getElementsType(3) == JCell::HEXAHEDRON)
            numInvertedHexsLineEdit->setText(QString::number(numInv));

        meshViewer->getFaceDraw()->setCulling(0);
        meshViewer->updateBuffers(mesh);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshUntangleDialog :: displayOffset()
{
    vector<double> offset = untangle.getOffset();
    int numnodes = offset.size();

    double maxval = *boost::max_element(offset );
    double minval = *boost::min_element(offset );

    maxOffsetLineEdit->setText(QString::number(maxval));

    Eigen::VectorXd K;
    K.resize(numnodes);
    for( size_t i = 0; i < numnodes; i++) {
        int val = 256*(offset[i]/(maxval-minval));
        K(i) = val;
    }

    Eigen::MatrixXd Color;
    igl::jet(K, true, Color);

    JNodeRenderPtr nAttrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = inflatedMesh->getNodeAt(i);
        vtx->getAttribute("Render", nAttrib);
        nAttrib->color[0] = Color(i,0);
        nAttrib->color[1] = Color(i,1);
        nAttrib->color[2] = Color(i,2);
    }

    meshViewer->updateBuffers(inflatedMesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshUntangleDialog :: execute()
{
    JStopWatch stopwatch;
    stopwatch.start();

    int numInverted = untangle.countInverted(srcMesh);
    if( numInverted == 0) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Input Mesh is already untangled");
        msg.setStandardButtons(QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok) return;
    }

    if( method == JAAL_METHOD)  useJaalMethod();
    if( method == SACHT_METHOD) useSachtMethod();

    stopwatch.stop();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshUntangleDialog :: useSachtMethod()
{

}
///////////////////////////////////////////////////////////////////////////////
void JMeshUntangleDialog :: useJaalMethod()
{
    if( inflatedMesh == nullptr) return;

    // Now show only the inflated model and hide the original mesh...
    JMeshRenderPtr mrender;
    srcMesh->getAttribute("Render", mrender);
    mrender->display = 0;

    displayInverted(inflatedMesh);

    QString qstr = projectEnergyComboBox->currentText();
    string  str  = StdString(qstr);
    int     eid  = JLocallyInjectiveMap::getEnergyType(str);

    int iter = numBoundaryModificationsSpinBox->value();
    untangle.setMaxInflationSteps(iter);

    QMessageBox msg;
    msg.setIcon(QMessageBox::Warning);
    msg.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel);

    int numInverted = 0;
    while(1) {
        untangle.startInflation();
        numInverted = untangle.countInverted(inflatedMesh);
        if( numInverted ) {
            displayInverted(inflatedMesh);
            msg.setText("Info: Mesh inflation not recovered inverted elements: Perform more ? ");
            int ret = msg.exec();
            if( ret == QMessageBox::Cancel) return;
        }
    }
    if( numInverted ) return;

    untangle.setEnergyType(eid);

    iter = numBackProjectionsSpinBox->value();
    untangle.setMaxBackProjectionSteps(iter);

    while(1) {
        untangle.startBackProjection();
        double maxDist = untangle.getMaxDistance();
        if( maxDist > 1.0E-06) {
            displayInverted(inflatedMesh);
            msg.setText("Info: Back Projection not recovered the original shape. Perform more ? ");
            int ret = msg.exec();
            if( ret == QMessageBox::Cancel) return;
        }
    }

    // First optimized the deflated mesh and then copy its coordinates back
    // to the source mesh.
    untangle.optimize();

    // Now the role of the infalted mesh is over. Delete is for display...
    inflatedMesh->getAttribute("Render", mrender);
    mrender->display = 0;

    /*
        JMeshNonlinearOptimization mopt;
        if( inflateOptCheckBox->isChecked() ) {
                mopt.setMesh(inflatedMesh);
                mopt.setNumIterations(100);
                mopt.setBoundaryPreserve(1);
                mopt.improveQuality();
         }
    */

    // Display the original mesh. If everything goes fine, it should be
    // inversion-free. Remember that at the end, Jacobian of some element
    // can still be negative for non-simplicial elements because of non-planarity
    // of some faces. But none of the elements should overlap or self-intersect
    // at this stage ..
    displayInverted(srcMesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshUntangleDialog :: closeDialog()
{
    if( meshViewer == nullptr) return;

    meshViewer->getFaceDraw()->setFrontFace(GL_CCW);
    meshViewer->getFaceDraw()->setCulling( GL_BACK );

    JNodeRenderPtr vAttrib;
    size_t numnodes =  srcMesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = srcMesh->getNodeAt(i);
        if( v->isActive() ) {
            v->getAttribute("Render", vAttrib);
            vAttrib->color = prevNodeColor[v];
        }
    }
    prevNodeColor.clear();

    JEdgeRenderPtr eAttrib;
    size_t numedges =  srcMesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &e = srcMesh->getEdgeAt(i);
        if( e->isActive() ) {
            e->getAttribute("Render", eAttrib);
            eAttrib->color   =  prevEdgeColor[e];
            eAttrib->display = 1;
        }
    }
    prevEdgeColor.clear();

    JFaceRenderPtr fAttrib;
    size_t numfaces =  srcMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = srcMesh->getFaceAt(i);
        if( f->isActive() ) {
            f->getAttribute("Render", fAttrib);
            fAttrib->color = prevFaceColor[f];
            fAttrib->display = 1;
        }
    }
    prevFaceColor.clear();
    meshViewer->updateBuffers(srcMesh);
    this->parentWidget()->show();
    this->close();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshUntangleDialog :: loadSrcMesh()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshUntangleDialog :: setOriginal()
{
    if( srcMesh == nullptr) return;
    srcMesh->getGeometry()->setCoordsArray(orgCoords,l2g);
    meshViewer->updateBuffers(srcMesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshUntangleDialog :: openMapQualityDialog()
{
    if( mapQualityDialog == nullptr )
        mapQualityDialog.reset(new JMeshMapQualityDialog(this));

    mapQualityDialog->setViewManager( viewManager );

    mapQualityDialog->setSourceMesh( inflatedMesh);
    mapQualityDialog->setTargetMesh( srcMesh);
    mapQualityDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshUntangleDialog :: makeConnections()
{
//    PushButton( searchPushButton,  [=] {searchInverted();});
    PushButton( resetPushButton,   [=] {setOriginal();});
    PushButton( sourcePushButton,  [=] {loadSrcMesh();});
    PushButton( applyPushButton,   [=] {execute();});
    PushButton( colorPushButton,   [=] {setColor();});
    PushButton( mapQualityPushButton, [=] {openMapQualityDialog();});
    PushButton( backProjectionPushButton, [=] {backProjection();});
    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
