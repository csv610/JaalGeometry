#include "MeshNodesDialog.hpp"

using namespace std;

#include <igl/jet.h>

///////////////////////////////////////////////////////////////////////////////

JMeshNodesDialog :: JMeshNodesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    drawNode    = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    numNodesLineEdit->setText( QString::number(0) );
    numBoundaryLineEdit->setText( QString::fromStdString("Unknown"));
    numVisibleLineEdit->setText( QString::number(0) );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;

    if( mesh == nullptr) return;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    mesh->pruneNodes();
    mesh->enumerate(0);

    size_t numNodes = mesh->getActiveSize(0);
    numNodesLineEdit->setText( QString::number(numNodes) );

    nodeIDSpinBox->setMaximum(numNodes-1);

    centerAtSpinBox->setMinimum( 0);
    centerAtSpinBox->setMaximum( numNodes-1);

    numNodes = mesh->getTopology()->getBoundarySize(0);
    numBoundaryLineEdit->setText( QString::number(numNodes) );

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    bool val = mrender->displayEntity[0];
    displayCheckBox->setChecked(val);

    setNumVisible();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    drawNode  = meshViewer->getNodeDraw();
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: openAttribListDialog()
{
    if( attribListDialog.get() == nullptr )
        attribListDialog.reset(new JMeshEntityAttribListDialog( this ));

    attribListDialog->setMesh(mesh, 0);
    attribListDialog->show();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: openPermuteDialog()
{
    if( permuteDialog.get() == nullptr )
        permuteDialog.reset(new JPermuteDialog( this ));

    permuteDialog->setViewManager( viewManager );
    permuteDialog->show();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: openNormalsDialog()
{
    if( normalsDialog.get() == nullptr )
        normalsDialog.reset(new JMeshNormalsDialog( this ));

    normalsDialog->setViewManager( viewManager );
    normalsDialog->setMesh(mesh, 0);
    normalsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: setNumVisible()
{
    if( meshViewer == nullptr ) return;
    
    size_t nCount = 0;
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender->displayEntity[0] )  {
    nCount = meshViewer->getNumVisible(0);
    }
    numVisibleLineEdit->setText( QString::number(nCount) );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: closeDialog()
{
    this->parentWidget()->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: checkNodes()
{
    if( mesh == nullptr ) return;

    double val = displayCheckBox->isChecked();
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = val;
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: setInternal()
{
    if( mesh == nullptr ) return;

    if( attribDialog.get() == nullptr )
        attribDialog.reset(new JNodeAttributesDialog( this ));

    attribDialog->setViewManager( viewManager );
    attribDialog->setMesh(mesh);

    JNodeSequence nodes;
    size_t numNodes = mesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &node = mesh->getNodeAt(i);
        if( !node->isBoundary() )
            nodes.push_back(node);
    }

    attribDialog->setNodes(nodes);
    attribDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: setBoundary()
{
    if( mesh == nullptr ) return;

    if( attribDialog.get() == nullptr )
        attribDialog.reset(new JNodeAttributesDialog( this ));

    attribDialog->setViewManager( viewManager );
    attribDialog->setMesh(mesh);

    JNodeSequence nodes;
    size_t numNodes = mesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &node = mesh->getNodeAt(i);
        if( node->isBoundary() )
            nodes.push_back(node);
    }

    attribDialog->setNodes(nodes);
    attribDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: checkDisplay()
{
    if( mesh == nullptr) return;

    bool val;
    val  = displayCheckBox->isChecked();

    size_t numVis   = 0;
    size_t numnodes = mesh->getSize(0);
    JNodeRenderPtr attrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", attrib);
            attrib->display = val;
            if( val ) numVis++;
        }
    }
    numVisibleLineEdit->setText( QString::number(numVis) );
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: displayIDs()
{
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayIDs[0] =  enumCheckBox->isChecked();
    meshViewer->refreshDisplay();

}
///////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: checkState()
{
    bool val;

    if( meshViewer == nullptr) return;
    if( drawNode == nullptr ) return;

    val  = antiAliasCheckBox->isChecked();
    drawNode->setAntiAliasing(val);

    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: getBoundary()
{
    if( mesh == nullptr ) return;
    mesh->getTopology()->searchBoundary();

    size_t numNodes = mesh->getTopology()->getBoundarySize(0);
    numBoundaryLineEdit->setText( QString::number(numNodes) );

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: setDefaultColorMethod()
{
    if( meshViewer == nullptr ) return;
    meshViewer->displayAll(0,1);
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: updateNodes()
{
}

////////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: centerAtNode()
{
    if( mesh->isEmpty() ) return;
    size_t id = centerAtSpinBox->value();
    const JNodePtr &vtx = mesh->getNodeAt(id);
    meshViewer->lookAt(vtx);
}
////////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: geomCenter()
{
    if( viewManager == nullptr || mesh == nullptr ) return;
    const Point3D &xyz = mesh->getGeometry()->getCenter();
    qglviewer::Vec vec;
    vec[0] = xyz[0];
    vec[1] = xyz[1];
    vec[2] = xyz[2];
    viewManager->setSceneCenter(vec);
    viewManager->refreshDisplay();

    /*
       Vec3F srcVec, dstVec, perpAxis;
       GLdouble mat[16];

       Vec3F normal;
       vtx->getAttribute("Normal", normal);

       Vec vec;
       vec = viewManager->camera()->position();
       double dist = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
       vec[0] = dist*normal[0];
       vec[1] = dist*normal[1];
       vec[2] = dist*normal[2];

       viewManager->camera()->setPosition(vec);

       const Point3D &p3d = vtx->getXYZCoords();
       vec[0] = p3d[0];
       vec[1] = p3d[1];
       vec[2] = p3d[2];

       viewManager->camera()->lookAt( vec );
       meshViewer->refreshDisplay();
    */
}

////////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: changeNodeID()
{
    if( mesh == nullptr) return;
    size_t id = nodeIDSpinBox->value();
    if( id > mesh->getSize(0) ) return;

    const Point3D &xyz = mesh->getNodeAt(id)->getXYZCoords();
    xCoordLineEdit->setText( QString::number(xyz[0]) );
    yCoordLineEdit->setText( QString::number(xyz[1]) );
    zCoordLineEdit->setText( QString::number(xyz[2]) );
}

////////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: showCoords()
{
    if( mesh == nullptr) return;
    size_t id = nodeIDSpinBox->value();
    if( id > mesh->getSize(0) ) return;

    const Point3D  &xyz = mesh->getNodeAt(id)->getXYZCoords();

    xCoordLineEdit->setText( QString::number(xyz[0]));
    yCoordLineEdit->setText( QString::number(xyz[1]));
    zCoordLineEdit->setText( QString::number(xyz[2]));
}
////////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: modifyCoords()
{
    if( mesh == nullptr) return;
    size_t id = nodeIDSpinBox->value();
    if( id > mesh->getSize(0) ) return;

    QString qstr;
    Point3D xyz;
    qstr = xCoordLineEdit->text();
    xyz[0] = qstr.toDouble();

    qstr = yCoordLineEdit->text();
    xyz[1] = qstr.toDouble();

    qstr = zCoordLineEdit->text();
    xyz[2] = qstr.toDouble();

    const JNodePtr &vtx  = mesh->getNodeAt(id);
    vtx->setXYZCoords(xyz);
    meshViewer->updateBuffers(mesh);
}

////////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: renumber()
{
    if( mesh == nullptr) return;
    mesh->enumerate(0);
    meshViewer->updateBuffers(mesh);
}
////////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: setColorScheme()
{
    if( mesh == nullptr) return;

    QString qstr = colorSchemeComboBox->currentText();
    Eigen::VectorXd  val;
    size_t numActive = mesh->getActiveSize(0);
    val.resize(numActive);

    size_t numnodes = mesh->getSize(0);
    size_t index = 0;

    if( qstr == "MinimumColors") {
        vector<int> nC = mesh->getTopology()->getMinNodesColor();
        for( size_t i = 0; i < numnodes; i++)
            val[i] = nC[i];
    }

    if( qstr == "XCoords") {
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() ) {
                const Point3D &p = v->getXYZCoords();
                val[index++] = p[0];
            }
        }
    }
    if( qstr == "YCoords") {
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() ) {
                const Point3D &p = v->getXYZCoords();
                val[index++] = p[1];
            }
        }
    }
    if( qstr == "ZCoords") {
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() ) {
                const Point3D &p = v->getXYZCoords();
                val[index++] = p[2];
            }
        }
    }

    Eigen::MatrixXd colors;
    igl::jet(val,true,colors);

    index = 0;
    JNodeRenderPtr vAttrib;
    JColor  color;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) {
            v->getAttribute("Render", vAttrib);
            color[0] = colors.coeff(index, 0);
            color[1] = colors.coeff(index, 1);
            color[2] = colors.coeff(index, 2);
            color[3] = 1.0;
            vAttrib->color = color;
            index++;
        }
    }
    meshViewer->updateBuffers(mesh);

}
////////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: openSamplesDialog()
{
    if( mesh == nullptr) return;
    if( samplesDialog == nullptr)
        samplesDialog.reset( new JMeshSamplesDialog(this));

    samplesDialog->setViewManager(viewManager);
    samplesDialog->setMesh(mesh);
    samplesDialog->show();
}

/////////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: getRCMOrdering()
{
    if( mesh == nullptr) return;
    JWaitCursor wCursor;
    wCursor.start();

    mesh->getTopology()->setRCMOrdering();
}

/////////////////////////////////////////////////////////////////////////////////
void JMeshNodesDialog :: makeConnections()
{
    PushButton( samplesPushButton,       [=] {openSamplesDialog();});
    PushButton( permutePushButton,       [=] {openPermuteDialog();});
    PushButton( normalsPushButton,       [=] {openNormalsDialog();});
    PushButton( enumeratePushButton,     [=] {renumber();});
    PushButton( attribListPushButton,    [=] {openAttribListDialog();});
    PushButton( geomCenterPushButton,    [=] {geomCenter();});
    PushButton( getBoundaryPushButton,   [=] {getBoundary();});
    PushButton( colorSchemePushButton,   [=] {setColorScheme();});
    PushButton( internalNodesPushButton, [=] {setInternal();});
    PushButton( boundaryNodesPushButton, [=] {setBoundary();});
    PushButton( rcmOrderPushButton,   [=] {getRCMOrdering();});
    PushButton( lookAtPushButton,        [=] { centerAtNode();});

    CheckBox( displayCheckBox, [=] {checkDisplay();});
    CheckBox( enumCheckBox,    [=] {displayIDs();});
    CheckBox( antiAliasCheckBox, [=] {checkState();});

    SpinBoxi( nodeIDSpinBox,    [=] { showCoords();});

    LineEdit( xCoordLineEdit, [=] { modifyCoords();});
    LineEdit( yCoordLineEdit, [=] { modifyCoords();});
    LineEdit( zCoordLineEdit, [=] { modifyCoords();});

    PushButton( closePushButton, [=] {closeDialog();});
}
///////////////////////////////////////////////////////////////////////////////
/*
void JMeshNodesDialog :: changeCenter( bool refresh )
{
     if( mesh == nullptr ) return;

     if( showOneCheckBox->isChecked() ) {
          if( changeCenterCheckBox->isChecked() ) {
               const Point3D &xyz = currNode->getXYZCoords();
               AffineTransform af(mesh);
               af.translate(-xyz[0], -xyz[1], -xyz[2]);
          }

          if( refresh ) meshViewer->refreshDisplay();
     }
}
*/
///////////////////////////////////////////////////////////////////////////////

/*
void JMeshNodesDialog :: getNodeID()
{
     if( mesh == nullptr ) return;

     size_t currNodeID  = nodeIDSlider->value();
     Vertex *vtx =  mesh->getNodeAt( currNodeID ) ;
     activate( vtx );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNodesDialog :: reset_neighs(bool val )
{
     if( mesh == nullptr ) return;
     if( currNode == nullptr ) return;

     if( cellNeighsCheckBox->isChecked() ) {
          if( mesh->getAdjTable(0,3) == 0) mesh->buildRelations(0,3);

          meshViewer->displayAll( 3, 1);
          JCellSequence cneighs;
          currNode->getRelations(cneighs);
          for( size_t i = 0; i < cneighs.size(); i++) {
               Cell *cell = cneighs[i];
               cell->setAttribute("Display", val);
          }

          if( faceNeighsCheckBox->isChecked() ) {
               meshViewer->displayAll( 1, 1);
               for( size_t i = 0; i < cneighs.size(); i++) {
                    Cell *cell = cneighs[i];
                    for( int j = 0; j < cell->getSize(2); j++) {
                         Face *face = cell->getFaceAt(j);
                         face->setAttribute("Display", val);
                    }
               }
          }

          if( edgeNeighsCheckBox->isChecked() ) {
               for( size_t i = 0; i < cneighs.size(); i++) {
                    Cell *cell = cneighs[i];
                    for( int j = 0; j < cell->getSize(1); j++) {
                         Edge *edge = cell->getEdgeAt(j);
                         edge->setAttribute("Display", val);
                    }
               }
          }

          if( nodeNeighsCheckBox->isChecked() ) {
               for( size_t i = 0; i < cneighs.size(); i++) {
                    Cell *cell = cneighs[i];
                    for( int j = 0; j < cell->getSize(0); j++) {
                         Vertex *vtx = cell->getNodeAt(j);
                         vtx->setAttribute("Display", val);
                    }
               }
          }
     }

     if( faceNeighsCheckBox->isChecked() ) {
          if( mesh->getAdjTable(0,2) == 0) mesh->buildRelations(0,2);
          meshViewer->displayAll( 2, 1);

          JFaceSequence fneighs;
          currNode->getRelations(fneighs);
          for( size_t i = 0; i < fneighs.size(); i++) {
               Face *face = fneighs[i];
               face->setAttribute("Display", val);
          }

          if( edgeNeighsCheckBox->isChecked() ) {
               for( size_t i = 0; i < fneighs.size(); i++) {
                    Face *face = fneighs[i];
                    for( int j = 0; j < face->getSize(1); j++) {
                         Edge *edge= face->getEdgeAt(j);
                         edge->setAttribute("Display", val);
                    }
               }
          }

          if( nodeNeighsCheckBox->isChecked() ) {
               for( size_t i = 0; i < fneighs.size(); i++) {
                    Face *face = fneighs[i];
                    for( int j = 0; j < face->getSize(0); j++) {
                         Vertex *vtx = face->getNodeAt(j);
                         vtx->setAttribute("Display", val);
                    }
               }
          }
     }


     if( edgeNeighsCheckBox->isChecked() ) {
          if( mesh->getAdjTable(0,1) == 0) mesh->buildRelations(0,1);

          meshViewer->displayAll(1, 1);
          JEdgeSequence eneighs;
          currNode->getRelations(eneighs);
          for( size_t i = 0; i < eneighs.size(); i++) {
               Edge *edge = eneighs[i];
               edge->setAttribute("Display", val);
          }

          if( nodeNeighsCheckBox->isChecked() ) {
               for( size_t i = 0; i < eneighs.size(); i++) {
                    Edge *edge = eneighs[i];
                    Vertex *v0 = edge->getNodeAt(0);
                    Vertex *v1 = edge->getNodeAt(1);
                    v0->setAttribute("Display", val);
                    v1->setAttribute("Display", val);
               }
          }
     }

     meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshNodesDialog :: showNeighs()
{
     if( mesh == nullptr ) return;

     if( !showOneCheckBox->isChecked() )  {
          meshViewer->displayAll(1);
          return;
     }

     int val = nodeNeighsCheckBox->isChecked() + edgeNeighsCheckBox->isChecked()  +
               faceNeighsCheckBox->isChecked() + cellNeighsCheckBox->isChecked();

     if( val ) {
          meshViewer->displayAll(0);
          reset_neighs(1);
          activate( currNode );
     }

     meshViewer->refreshDisplay();
}

void JMeshNodesDialog :: alignNormal( bool refresh )
{
     if( mesh == nullptr ) return;

     if( !showOneCheckBox->isChecked() ) return;
     if( !alignCheckBox->isChecked()  ) return;

     if( currNode == nullptr ) return;

     AffineTransform af(mesh);
     Point3D p3d = currNode->getXYZCoords();
     af.translate(-p3d[0], -p3d[1], -p3d[2] );

     mesh->getGeometry()->setFacesNormal();
     mesh->getGeometry()->setNodesNormal();

     Vec3D zAxis;
     zAxis[0] = 0.0;
     zAxis[1] = 0.0;
     zAxis[2] = 1.0;

     Point3F nr;
     currNode->getAttribute("Normal", nr);

     Point3D normal;
     normal[0] = nr[0];
     normal[1] = nr[1];
     normal[2] = nr[2];

     Vec3D  perpAxis;
     JMath::cross_product( zAxis, normal, perpAxis);

     double angle = JMath::getVecAngle(normal, zAxis, ANGLE_IN_RADIANS);

     qglviewer::Vec axis(perpAxis[0], perpAxis[1], perpAxis[2] );

     qglviewer::Quaternion quaternion(axis, -1.0*angle);

     qglviewer::Vec prot;

     size_t numNodes = mesh->getSize(0);
     for( size_t i = 0; i < numNodes; i++) {
          Vertex *v = mesh->getNodeAt(i);
          Point3D p3d = v->getXYZCoords();
          prot[0] = p3d[0];
          prot[1] = p3d[1];
          prot[2] = p3d[2];
          prot    = quaternion.rotate(prot);
          p3d[0]  = prot[0];
          p3d[1]  = prot[1];
          p3d[2]  = prot[2];
          v->setXYZCoords(p3d);
     }

     mesh->getGeometry()->setFacesNormal();
     mesh->getGeometry()->setNodesNormal();
     currNode->getAttribute("Normal", nr);

     if( refresh ) meshViewer->refreshDisplay();
}
*/

