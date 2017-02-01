#include "MeshGeodesicsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshGeodesicsDialog :: JMeshGeodesicsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    oneSrcDst= 1;
    timeStepLineEdit->setText( QString::number(1.0) );
}

///////////////////////////////////////////////////////////////////////////////

JMeshGeodesicsDialog :: ~JMeshGeodesicsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshGeodesicsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

    picker = meshViewer->getEntityPicker();

    if( picker ) picker->setMode(1);

    viewManager->attach( this );
    viewManager->setMouseTracking(1);
    viewManager->setSelectRegionHeight(10);
    viewManager->setSelectRegionWidth(10);

    setMesh( meshViewer->getCurrentMesh() );
}
///////////////////////////////////////////////////////////////////////////////

void JMeshGeodesicsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str()));

    jgeodesic.reset( new JTriMeshGeodesics());
    jgeodesic->setMesh( mesh );
//  jgeodesic->initialize();
    srcSpinBox->setMinimum( 0 );
    dstSpinBox->setMinimum( 0 );
    srcSpinBox->setMaximum( mesh->getSize(0));
    dstSpinBox->setMaximum( mesh->getSize(0));

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 0;
    mrender->setSurfaceShade(JRender::SMOOTH_SHADE);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshGeodesicsDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }

    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshGeodesicsDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( !oneSrcDst ) return;
    if( picker == nullptr ) return;

    JNodeSequence  nodes = picker->getPickedNodes();
    if( nodes.empty() ) return;

    cout << "Node Picked " << endl;

/*
    JNodePtr currPicked = nodes.back();
    JNodeRenderPtr  nAttrib;
    currPicked->getAttribute("Render", nAttrib);
    nAttrib->glyph  = JNodeRender::NODE_AS_SPHERE;

    // If the node already picked earlier, we deselect it ...
    if( srcPickingRadioButton->isChecked() ) {
        srcSpinBox->setValue( currPicked->getID() );
    if( srcNodes.find(currPicked) != srcNodes.end()) {
        srcNodes.erase(currPicked);
        nAttrib->scale  = 1.0;
        nAttrib->glyph  = JNodeRender::NODE_AS_POINT;
        picker->clearAll();
    }
        srcNodes.insert( nodes.back() );
        nAttrib->color = JEntityColor::getColor("Red");
    }

    if( dstPickingRadioButton->isChecked() ) {
        dstSpinBox->setValue( currPicked->getID() );
    if( dstNodes.find(currPicked) != dstNodes.end()) {
        dstNodes.erase(currPicked);
        nAttrib->scale  = 1.0;
        nAttrib->glyph  = JNodeRender::NODE_AS_POINT;
        picker->clearAll();
    }
        dstNodes.insert( nodes.back() );
        nAttrib->color = JEntityColor::getColor("Blue");
    }

    picker->clearAll();
    meshViewer->getNodeDraw()->updateBuffers(mesh);
    meshViewer->refreshDisplay();
*/
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeodesicsDialog :: setSrcDstColor( const JNodePtr &vertex, int src_or_dst)
{
    JColor clr;
    clr[0] = 0.0;
    clr[1] = 0.0;
    clr[2] = 0.0;
    clr[3] = 0.0;
    if( src_or_dst == 0) clr[0] = 1.0;
    if( src_or_dst == 1) clr[2] = 1.0;

    JNodeRenderPtr nAttrib;
    double radius = meshViewer->getNodeDraw()->getBallRadius();
    vertex->getAttribute("Render", nAttrib);
    nAttrib->scale  = 1.5;
    nAttrib->glyph  = JNodeRender::NODE_AS_SPHERE;
    nAttrib->color  = clr;
}

///////////////////////////////////////////////////////////////////////////////
void JMeshGeodesicsDialog :: getPath( const JNodePtr &src, const JNodePtr &dst)
{
/*
    if( src == dst) return;
        JEdgeSequence geoedges;
        jgeodesic->getPath( src, dst, geoedges);
        if( geoedges.empty() ) return;

        pathEdges.push_back(geoedges);
        double geolen = JEdgeGeometry::getLength(geoedges);
        lengthLineEdit->setText( QString::number(geolen) );

        double elen = JNodeGeometry::getLength(src, dst);
        double dilate = geolen/elen;
        dilationLineEdit->setText( QString::number(dilate) );

        numSegmentsLineEdit->setText( QString::number(geoedges.size() ) );

        setSrcDstColor( src, 0);
        setSrcDstColor( dst, 1);

        JEdgeRenderPtr eAttrib;
        JColor clr;
        clr[0] = 1.0;
        clr[1] = 1.0;
        clr[2] = 1.0;
        clr[3] = 1.0;
        bool val = 1;
        for( size_t i = 0; i < geoedges.size(); i++) {
            geoedges[i]->getAttribute("Render", eAttrib);
            geoedges[i]->setAttribute("Interface", val);
            eAttrib->scale = 1.5;
            eAttrib->color = clr;
            eAttrib->display = 1;
        }
        meshViewer->updateBuffers(mesh);
*/
}
///////////////////////////////////////////////////////////////////////////////

void JMeshGeodesicsDialog :: getDijkstraPath()
{
    if( meshViewer == nullptr ) return;

    /*
        JNodeSequence srcnodes, dstnodes;
        std::copy(srcNodes.begin(), srcNodes.end(), back_inserter(srcnodes));
        dist = jgeodesic->getDistanceField( srcnodes);
        JNodeColor::assign(mesh, dist);
    */

    /*
            if( distanceMapRadioButton->isChecked() ) {
                if( farthestPointsCheckBox->isChecked() )
                    farthestPointSampling();
                for( size_t i = 0; i < srcnodes.size(); i++)
                    setSrcDstColor( srcnodes[i], 0);
            }

            if( voronoiRadioButton->isChecked() ) {
                meshViewer->getDrawFace()->setShade( DrawFace::SMOOTH_SHADE);
                QString qstr = numFarthestPointsLineEdit->text() ;
                int  numRegions  = qstr.toInt();
                if( numRegions > 1) {
                    jgeodesic->getVoronoi( numRegions, srcnodes);
                    NodePartitionColor nodeColor;
                    NodeColor::assign(mesh, &nodeColor);
                }
                for( size_t i = 0; i < srcnodes.size(); i++)
                    setSrcDstColor( srcnodes[i], 0);
            }


            if( boundSrcRadioButton->isChecked() ) {
            }

            if( extremaSrcRadioButton->isChecked() ) {
            }

            if( infyDstRadioButton->isChecked() ) {
            }

            QMessageBox msg;
            if( oneSrcOneDstRadioButton->isChecked() ) {
                if( extremeSrcDstRadioButton->isChecked() ) {
                    Vertex *vsrc = nullptr;
                    Vertex *vdst = nullptr;
                    jgeodesic->getExtremes( vsrc, vdst);
                    if( vsrc ) usrSrcSet.insert(vsrc);
                    if( vdst ) usrDstSet.insert(vdst);
                }

                meshViewer->getDrawFace()->setShade( DrawFace::FLAT_SHADE);
                if( usrSrcSet.empty()) {
                    msg.setIcon(QMessageBox::Warning);
                    msg.setText("No source node selected: Select one node from the model");
                    msg.setStandardButtons( QMessageBox::Ok);
                    int ret = msg.exec();
                    if( ret == QMessageBox::Ok) return;
                }

                if( usrSrcSet.size() > 1 ) {
                    msg.setIcon(QMessageBox::Warning);
                    msg.setText("Only one source is required: The last one will be selected");
                    msg.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel );
                    int ret = msg.exec();
                    if( ret == QMessageBox::Cancel) return;
                }

                if( usrDstSet.empty()) {
                    msg.setIcon(QMessageBox::Warning);
                    msg.setText("No destination node selected: Select one node from the model");
                    msg.setStandardButtons( QMessageBox::Ok);
                    int ret = msg.exec();
                    if( ret == QMessageBox::Ok) return;
                }

                if( usrDstSet.size() > 1) {
                    msg.setIcon(QMessageBox::Warning);
                    msg.setText("Only one destination node is required: The last one will be selected");
                    msg.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel );
                    int ret = msg.exec();
                    if( ret == QMessageBox::Cancel) return;
                    return;
                }

                std::copy(usrSrcSet.begin(), usrSrcSet.end(), back_inserter(srcnodes));
                std::copy(usrDstSet.begin(), usrDstSet.end(), back_inserter(dstnodes));
                JEdgeSequence geoedges;
                jgeodesic->getPath( srcnodes.back(), dstnodes.back(), geoedges);

                for( size_t i = 0; i < srcnodes.size(); i++)
                    setSrcDstColor( srcnodes[i], 0);

                for( size_t i = 0; i < dstnodes.size(); i++)
                    setSrcDstColor( dstnodes[i], 1);

                RenderEdgeAttribute *eAttrib = nullptr;
                Color clr;
                clr[0] = 1.0;
                clr[1] = 1.0;
                clr[2] = 1.0;
                clr[3] = 1.0;
                bool val = 1;
                for( size_t i = 0; i < geoedges.size(); i++) {
                    geoedges[i]->getAttribute("Render", eAttrib);
                    geoedges[i]->setAttribute("Interface", val);
                    eAttrib->scale = 1.5;
                    eAttrib->color = clr;
                }
                usrSrcSet.clear();
                usrDstSet.clear();
            }

            if( manySrcOneDstRadioButton->isChecked() ) {
                meshViewer->getDrawFace()->setShade( DrawFace::FLAT_SHADE);
                if( usrSrcSet.empty()) {
                    msg.setIcon(QMessageBox::Warning);
                    msg.setText("No source node selected: Select at least one node from the model");
                    msg.setStandardButtons( QMessageBox::Ok);
                    int ret = msg.exec();
                    if( ret == QMessageBox::Ok) return;
                }

                if( usrDstSet.empty()) {
                    msg.setIcon(QMessageBox::Warning);
                    msg.setText("No destination node selected: Select one node from the model");
                    msg.setStandardButtons( QMessageBox::Ok);
                    int ret = msg.exec();
                    if( ret == QMessageBox::Ok) return;
                }

                if( usrDstSet.size() > 1) {
                    msg.setIcon(QMessageBox::Warning);
                    msg.setText("Only one destination node is required: The last one will be selected");
                    msg.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel );
                    int ret = msg.exec();
                    if( ret == QMessageBox::Cancel) return;
                    return;
                }

                std::copy(usrSrcSet.begin(), usrSrcSet.end(), back_inserter(srcnodes));
                std::copy(usrDstSet.begin(), usrDstSet.end(), back_inserter(dstnodes));
                for( size_t j = 0; j < srcnodes.size(); j++) {
                    JEdgeSequence geoedges;
                    jgeodesic->getPath( srcnodes[j], dstnodes.back(), geoedges);
                    RenderEdgeAttribute *eAttrib = nullptr;
                    Color clr;
                    clr[0] = 1.0;
                    clr[1] = 1.0;
                    clr[2] = 1.0;
                    clr[3] = 1.0;
                    for( size_t i = 0; i <  geoedges.size(); i++) {
                        geoedges[i]->getAttribute("Render", eAttrib);
                        eAttrib->scale = 2.0;
                        eAttrib->color = clr;
                    }
                }
                usrSrcSet.clear();
                usrDstSet.clear();
            }
        */
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeodesicsDialog :: donePicking()
{
    /*
        if( meshViewer == nullptr ) return;
            JNodeSequence nodes;

    //        JNodeRenderPtr nAttrib = nullptr;
    //        Color color;
            if( srcPickingRadioButton->isChecked() ) {
                nodes = picker->getPickedNodes();
                for( size_t i = 0; i < nodes.size(); i++) {
                    usrSrcSet.insert( nodes[i] );
                    setSrcDstColor( nodes[i], 0);
                }
                picker->clearAll();

                if( !usrSrcSet.empty() )
                    userSrcRadioButton->setChecked( true );
            }

            if( dstPickingRadioButton->isChecked() ) {
                nodes = picker->getPickedNodes();
                for( size_t i = 0; i < nodes.size(); i++) {
                    usrDstSet.insert( nodes[i] );
                    setSrcDstColor( nodes[i], 1);
                }
                if( !usrDstSet.empty() )
                    usrDstRadioButton->setChecked( true );
            }
    //     meshViewer->refreshDisplay();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeodesicsDialog :: deleteLastSegment()
{
    if( meshViewer == nullptr ) return;

    if( pathEdges.empty() ) return;

    JEdgeSequence geoedges = pathEdges.back();

    pathEdges.pop_back();

    JEdgeRenderPtr eAttrib;
    JColor clr;
    clr[0] = 0.0;
    clr[1] = 0.0;
    clr[2] = 0.0;
    clr[3] = 1.0;
    bool val = 1;
    for( size_t i = 0; i < geoedges.size(); i++) {
        geoedges[i]->getAttribute("Render", eAttrib);
        geoedges[i]->setAttribute("Interface", val);
        eAttrib->scale = 1.0;
        eAttrib->color = clr;
        eAttrib->display = 1;
    }

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshGeodesicsDialog :: deleteAll()
{
    if( meshViewer == nullptr ) return;
    pickedNodes.clear();

    if( pathEdges.empty() ) return;

    int nSegments = pathEdges.size();
    for( int j = 0; j < nSegments; j++) {
        JEdgeSequence geoedges = pathEdges.back();
        pathEdges.pop_back();

        JEdgeRenderPtr eAttrib;
        JColor clr;
        clr[0] = 0.0;
        clr[1] = 0.0;
        clr[2] = 0.0;
        clr[3] = 1.0;
        bool val = 1;
        for( size_t i = 0; i < geoedges.size(); i++) {
            geoedges[i]->getAttribute("Render", eAttrib);
            geoedges[i]->setAttribute("Interface", val);
            eAttrib->scale = 1.0;
            eAttrib->color = clr;
            eAttrib->display = 1;
        }
    }

    meshViewer->updateBuffers(mesh);
}

/*
void JMeshGeodesicsDialog :: farthestPointSampling()
{
    JNodeSequence nodes;

    QString qstr = numFarthestPointsLineEdit->text() ;
    int  numSamples  = qstr.toInt();
    if( numSamples ) {
        jgeodesic->getFarthestSamples( numSamples, nodes, dist);
        for( size_t i = 0; i < nodes.size(); i++)
            srcSet.insert( nodes[i] );
    }
}
*/

///////////////////////////////////////////////////////////////////////////////
void JMeshGeodesicsDialog :: openNodeAttribDialog()
{
    if( nodeAttribDialog == nullptr )
        nodeAttribDialog.reset(new JNodeAttributesDialog());

    JNodeSet nset;
    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++)  {
        JEdgePtr edge = mesh->getEdgeAt(i);
        if( edge->hasAttribute("Interface") ) {
            nset.insert( edge->getNodeAt(0));
            nset.insert( edge->getNodeAt(1));
        }
    }

    if( nset.size() ) {
        JNodeSequence nodes( nset.size() );
        std::copy( nset.begin(), nset.end(), nodes.begin() );
        nodeAttribDialog->setNodes(nodes);
        nodeAttribDialog->show();
        meshViewer->refreshDisplay();
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeodesicsDialog :: openEdgeAttribDialog()
{

    if( edgeAttribDialog == nullptr )
        edgeAttribDialog.reset( new JEdgeAttributesDialog());

    JEdgeSequence edges;
    size_t numEdges = mesh->getSize(1);
    JEdgeRenderPtr eAttrib;
    for( size_t i = 0; i < numEdges; i++)  {
        JEdgePtr edge = mesh->getEdgeAt(i);
        if( edge->hasAttribute("Interface") ) {
            edges.push_back(edge);
            int err = edge->getAttribute("Render", eAttrib);
            if( !err ) eAttrib->display = 1;
        }
    }

    edgeAttribDialog->setEdges(edges);
    edgeAttribDialog->show();
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeodesicsDialog :: checkState()
{
    /*
        if( meshViewer == nullptr ) return;
        bool val = farthestPointsCheckBox->isChecked();
        distanceMapRadioButton->setChecked( val);

        val = oneSrcOneDstRadioButton->isChecked();
        if( val ) meshViewer->getDrawFace()->setShade( DrawFace::FLAT_SHADE);

        val = manySrcOneDstRadioButton->isChecked();
        if( val ) meshViewer->getDrawFace()->setShade( DrawFace::FLAT_SHADE);

        val = distanceMapRadioButton->isChecked();
        if( val ) meshViewer->getDrawFace()->setShade( DrawFace::SMOOTH_SHADE);
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshGeodesicsDialog :: startNewPath()
{
    if( mesh == nullptr ) return;

    int srcid = srcSpinBox->value();  // Select source node ...
    int dstid = dstSpinBox->value();  // Select target node ...
    JNodePtr vsrc = mesh->getNodeAt(srcid);
    JNodePtr vdst = mesh->getNodeAt(dstid);

    getPath(vsrc,vdst);
    picker->clearAll();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshGeodesicsDialog :: closeDialog()
{
    if( mesh ) {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->pickableEntity = -1;
    }

    if( viewManager) {
        viewManager->detach( this );
        viewManager->setMouseTracking(0);
        if( meshViewer ) {
            JMeshEntityPickerPtr picker = meshViewer->getEntityPicker();
            if( picker ) picker->clearAll();
        }
    }
    this->hide();
    parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////
/*
void JMeshGeodesicsDialog :: getFarthestNode()
{
   if( mesh == nullptr) return;

   int id = srcSpinBox->value();
   const JNodePtr &src =  mesh->getNodeAt(id);
   JNodePtr dst;
// dst = jgeodesic->getFarthestNode(src);
   if( dst ) dstSpinBox->setValue( dst->getID() );
   getPath(src, dst);


   double val;
   vector<double>  darray;
   size_t numnodes = mesh->getSize(0);
   darray.reserve(numnodes);
   for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Distance", val);
            darray.push_back(val);
        }
   }

   vector<JColor> colors;
   JColorMap::jet(darray, colors);
   JNodeRenderPtr nattrib;

   int index = 0;
   for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", nattrib);
            nattrib->color = colors[index++];
        }
   }
   meshViewer->getFaceDraw()->setRenderMode( JRenderMode::SMOOTH_SHADE);
   meshViewer->updateBuffers(mesh);
}
*/
///////////////////////////////////////////////////////////////////////////////
void JMeshGeodesicsDialog :: setNodesDistanceColor()
{
    if( mesh == nullptr) return;

    size_t numnodes = mesh->getSize(0);

    Eigen::VectorXd K;
    K.resize(numnodes);
    size_t index = 0;
    double val;

    vector<double> vec;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Distance", val);
            K(index++) = val;
            vec.push_back(val);
        }
    }
    sort( vec.begin(), vec.end() );

    Eigen::MatrixXd Color;
    igl::jet(K, true, Color);

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->setSurfaceShade(JRender::SMOOTH_SHADE);
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshGeodesicsDialog :: heatFlow()
{
    if( mesh == nullptr) return;

    JNodeSequence srcnodes = picker->getPickedNodes();

    QMessageBox msg;
    if( srcnodes.empty()) {
        msg.setIcon(QMessageBox::Warning);
        msg.setText("No source node selected: Select one node from the model");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok) return;
    }
    QString qstr = timeStepLineEdit->text();
    double  dt   = qstr.toDouble();
    DDG::JMeshGeodesics ddg;
    ddg.setMesh(mesh);
    ddg.setTimeStep(dt);
    ddg.setSource(srcnodes);
    ddg.execute();

//  setNodesColor();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshGeodesicsDialog :: makeConnections()
{
    PushButton( dijkstraPushButton, [=] {getDijkstraPath();});
    PushButton( deleteLastPushButton, [=] {deleteLastSegment() ;});
    PushButton( deleteAllPushButton, [=] {deleteAll();});
    PushButton( closePushButton, [=] {closeDialog();});
    PushButton( applyHeatFlowPushButton, [=] {heatFlow();});
}
///////////////////////////////////////////////////////////////////////////////
