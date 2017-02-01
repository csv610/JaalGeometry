#include "MeshFacesDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshFacesDialog :: JMeshFacesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    drawFace = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    numFacesLineEdit->setText( QString::number(0) );
    numBoundaryLineEdit->setText( QString::fromStdString("Unknown"));
    numVisibleLineEdit->setText( QString::number(0) );

    /*
         screenPatternSpinBox->setMinimum(0);
         screenPatternSpinBox->setMaximum(17);
         screenPatternSpinBox->setValue(1);
         alphaSpinBox->setValue(0.5);
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFacesDialog :: displayFaceType()
{
    if( mesh == nullptr) return;
    int err, nCount = 0;

    bool displayTriangles = displayTrianglesCheckBox->isChecked();
    bool displayQuads     = displayQuadsCheckBox->isChecked();
    bool displayPolygons  = displayPolygonsCheckBox->isChecked();

    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        if( f->isActive() ) {
            err = f->getAttribute("Render", fAttrib);
            if( fAttrib) {
                if( f->getSize(0) == 3) fAttrib->display = displayTriangles;
                if( f->getSize(0) == 4) fAttrib->display = displayQuads;
                if( f->getSize(0) >  4) fAttrib->display = displayPolygons;
            }
        }
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshFacesDialog :: displayFaces()
{
    if( mesh == nullptr) return;
    int err, nCount = 0;

    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);

    if( showAllRadioButton->isChecked() ) {
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() ) {
                err = f->getAttribute("Render", fAttrib);
                if( fAttrib) fAttrib->display = 1;
            }
        }
        if( meshViewer ) meshViewer->updateBuffers(mesh);
        return;
    }

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        if( f->isActive() ) {
            err = f->getAttribute("Render", fAttrib);
            if( fAttrib) fAttrib->display = 0;
        }
    }

    if( showInvertedRadioButton->isChecked() ) {
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() ) {
                if( JFaceGeometry::isInverted(f) ) {
                    err = f->getAttribute("Render", fAttrib);
                    if( fAttrib) {
                        fAttrib->display = 1;
                        nCount++;
                    }
                }
            }
        }
        if( meshViewer ) meshViewer->updateBuffers(mesh);
        return;
    }

    if( showConvexRadioButton->isChecked() ) {
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() ) {
                if( JFaceGeometry::isConvex(f) ) {
                    err = f->getAttribute("Render", fAttrib);
                    if( fAttrib) {
                        fAttrib->display = 1;
                        nCount++;
                    }
                }
            }
        }
        if( meshViewer ) meshViewer->updateBuffers(mesh);
        return;
    }

    if( showConcaveRadioButton->isChecked() ) {
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() ) {
                if( !JFaceGeometry::isConvex(f) ) {
                    err = f->getAttribute("Render", fAttrib);
                    if( fAttrib) {
                        fAttrib->display = 1;
                        nCount++;
                    }
                }
            }
        }
        if( meshViewer ) meshViewer->updateBuffers(mesh);
        return;
    }

}
///////////////////////////////////////////////////////////////////////////////
void JMeshFacesDialog :: setMesh(const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    lookAtSpinBox->setMinimum(0);
    lookAtSpinBox->setMaximum(mesh->getSize(2) );

    size_t numfaces = mesh->getActiveSize(2);
    numFacesLineEdit->setText( QString::number(numfaces) );

    numfaces = mesh->getTopology()->getBoundarySize(2);
    numBoundaryLineEdit->setText( QString::number(numfaces) );
    setNumVisible();

    size_t numtris  = 0;
    size_t numquads = 0;
    size_t numpolys = 0;

    numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        int n = f->getSize(0);
        switch(n)
        {
        case 3:
            numtris++;
            break;
        case 4:
            numquads++;
            break;
        default:
            numpolys++;
            break;
        }
    }
    numTrianglesLineEdit->setText( QString::number(numtris));
    numQuadsLineEdit->setText( QString::number(numquads));
    numPolygonsLineEdit->setText( QString::number(numpolys));

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    bool val = mrender->displayEntity[1];
    displayCheckBox->setChecked(val);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

    drawFace  = meshViewer->getFaceDraw();
    /*
    faceColor = (MeshFaceColor *)drawFace->getDefaultColorMethod();

    currFace = nullptr;
         faceIDSlider->setMinimum(0);
         faceIDSlider->setMaximum(numfaces-1);
    //       faceIDSlider->setValue(0);
         bool val = meshViewer->isEnabled(2);
         displayFacesCheckBox->setChecked( val );
    */

}
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: lookAt()
{
    if( viewManager == nullptr || mesh == nullptr ) return;

    size_t id = lookAtSpinBox->value();
    const JFacePtr &face = mesh->getFaceAt(id);
    const JNodePtr &vtx =  face->getNodeAt(0);
    meshViewer->lookAt(vtx);

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

///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: openAttribListDialog()
{
    if( attribListDialog.get() == nullptr )
        attribListDialog.reset(new JMeshEntityAttribListDialog( this ));

    attribListDialog->setMesh( mesh, 2);
    attribListDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: checkLights()
{
    if( meshViewer == nullptr ) return;
    bool val = lightsCheckBox->isChecked();
    meshViewer->getFaceDraw()->setLights(val);
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFacesDialog :: setInternal()
{
    if( mesh == nullptr) return;

    if( attribDialog.get() == nullptr )
        attribDialog.reset(new JFaceAttributesDialog( this ));

    attribDialog->setViewManager(viewManager);
    attribDialog->setMesh(mesh);

    JFaceSequence seq;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr & f = mesh->getFaceAt(i);
        if( f->isActive() && !f->isBoundary()  ) seq.push_back(f);
    }
    attribDialog->setFaces(seq);
    attribDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFacesDialog :: setInterface()
{
    if( mesh == nullptr) return;

    if( attribDialog.get()== nullptr )
        attribDialog.reset(new JFaceAttributesDialog( this ));

    attribDialog->setViewManager(viewManager);
    attribDialog->setMesh( mesh );

    JFaceSequence seq;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        if( f->isActive() && f->hasAttribute("Interface") ) seq.push_back(f);
    }
    attribDialog->setFaces(seq);
    attribDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: setBoundary()
{
    if( mesh == nullptr) return;

    if( attribDialog.get() == nullptr )
        attribDialog.reset(new JFaceAttributesDialog( this ));

    attribDialog->setViewManager(viewManager);
    attribDialog->setMesh(mesh);

    JFaceSequence seq;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        if( f->isActive() && f->isBoundary()  ) seq.push_back(f);
    }
    attribDialog->setFaces(seq);
    attribDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: setBoundaryMaterial()
{
JMeshMaterialDialog *matdialog = new JMeshMaterialDialog();
matdialog->show();

     QColor color = QColorDialog::getColor();

     if( meshViewer == nullptr ) return;

     float rgb[3];
     rgb[0] = color.red()/255.0;
     rgb[1] = color.green()/255.0;
     rgb[2] = color.blue()/255.0;
     if( faceColor ) {
          faceColor->setBoundaryColor(rgb);
          if( uniform_color )
               faceColor->setInternalColor(rgb);
     }
     meshViewer->refreshDisplay();
}
*/
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: closeDialog()
{
    this->parentWidget()->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: setTransparency()
{
     if( meshViewer == nullptr ) return;

     if( opaqueRadioButton->isChecked() ) {
          if( drawFace )
               drawFace->setTransparencyMethod( DrawFace::OPAQUE );
     }

     if( screendoorRadioButton->isChecked() ) {
          if( drawFace ) {
               drawFace->setTransparencyMethod( DrawFace::SCREENDOOR );
               int p = screenPatternSpinBox->value();
               drawFace->setScreendoorPattern(p);
          }
     }

     if( blendRadioButton->isChecked() ) {
          if( drawFace ) {
               drawFace->setTransparencyMethod( DrawFace::BLENDING );
               double a  = alphaSpinBox->value();
               if( faceColor ) faceColor->setAlpha( a);
          }
     }
     meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: remesh()
{
    if( mesh == nullptr ) return;
    if( trimesherDialog.get() == nullptr )
        trimesherDialog.reset(new JTriMesherDialog( this ));

    trimesherDialog->show();
    this->hide();

    return;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: genVolMesh()
{
    if( mesh == nullptr ) return;

    if( tetmesherDialog.get() == nullptr )
        tetmesherDialog.reset(new JTetMesherDialog( this ));

         if( tetmeshRadioButton->isChecked() ) {
                    tetmesherDialog->show();
              this->hide();
         }
    return;
}

///////////////////////////////////////////////////////////////////////////////
*/

void JMeshFacesDialog :: setNumVisible()
{
    if( meshViewer == nullptr ) return;
    size_t nCount = meshViewer->getNumVisible(2);
    numVisibleLineEdit->setText( QString::number(nCount) );
    bool val = nCount == 0 ? 0:1;
    displayCheckBox->setChecked(val);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFacesDialog :: checkDisplay()
{
    if( mesh == nullptr) return;

    bool val;
    val  = displayCheckBox->isChecked();

    size_t numfaces = mesh->getSize(2);
    size_t numVis   = 0;
    JFaceRenderPtr attrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int err  = face->getAttribute("Render", attrib);
            if( !err) {
                attrib->display = val;
                if( val ) numVis++;
            }
        }
    }
    meshViewer->updateBuffers(mesh);
    numVisibleLineEdit->setText( QString::number(numVis) );
}
///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: checkFaces()
{
     if( mesh == nullptr ) return;

     bool val;
     val = displayFacesCheckBox->isChecked();
     meshViewer->displayAll(2,val);
     size_t numFaces = mesh->getSize(2);
     for( size_t i = 0; i < numFaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setAttribute("Display", val);
     }

     setNumVisible();
     meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: checkInternal()
{
     if( mesh == nullptr ) return;

     bool val;
     val = internalCheckBox->isChecked();
     faceColor->displayInternal(val);

     size_t numFaces = mesh->getSize(2);
     for( size_t i = 0; i < numFaces; i++) {
          Face *face = mesh->getFaceAt(i);
          if( !face->isBoundary() )
               face->setAttribute("Display", val);
     }

     setNumVisible();
     meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: checkBoundary()
{
     if( mesh == nullptr ) return;

     bool val;
     val = boundaryCheckBox->isChecked();
     faceColor->displayBoundary(val);

     size_t numFaces = mesh->getSize(2);
     for( size_t i = 0; i < numFaces; i++) {
          Face *face = mesh->getFaceAt(i);
          if( face->isBoundary() )
               face->setAttribute("Display", val);
     }
     setNumVisible();
     meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: checkState()
{
    if( meshViewer == nullptr ) return;

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayIDs[2] =  enumCheckBox->isChecked();

    /*
    bool val;
    val = displayCheckBox->isChecked();
    meshViewer->setDisplayEntity(2,val);

    if( drawFace == nullptr) return;

            val = filledCheckBox->isChecked();


            if( val )
                drawFace->setStyle( DrawFace::FACE_FILL);
            else
                drawFace->setStyle( DrawFace::FACE_LINES);

            bool cull = cullCheckBox->isChecked();
            drawFace->setBackfaceCull( cull );


                      val  = materialCheckBox->isChecked();
                      drawFace->setMaterial(val);


    val  = shadingCheckBox->isChecked();
    if( val )
        drawFace->setShade(JRenderMode::SMOOTH_SHADE);
    else
        drawFace->setShade(JRenderMode::FLAT_SHADE);
    */

    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: reverseAll()
{
    if( mesh == nullptr) return;
    mesh->getTopology()->reverseAll();
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshFacesDialog :: deleteMesh()
{
    if( mesh == nullptr ) return;

    /*
         if( deleteFacesKeepNodesRadioButton->isChecked() ) {
              mesh->deleteFaces( MeshEntity::ANY_ENTITY);
         }

         if( deleteFacesKeepBoundaryRadioButton->isChecked() ) {
              mesh->deleteFaces( MeshEntity::INTERNAL_ENTITY);
         }

         if( deleteAllRadioButton->isChecked() )  {
              mesh->deleteAll();
         }
    */

    size_t numfaces = mesh->getSize(2);
    numFacesLineEdit->setText( QString::number(numfaces) );

    meshViewer->refreshDisplay();

}
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: getBoundary()
{
    if( mesh == nullptr ) return;

    mesh->getTopology()->searchBoundary();

    size_t numfaces = mesh->getTopology()->getBoundarySize(2);
    numBoundaryLineEdit->setText( QString::number(numfaces) );

    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void JMeshFacesDialog :: openNormalsDialog()
{
    if( normalsDialog.get() == nullptr )
        normalsDialog.reset(new JMeshNormalsDialog( this ));

    normalsDialog->setViewManager( viewManager );
    normalsDialog->setMesh(mesh, 2);
    normalsDialog->show();
    this->hide();
}

////////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: openCongruentDialog()
{
    if( congruentSetDialog.get() == nullptr )
        congruentSetDialog.reset(new JCongruentSetDialog( this ));

    congruentSetDialog->setViewManager( viewManager );
    congruentSetDialog->show();
    this->hide();
}
*/

////////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: makeConnections()
{
    PushButton( attribListPushButton,  [=] {openAttribListDialog();});
    PushButton( getBoundaryPushButton, [=] {getBoundary();});
    PushButton( normalsPushButton,     [=] {openNormalsDialog(); });
    PushButton( internalPushButton,    [=] {setInternal(); });
    PushButton( boundaryPushButton,    [=] {setBoundary(); });
    PushButton( interfacePushButton,   [=] {setInterface(); });
    PushButton( lookAtPushButton,      [=] {lookAt(); });
    PushButton( reverseAllPushButton,  [=] {reverseAll(); });

    CheckBox( enumCheckBox, [=] {checkState();});
    CheckBox( cullCheckBox,  [=] {checkState();});
    CheckBox( filledCheckBox, [=] { checkState();});
    CheckBox( lightsCheckBox, [=] { checkLights();});
    CheckBox( shadingCheckBox, [=] { checkState();});
    CheckBox( displayCheckBox, [=] {checkDisplay();});
    CheckBox( displayQuadsCheckBox, [=] {displayFaceType();});
    CheckBox( displayPolygonsCheckBox, [=] {displayFaceType();});
    CheckBox( displayTrianglesCheckBox, [=] {displayFaceType(); });

    RadioButton( showAllRadioButton,      [=] { displayFaces();});
    RadioButton( showConvexRadioButton,   [=] { displayFaces(); });
    RadioButton( showInvertedRadioButton, [=] { displayFaces();});
    RadioButton( showConcaveRadioButton,  [=] { displayFaces();});

    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: buildEdges()
{
     if( mesh == nullptr ) return;
     mesh->getTopology()->collect_edges();

     meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshFacesDialog :: buildDualGraph()
{
    if( mesh == nullptr ) return;

    Mesh *dgraph = mesh->getTopology()->getDualGraph();
    meshViewer->setDualGraph( dgraph );

    bool val = 1;
    size_t nSize = dgraph->getSize(0);
    for( size_t i = 0; i < nSize; i++) {
         Vertex *vtx = dgraph->getNodeAt(i);
         vtx->setAttribute("Display", val);
    }

    nSize = dgraph->getSize(1);
    for( size_t i = 0; i < nSize; i++) {
         Edge *edge = dgraph->getEdgeAt(i);
         edge->setAttribute("Display", val);
    }
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshFacesDialog :: getConsistent()
{
     if( mesh == nullptr ) return;
     mesh->getTopology()->getConsistent();
     meshViewer->refreshDisplay();
}

void JMeshFacesDialog :: activate( Face *face)
{
     if( face == nullptr ) return;

     bool val = 0;
     if( currFace )
          currFace->setAttribute("Display", val);

     val = 1;
     currFace = face;
     currFace->setAttribute("Display", val);

     if( changeCenterCheckBox->isChecked() ) {
          changeCenter(0);
     }

     if( alignCheckBox->isChecked() ) {
          alignAlongAxis(0);
     }

     meshViewer->refreshDisplay();
}
*/
///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: showNeighs()
{
     if( mesh == nullptr ) return;

     if( showOneCheckBox->isChecked() )  {
          int val = nodeNeighsCheckBox->isChecked() +
                    edgeNeighsCheckBox->isChecked() +
                    faceNeighsCheckBox->isChecked() +
                    cellNeighsCheckBox->isChecked();

          if( val ) {
               meshViewer->displayAll(0);
               reset_neighs(1);
               activate( currFace );
          }
     } else {
//        meshViewer->displayAll(2,0);
          meshViewer->displayAll(1);
     }

     meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: reset_neighs(bool val )
{
     if( mesh == nullptr ) return;
     if( currFace == nullptr ) return;

     if( cellNeighsCheckBox->isChecked() ) {
          if( mesh->getAdjTable(2,3) == 0) mesh->buildRelations(2,3);

          meshViewer->displayAll( 3, 1);
          JCellSequence cneighs;
          currFace->getRelations(cneighs);
          for( size_t i = 0; i < cneighs.size(); i++) {
               Cell *cell = cneighs[i];
               cell->setAttribute("Display", val);
          }

          if( faceNeighsCheckBox->isChecked() ) {
               meshViewer->displayAll( 2, 1);
               for( size_t i = 0; i < cneighs.size(); i++) {
                    Cell *cell = cneighs[i];
                    for( int j = 0; j < cell->getSize(2); j++) {
                         Face *face = cell->getFaceAt(j);
                         face->setAttribute("Display", val);
                    }
               }
          }

          if( edgeNeighsCheckBox->isChecked() ) {
               meshViewer->displayAll( 1, 1);
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

     if( edgeNeighsCheckBox->isChecked() ) {
          meshViewer->displayAll(1, 1);
          for( int  i = 0; i < currFace->getSize(1); i++) {
               Edge *edge = currFace->getEdgeAt(i);
               edge->setAttribute("Display", val);
          }
     }

     if( nodeNeighsCheckBox->isChecked() ) {
          meshViewer->displayAll(0, 1);
          for( int i = 0; i < currFace->getSize(0); i++) {
               Vertex *vtx = currFace->getNodeAt(i);
               vtx->setAttribute("Display", val);
          }
     }
     meshViewer->refreshDisplay();
}
*/

/////////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: nextFace()
{
     if( mesh == nullptr ) return;
     if( !showOneCheckBox->isChecked() ) return;

     reset_neighs(0);

     size_t id  = faceIDSlider->value() + 1;

     if( id >= mesh->getSize(2) ) return;

     faceIDSlider->setValue(id);
     Face *face =  mesh->getFaceAt( id ) ;
     activate( face );

     reset_neighs(1);

     meshViewer->refreshDisplay();
}
*/
///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: showOne()
{
     if( mesh == nullptr ) return;

     if( showOneCheckBox->isChecked() ) {
          meshViewer->displayAll(2, 0);
          meshViewer->displayAll(3, 0);
          size_t currFaceID = faceIDSlider->value();
          if( currFaceID >= mesh->getSize(2) ) return;
          currFace =  mesh->getFaceAt( currFaceID ) ;
          if( currFace ) {
               bool val = 1;
               currFace->setAttribute("Display", val);
//             reset_neighs(1);
          }
     } else {
          meshViewer->displayAll(2, 1);
     }

     setNumVisible();
     meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////
/*
void JMeshFacesDialog :: getFaceID()
{
     if( mesh == nullptr ) return;

     size_t currFaceID  = faceIDSlider->value();
     Face *face =  mesh->getFaceAt( currFaceID ) ;
     activate(face );
}

void JMeshFacesDialog :: changeCenter( bool refresh )
{
     const Point3D &xyz = currFace->getNodeAt(0)->getXYZCoords();
//        Vertex::mid_point(currEdge->getNodeAt(0), currEdge->getNodeAt(1), xyz);
     AffineTransform af(mesh);
     af.translate(-xyz[0], -xyz[1], -xyz[2]);

     if( refresh ) meshViewer->refreshDisplay();
}

void JMeshFacesDialog :: checkLower()
{
    bool val;

    val = lowerNodesCheckBox->isChecked();
    drawFace->display_lower_entity(0, val);

    val = lowerEdgesCheckBox->isChecked();
    drawFace->display_lower_entity(1, val);
    meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: setUniformColor()
{
     if( uniformColorCheckBox->isChecked() ) {
          uniform_color = 1;
          Color clr = faceColor->getInternalColor();
          faceColor->setBoundaryColor( clr );
     } else
          uniform_color = 0;
}
*/
///////////////////////////////////////////////////////////////////////////////

/*
void JMeshFacesDialog :: setDefaultColorMethod()
{
if( meshViewer == nullptr ) return;
meshViewer->displayAll(2,1);
     DrawFace *drawFace = meshViewer->getDrawFace();
     FaceColor *nColor  = drawFace->getDefaultColorMethod();

     string str = nColor->getName();
//   colorMethodLineEdit->setText( QString::fromStdString(str));
     drawFace->setColorMethod( nColor );
     meshViewer->refreshDisplay();

}
*/
#ifdef DFDF
void JMeshFacesDialog :: alignAlongAxis( bool refresh )
{
    if( !alignCheckBox->isChecked()  ) return;

    meshViewer->alignAlong(currFace, 2, 0);
    if( refresh ) meshViewer->refreshDisplay();

    /*
         if( currFace == nullptr ) return;

         const Point3D &p1 = currFace->getNodeAt(0)->getXYZCoords();
         const Point3D &p2 = currFace->getNodeAt(1)->getXYZCoords();
         const Point3D &p3 = currFace->getNodeAt(2)->getXYZCoords();

         AffineTransform af(mesh);
         af.translate(-p1[0], -p1[1], -p1[2] );

         Vec3D xAxis, yAxis;
         xAxis[0] = 1.0;
         xAxis[1] = 0.0;
         xAxis[2] = 0.0;

         yAxis[0] = 0.0;
         yAxis[1] = 1.0;
         yAxis[2] = 0.0;

         Point3D edgevec;
         edgevec[0] = p2[0] - p1[0];
         edgevec[1] = p2[1] - p1[1];
         edgevec[2] = p2[2] - p1[2];

         Vec3D  perpAxis;
         qglviewer::Vec prot;

         JMath::cross_product( xAxis, edgevec, perpAxis);
         double angle = JMath::getVecAngle(edgevec, xAxis, ANGLE_IN_RADIANS);
         qglviewer::Vec axis1(perpAxis[0], perpAxis[1], perpAxis[2] );
         qglviewer::Quaternion quaternion1(axis1, -1.0*angle);

         edgevec[0] = p3[0] - p2[0];
         edgevec[1] = p3[1] - p2[1];
         edgevec[2] = p3[2] - p2[2];

         angle = JMath::getVecAngle(edgevec, yAxis, ANGLE_IN_RADIANS);
         qglviewer::Vec axis2( 1.0, 0.0, 0.0 );
         qglviewer::Quaternion quaternion2(axis2, -1.0*angle);

    //     qglviewer::Quaternion quaternion = quaternion2*quaternion1;
         qglviewer::Quaternion quaternion = quaternion1;

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

         if( refresh ) meshViewer->refreshDisplay();
    */
}
#endif

