#include "MeshSurfaceVectorFieldDialog.hpp"

//////////////////////////////////////////////////////////////////////////////////////
void JVectorFieldViewer :: sampleField( const JFacePtr &face)
{
    Point3D pc;
    face->getAvgXYZ(pc);
    const Point3D & p0 = face->getNodeAt(0)->getXYZCoords();
    double len = JMath::length(p0, pc);

    double vlen = len;

    size_t faceid = face->getID();
    const JEdgePtr &fieldDir = vecField->getEdgeAt(faceid);

    const Point3D &tail = fieldDir->getNodeAt(0)->getXYZCoords();
    const Point3D &head = fieldDir->getNodeAt(1)->getXYZCoords();

    double dx = head[0] - tail[0];
    double dy = head[1] - tail[1];
    double dz = head[2] - tail[2];
    double dl = sqrt(dx*dx + dy*dy + dz*dz);

    double vmin = 0.0;
    double vmax = 1.0;

    const Point3D &p1 = face->getNodeAt(0)->getXYZCoords();
    const Point3D &p2 = face->getNodeAt(1)->getXYZCoords();
    const Point3D &p3 = face->getNodeAt(2)->getXYZCoords();

    FieldVec fv;
    Point3D  pn;
    for( int i = 0; i < numSamples; i++) {
        double a =   JMath::random_value(vmin, vmax);
        double b =   JMath::random_value(vmin, 1.0-a);
        double c =   1.0 - a - b;
        assert( a + b + c == 1.0);
        fv.tail[0] = a*p1[0] + b*p2[0] + c*p3[0];
        fv.tail[1] = a*p1[1] + b*p2[1] + c*p3[1];
        fv.tail[2] = a*p1[2] + b*p2[2] + c*p3[2];
        double  r  = JMath::random_value(0.2, 1.0);
        fv.head[0] = fv.tail[0] + r*vlen*dx/dl;
        fv.head[1] = fv.tail[1] + r*vlen*dy/dl;
        fv.head[2] = fv.tail[2] + r*vlen*dz/dl;
        fieldVec.push_back(fv);
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void JVectorFieldViewer :: genField()
{
    fieldVec.clear();
    if( mesh == nullptr || vecField == nullptr) return;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) sampleField(face);
    }
}

//////////////////////////////////////////////////////////////////////////////////////
void JVectorFieldViewer :: draw()
{
    glEnable( GL_BLEND);
    glDisable( GL_LIGHTING);

    size_t n = fieldVec.size();

    for( int i  = 0; i < n; i++) {
        glBegin(GL_LINES);
        glColor3f( 0.0, 0.0, 0.0);
        glVertex3fv( fieldVec[i].tail);
        glColor3f( 0.9, 0.9, 0.9);
        glVertex3fv( fieldVec[i].head);
        glEnd();
    }
    glEnable( GL_LIGHTING);
    glDisable( GL_BLEND);
}
//////////////////////////////////////////////////////////////////////////////////////
JMeshSurfaceVectorFieldDialog :: JMeshSurfaceVectorFieldDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    length      = 1.0;
    numRandomVectorsLineEdit->setText(QString(QString::number(1)));
    currFieldGenerator = 0;
}

///////////////////////////////////////////////////////////////////////////////

JMeshSurfaceVectorFieldDialog :: ~JMeshSurfaceVectorFieldDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSurfaceVectorFieldDialog :: init()
{
    if( viewManager == nullptr ) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    if( c == nullptr) return;

    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    picker = meshViewer->getEntityPicker();
    if( picker ) picker->setMode(1);

    viewManager->attach( this );
    viewManager->setMouseTracking(1);
    viewManager->setSelectRegionHeight(10);
    viewManager->setSelectRegionWidth(10);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSurfaceVectorFieldDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string  name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 2;

    /*
        constraintMesh = JMesh::newObject();
        meshViewer->addObject(constraintMesh);
        meshViewer->setCurrentMesh(mesh); // You need to pick faces this mesh...
    */

    picker->setMesh(mesh);
    mesh->deleteFaceAttribute("ConstraintVector");
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSurfaceVectorFieldDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }

    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( picker == nullptr ) return;

    /*
        JFaceSequence  faces = picker->getPickedFaces();
        if( faces.empty() ) return;
        pickedFaces.push_back( faces[0] );

        currPickedFace = faces[0];

        JColor highlight;
        Point3D xyz, p0, p1;
        if( !currPickedFace->hasAttribute("ConstraintVector") )  {
            JNodePtr v0 = JNode::newObject();
            JNodePtr v1 = JNode::newObject();
            currPickedFace->getAvgXYZ(xyz);
            v0->setXYZCoords(xyz);
            p0 = currPickedFace->getNodeAt(0)->getXYZCoords();
            p1 = currPickedFace->getNodeAt(1)->getXYZCoords();
            xyz[0] += p1[0] - p0[0];
            xyz[1] += p1[1] - p0[1];
            xyz[2] += p1[2] - p0[2];
            v1->setXYZCoords(xyz);
            JEdgePtr edge = JEdge::newObject(v0,v1);
            constrainedMesh->addObject(v0);
            constrainedMesh->addObject(v1);
            constrainedMesh->addObject(edge);
            currPickedFace->setAttribute("ConstraintVector", edge);
            meshViewer->updateBuffers(constrainedMesh);
            JEdgeRenderPtr eAttrib;
            edge->getAttribute("Render", eAttrib);
            highlight[0] = 0.2;
            highlight[1] = 0.2;
            highlight[2] = 1.0;
            highlight[3] = 1.0;
            eAttrib->color = highlight;
            eAttrib->scale = 1.5;
        }

        JFaceRenderPtr fAttrib;
        highlight[0] = 1.0;
        highlight[1] = 0.2;
        highlight[2] = 0.2;
        highlight[3] = 1.0;
        faces[0]->getAttribute("Render", fAttrib);
        fAttrib->color = highlight;
        meshViewer->updateBuffers(mesh);
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: getVecField()
{
    if( mesh == nullptr) return;

    if( mesh->getTopology()->getDimension() != 2) {
        cout << "Warning: Vector fields only for surfaces " << endl;
        return;
    }

    if( mesh->getTopology()->isClosed() ) {
        size_t nCount = mesh->getNumAttributes("ConstraintVector", 2);
        if( nCount == 0  && !randomConstraintsCheckBox->isChecked() ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Warning At least one of the face must have constraint vector: Select one face");
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) return;
        }
    }

    for(const JMeshPtr &m: vecFields)
        meshViewer->removeObject(m);

    JWaitCursor waitCursor;
    waitCursor.start();
    if( nrosyFieldRadioButton->isChecked() ) getNRoSyField();
    if( polyvectorsRadioButton->isChecked()) getPolyVecField();
    if( conjugateFieldRadioButton->isChecked()) getConjugateField();
    if( integrableFieldRadioButton->isChecked()) getIntegrableField();

    vecComponentSpinBox->setMaximum( vecFields.size()-1);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void JMeshSurfaceVectorFieldDialog :: getNRoSyField()
{

    JNRoSyField *nRoSyField = new JNRoSyField;
    surfVecFieldPtr.reset( nRoSyField );

    nRoSyField->setMesh(mesh);
    int nd = numVecPerFaceSpinBox->value();
    nRoSyField->setNumVectorsPerFace(nd);

    if( randomConstraintsCheckBox->isChecked() ) {
        QString qstr = numRandomVectorsLineEdit->text();
        int n = qstr.toInt();
        nRoSyField->genRandomConstraints(n);
    }

    if( readConstraintsCheckBox->isChecked() ) {
        nRoSyField->readConstraints( constraintFiles[0], constraintFiles[1] );
    }

    int err = nRoSyField->genField();

    if( err ) return;

    vecFields  = nRoSyField->getVecFields();

    int index = 0;
    for( const JMeshPtr &m: vecFields) {
        vecFields[index]->setName("NRosyField" + to_string(index));
        meshViewer->addObject(vecFields[index]);
        index++;
    }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: getPolyVecField()
{
    JPolyVectorsField *polyVecField = new JPolyVectorsField;
    surfVecFieldPtr.reset( polyVecField );

    int nd = numVecPerFaceSpinBox->value();
    polyVecField->setNumVectorsPerFace(nd);

    polyVecField->setMesh(mesh );

    if( randomConstraintsCheckBox->isChecked() ) {
        QString qstr = numRandomVectorsLineEdit->text();
        int n = qstr.toInt();
        polyVecField->genRandomConstraints(n);
    }

    if( readConstraintsCheckBox->isChecked() ) {
        polyVecField->readConstraints( constraintFiles[0], constraintFiles[1] );
    }

    int err = polyVecField->genField();
    if( err ) return;

    vecFields  = polyVecField->getVecFields();

    int index = 0;
    for( const JMeshPtr &m: vecFields) {
        vecFields[index]->setName("PolyVecField" + to_string(index));
        meshViewer->addObject(vecFields[index]);
        index++;
    }

//  displayField( polyVecField );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSurfaceVectorFieldDialog :: getConjugateField()
{
    /*
        meshViewer->removeObject(vecField);

        if( conjugateFieldPtr == nullptr)
            conjugateFieldPtr.reset( new JConjugateField );

        conjugateFieldPtr->setMesh(mesh );

        if( randomConstraintsCheckBox->isChecked() ) {
            QString qstr = numRandomVectorsLineEdit->text();
            int n = qstr.toInt();
            conjugateFieldPtr->genRandomConstraints(n);
        }
        if( readConstraintsCheckBox->isChecked() ) {
            conjugateFieldPtr->readConstraints( constraintFiles[0], constraintFiles[1] );
        }

        int err = conjugateFieldPtr->genField();
        if( err ) return;

        vecField  = conjugateFieldPtr->getVecField();
        vecField->setName("ConjugateField");
        meshViewer->addObject(vecField);

        displayField( conjugateFieldPtr.get() );
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: getIntegrableField()
{
    /*
        meshViewer->removeObject(vecField);

        JIntergrableField *integrableField = new JIntegrableField;

        integrableField->setMesh(mesh );

        if( randomConstraintsCheckBox->isChecked() ) {
            QString qstr = numRandomVectorsLineEdit->text();
            int n = qstr.toInt();
            integrableField->genRandomConstraints(n);
        }
        if( readConstraintsCheckBox->isChecked() ) {
            integrableField->readConstraints( constraintFiles[0], constraintFiles[1] );
        }

        int err = integrableField->genField();
        if( err ) return;

        vecField  = integrableField->getVecField();
        vecField->setName("ConjugateField");
        meshViewer->addObject(vecField);
        displayField( integrableField );
    */
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSurfaceVectorFieldDialog :: openEdgeAttribsDialog()
{
    if( edgeAttribsDialog == nullptr)
        edgeAttribsDialog.reset( new JEdgeAttributesDialog(this));

    edgeAttribsDialog->setViewManager( viewManager );
    edgeAttribsDialog->show();
    this->hide();

    if( vecFields.empty() ) return;
    int id = vecComponentSpinBox->value();
    edgeAttribsDialog->setMesh(vecFields[id]);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: clearField()
{
    for( const JMeshPtr &m : vecFields)
        meshViewer->removeObject(m);
    vecFields.clear();
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: readFixedFacesID()
{
    /*
        QString qstr = QFileDialog::getOpenFileName(this,
                       *new QString("Select Mesh File "),
                       lastSelectedDirectory, *new QString( "esh Format (*.mat)"));

        constraintFiles[0] = qstr.toUtf8().constData();
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: readFixedFacesVectors()
{
    /*
        QString qstr = QFileDialog::getOpenFileName(this,
                       *new QString("Select Mesh File "),
                       lastSelectedDirectory, *new QString( "esh Format (*.mat)"));

        constraintFiles[1] = qstr.toUtf8().constData();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSurfaceVectorFieldDialog :: setConstraints()
{
    if( mesh == nullptr) return;
    mesh->deleteFaceAttribute("ConstraintVector");

    Array3D val;
    size_t  numfaces = mesh->getSize(2);

    int nSamples;
    double minCutOff, maxCutOff;

    JFaceSet fixedFaces;
    vector<double> lowVal, highVal;
    vector<double> vec;
    if( minCurvatureCheckBox->isChecked() || maxCurvatureCheckBox->isChecked() ) {
        JMeshCurvature mcurv;
        mcurv.setMesh(mesh);
        mcurv.setEvalPos(2);
        mcurv.setMeanCurvature();
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            f->getAttribute("Curvature", val);
            lowVal.push_back( val[1] );
            highVal.push_back( val[2] );
        }

        if( maxCurvatureCheckBox->isChecked() )  {
            sort( highVal.begin(), highVal.end() );
            nSamples = 0.01*maxCurvatureSpinBox->value()*numfaces;
            maxCutOff = highVal[numfaces - 1 - nSamples];
            size_t index = 0;
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &f = mesh->getFaceAt(i);
                f->getAttribute("Curvature", val);
                if( val[2] > maxCutOff) {
                    f->getAttribute("MaxCurvatureDir", vec);
                    f->setAttribute("ConstraintVector", vec);
                    fixedFaces.insert(f);
                    index++;
                }
                if( index  > nSamples) break;

            }
        }

        if( minCurvatureCheckBox->isChecked() )  {
            sort( lowVal.begin(), lowVal.end() );
            nSamples  = 0.01*minCurvatureSpinBox->value()*numfaces;
            minCutOff = lowVal[numfaces - 1 - nSamples];
            size_t index = 0;
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &f = mesh->getFaceAt(i);
                f->getAttribute("Curvature", val);
                if( val[1] > minCutOff) {
                    f->getAttribute("MinCurvatureDir", vec);
                    f->setAttribute("ConstraintVector", vec);
                    fixedFaces.insert(f);
                    index++;
                }
                if( index > nSamples) break;
            }
        }
    }

    JFaceRenderPtr fAttrib;
    JColor color = JEntityColor::getColor("Red");

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->hasAttribute("ConstraintVector") ) {
            face->getAttribute("Render", fAttrib);
            fAttrib->color  = color;
            fAttrib->display = 1;
        }
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: closeDialog()
{
    this->close();
    parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: showAllVectors()
{
    JMeshRenderPtr mrender;
    for( const JMeshPtr &mesh : vecFields) {
        mesh->getAttribute("Render", mrender);
        mrender->display = 1;
    }
    meshViewer->refreshDisplay();

}
///////////////////////////////////////////////////////////////////////////////
void JMeshSurfaceVectorFieldDialog :: showVecComponent()
{
    if( surfVecFieldPtr == nullptr) return;

    if( vecFields.empty() ) return;


    JMeshRenderPtr mrender;
    for( const JMeshPtr &mesh : vecFields) {
        mesh->getAttribute("Render", mrender);
        mrender->display = 0;
    }
    int id = vecComponentSpinBox->value();
    if( id < vecFields.size() ) {
        vecFields[id]->getAttribute("Render", mrender);
        mrender->display = 1;
    }
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSurfaceVectorFieldDialog :: makeConnections()
{
    PushButton( applyPushButton,  [=] {getVecField();});

    PushButton( clearPushButton, [=] {clearField();});
    PushButton( constraintsPushButton, [=] {setConstraints(); });
    PushButton( edgeAttribsPushButton, [=] {openEdgeAttribsDialog();});
    PushButton( readFixedFacesIDPushButton, [=] {readFixedFacesID();});
    PushButton( readFixedFacesVectorsPushButton, [=] {readFixedFacesVectors(); });
    PushButton( showAllPushButton, [=] { showAllVectors(); });
    PushButton( showComponentPushButton, [=] { showVecComponent(); });

    PushButton( closePushButton, [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////

/*
void JSurfaceVectorFieldDialog :: rotateVector()
{
    if( currPickedFace == nullptr) return;

    Vec3F normal;
    currPickedFace->getAttribute("Normal", normal);

    double angle = M_PI*rotateVectorSpinBox->value()/180.0;

    qglviewer::Vec axis(normal[0], normal[1], normal[2] );
    qglviewer::Quaternion quaternion(axis, -1.0*angle);

    const Point3D p0 = currPickedFace->getNodeAt(0)->getXYZCoords();
    const Point3D p1 = currPickedFace->getNodeAt(1)->getXYZCoords();

    double dx = p1[0] - p0[0];
    double dy = p1[1] - p0[1];
    double dz = p1[2] - p0[2];
    double dl = sqrt(dx*dx + dy*dy + dz*dz);

    qglviewer::Vec edgeVec(dx/dl, dy/dl, dz/dl);

    qglviewer::Vec prot = quaternion.rotate(edgeVec);

    JEdgePtr edge;
    int err = currPickedFace->getAttribute("ConstraintVector", edge);
    if( err) return;
    Point3D xyz;
    currPickedFace->getAvgXYZ(xyz);

    xyz[0] += length*prot[0];
    xyz[1] += length*prot[1];
    xyz[2] += length*prot[2];
    edge->getNodeAt(1)->setXYZCoords(xyz);
    meshViewer->updateBuffers(constraintMesh);
}
void JSurfaceVectorFieldDialog :: openNodeAttribsDialog()
{
    if( vecField == nullptr) return;

    if( nodeAttribsDialog == nullptr)
        nodeAttribsDialog.reset( new JNodeAttributesDialog(this));

    nodeAttribsDialog->setViewManager(viewManager );
    nodeAttribsDialog->setMesh(mesh);
    nodeAttribsDialog->setNodes(singularNodes);
    nodeAttribsDialog->show();
}
*/
void JMeshSurfaceVectorFieldDialog :: displayField(JMeshSurfaceVectorField *field)
{
    /*
        int numSingular;
        numSingular = fieldPtr->getNumSingularities();
        numSingularitiesLineEdit->setText(QString(QString::number(numSingular)));

        int numFixed;
        numFixed = fieldPtr->getNumConstraints();
        numConstraintsLineEdit->setText(QString(QString::number(numFixed)));

        JColor color;
        JNodeRenderPtr nAttrib;
        if( displaySingularitiesCheckBox->isChecked() ) {
            singularNodes = fieldPtr->getPositiveSingularNodes();
            color = JEntityColor::getColor("Red");
            for( const JNodePtr &vtx : singularNodes) {
                vtx->getAttribute("Render", nAttrib);
                nAttrib->glyph = 1;
                nAttrib->display = 1;
                nAttrib->color   = color;
            }
            color = JEntityColor::getColor("Green");
            singularNodes = fieldPtr->getNegativeSingularNodes();
            for( const JNodePtr &vtx : singularNodes) {
                vtx->getAttribute("Render", nAttrib);
                nAttrib->glyph   = 1;
                nAttrib->display = 1;
                nAttrib->color   = color;
            }
        }

        JFaceRenderPtr fAttrib;
        if( displayConstrainedFacesCheckBox->isChecked() ) {
            size_t numfaces = mesh->getSize(2);
            color = JEntityColor::getColor("Gray");
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &f = mesh->getFaceAt(i);
                if( f->isActive() ) {
                    f->getAttribute("Render", fAttrib);
                    fAttrib->color = color;
                }
            }
            color = JEntityColor::getColor("Red");
            constrainedFaces = fieldPtr->getConstrainedFaces();
            for( const JFacePtr &face : constrainedFaces) {
                face->getAttribute("Render", fAttrib);
                fAttrib->display = 1;
                fAttrib->color   = color;
            }
        }
        meshViewer->updateBuffers(mesh);

        vecFieldViewer->setMesh(mesh);
        vecFieldViewer->setField(vecField);
        vecFieldViewer->setNumSamplesPerFace(5);
        vecFieldViewer->genField();
    */
}


void JMeshSurfaceVectorFieldDialog :: setVecLength()
{
    /*
        if( currPickedFace == nullptr) return;

        JEdgePtr edge;
        int err = currPickedFace->getAttribute("ConstraintVector", edge);
        if( err) return;

        const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
        const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();
        double dx = p1[0] - p0[0];
        double dy = p1[1] - p0[1];
        double dz = p1[2] - p0[2];
        double dl = sqrt(dx*dx + dy*dy + dz*dz);

        QString qstr = scaleLineEdit->text();
        double len   = qstr.toDouble();

        len = max(len, 1.0E-06);

        Point3D pnew;
        pnew[0] = p0[0] + len*dx;
        pnew[1] = p0[1] + len*dy;
        pnew[2] = p0[2] + len*dz;
        edge->getNodeAt(1)->setXYZCoords(pnew);
        meshViewer->updateBuffers(constraintMesh);
    */
}

