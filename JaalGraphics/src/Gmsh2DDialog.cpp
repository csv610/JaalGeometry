#include "Gmsh2DDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JGmsh2DDialog :: JGmsh2DDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    lengthFactorSpinBox->setValue(lengthFactor);
    maxAnisotrophyLineEdit->setText( QString::number(maxAnisotrophy));
    maxElementSizeLineEdit->setText( QString::number(maxElementSize));
    minElementSizeLineEdit->setText( QString::number(minElementSize));
}

///////////////////////////////////////////////////////////////////////////////

JGmsh2DDialog :: ~JGmsh2DDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JGmsh2DDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JGmsh2DDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    JNodeSequence nodes;
    mesh->getTopology()->getBoundary(nodes);

    for( const JNodePtr &vtx : nodes) {
        const Point3D &xyz = vtx->getXYZCoords();
        vtx->setAttribute("TargetPos", xyz);
    }
    modifyBoundary = 0;
}

///////////////////////////////////////////////////////////////////////////////
void JGmsh2DDialog :: showEvent( QShowEvent *event)
{
     QDialog::showEvent(event);
}
///////////////////////////////////////////////////////////////////////////////

void JGmsh2DDialog :: recoverBoundary()
{
    if( newQuadMesh == nullptr || mesh == nullptr || modifyBoundary == 0) return;

    JNodeSequence srcnodes, dstnodes;
    mesh->getTopology()->getBoundary(srcnodes);
    newQuadMesh->getTopology()->getBoundary(dstnodes);

    /*
       if( srcnodes.size() != dstnodes.size() ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("#Boundary nodes on source and destination objects differ" ;
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) return;
       }
    */
}
///////////////////////////////////////////////////////////////////////////////
void JGmsh2DDialog :: getOptions()
{
    QString qstr;
    string  str;

    algorithm = 0;
    qstr = algorithmComboBox->currentText();
    str  = StdString(qstr);
    if( str == "MeshAdapt")      algorithm = 1;
    if( str == "Automatic")      algorithm = 2;
    if( str == "Delaunay")       algorithm = 5;
    if( str == "Frontal")        algorithm = 6;
    if( str == "DelaunayQuads")  algorithm = 8;
    if( str == "PackingQuads")   algorithm = 9;

    recombinationAlgorithm = 0;
    if( blossomRadioButton->isChecked()  )
        recombinationAlgorithm = 1;

    recombineAll = 0;
    if( recombineAllCheckBox->isChecked()  )
        recombineAll = 1;

    smoothSteps = smoothStepsSpinBox->value();

    lloydSteps = lloydStepsSpinBox->value();

    crossFieldSmoothSteps = crossFieldSmoothStepsSpinBox->value();

    sizeFromPoints    = sizeFromPointsCheckBox->isChecked();

    sizeFromCurvature = sizeFromCurvatureCheckBox->isChecked();

    lengthFactor      = lengthFactorSpinBox->value();

    qstr = minElementSizeLineEdit->text();
    minElementSize = qstr.toDouble();

    qstr = maxElementSizeLineEdit->text();
    maxElementSize = qstr.toDouble();
}

///////////////////////////////////////////////////////////////////////////////

void JGmsh2DDialog :: reparamBoundary()
{
    if( mesh == nullptr) return;

    vector<JEdgeSequence> boundCurves;

    mesh->getTopology()->getBoundary(boundCurves);
    int ncurves = boundCurves.size();

    for( int i = 0; i < ncurves; i++)
        JEdgeGeometry::makeUniform( boundCurves[i] );

    meshViewer->updateBuffers(mesh);
    modifyBoundary = 1;
}

///////////////////////////////////////////////////////////////////////////////

void JGmsh2DDialog :: smoothBoundary()
{
    if( mesh == nullptr) return;

    int niter = smoothCurvesSpinBox->value();

    vector<JEdgeSequence> boundCurves;

    mesh->getTopology()->getBoundary(boundCurves);
    int ncurves = boundCurves.size();

    for( int i = 0; i < ncurves; i++)
        JEdgeGeometry::smooth( boundCurves[i], niter);
    meshViewer->updateBuffers(mesh);

    modifyBoundary = 1;
}
///////////////////////////////////////////////////////////////////////////////

void JGmsh2DDialog :: generate()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    getOptions();

    JMeshPtr halfmesh = mesh;

    vector<Point3D> midPoints;
    JNodeSequence skippedNodes;
    Point3D pmid;

    if( halfSamplesCheckBox->isChecked() ) {
        if( recombineAll) {
            halfmesh = JMesh::newObject();

            vector<JEdgeSequence> boundloops;
            mesh->getTopology()->getBoundary(boundloops);

            JNodeSequence boundnodes;
            int numloops = boundloops.size();
            int indx = 0;

            for( int i = 0; i < numloops; i++) {
                JEdgeTopology::getChainNodes(boundloops[i], boundnodes);
                int nnodes = boundnodes.size();
                if( nnodes%2 ) {
                    QMessageBox msg;
                    msg.setIcon(QMessageBox::Critical);
                    msg.setText("#Boundary nodes is not even");
                    msg.setStandardButtons(QMessageBox::Ok);
                    int ret = msg.exec();
                    if( ret == QMessageBox::Ok ) return;
                }
                assert( nnodes%2 == 0);
                for( int j = 0; j < nnodes/2; j++) {
                    halfmesh->addObject(boundnodes[2*j] );
                    skippedNodes.push_back( boundnodes[2*j+1]);
                    boundnodes[2*j]->setID(indx++);
                }
                for( int j = 0; j < nnodes/2; j++) {
                    const JNodePtr &v1 = boundnodes[2*j];
                    const JNodePtr &v2 = boundnodes[(2*j+2)%nnodes];
                    pmid = JNodeGeometry::getMidPoint(v1,v2);
                    midPoints.push_back(pmid);
                    JEdgePtr edge = JEdge::newObject(v1,v2);
                    halfmesh->addObject( edge );
                }
            }
        }
    }

    // First write the geometry into "geo" format"
    JMeshGmshExporter  mexp;
    mexp.writeFile(halfmesh, "model.geo");

    ofstream ofile("model.geo", fstream::out | fstream::app);
    if( !ofile.fail() ) {
        ofile << "Mesh.Algorithm = " << algorithm << ";" << endl;
        ofile << "Mesh.RecombinationAlgorithm =  " << recombinationAlgorithm << ";" << endl;
        ofile << "Mesh.RecombineAll = " << recombineAll << ";" << endl;
        ofile << "Mesh.Smoothing = " << smoothSteps << ";" << endl;
        ofile << "Mesh.Lloyd = " << lloydSteps << ";" << endl;
        ofile << "Mesh.SmoothCrossField = " << crossFieldSmoothSteps << ";" << endl;
        ofile << "Mesh.CharacteristicLengthFromPoints = " << sizeFromPoints << ";" << endl;
        ofile << "Mesh.CharacteristicLengthExtendFromBoundary = 1; " << endl;
        ofile << "Mesh.CharacteristicLengthFromCurvature = " << sizeFromCurvature << ";" << endl;
        ofile << "Mesh.AnisoMax = " << maxAnisotrophy << ";" << endl;
//      ofile << "Mesh.Bunin = 1;" << endl;
//      ofile << "Mesh.CharacteristicLengthMin = " << minElementSize << ";" << endl;
//      ofile << "Mesh.CharacteristicLengthMax = " << maxElementSize << ";" << endl;
        ofile << "Mesh.CharacteristicLengthFactor = " << lengthFactor << ";" << endl;
        ofile.close();
    }

    system("gmsh model.geo -2");

    if( newQuadMesh ) meshViewer->removeObject(newQuadMesh);

    // First write the geometry into "geo" format"
    newQuadMesh = JMeshIO::readFile( "model.msh");

    if( newQuadMesh  == nullptr) return;

    if( halfSamplesCheckBox->isChecked() ) {
        if( recombineAll) {
            JNodeSequence newBoundNodes;
            newQuadMesh->getTopology()->getBoundary( newBoundNodes);
            int nmid = midPoints.size();

            for( int i = 0; i < nmid; i++) {
                int found = 0;
                for( const JNodePtr &vtx : newBoundNodes) {
                    const Point3D &p1 = vtx->getXYZCoords();
                    double dist  = JMath::length(midPoints[i], p1);
                    if( dist < 1.0E-10) {
                        const Point3D &pold =  skippedNodes[i]->getXYZCoords();
                        vtx->setXYZCoords( pold );
                        found = 1;
                        break;
                    }
                }
                assert( found );
            }
        }
    }

    newQuadMesh->getTopology()->reverseAll();

    JMeshQuality mq;
    mq.setMesh(newQuadMesh);
    vector<double> faceArea;
    faceArea = mq.getFacesQuality(JMeshQuality::AREA, 0, 1);

    double sumArea   = std::accumulate(faceArea.begin(), faceArea.end(), 0.0);
    double elen      = newQuadMesh->getGeometry()->getMeanEdgeLength();
    int    nexpected = sumArea/(elen*elen);

    expectedQuadsLineEdit->setText( QString::number(nexpected));
    generatedQuadsLineEdit->setText( QString::number(newQuadMesh->getSize(2) ));

    JMSTQuadMesher mstMesher;
    mstMesher.setMesh(newQuadMesh);
    JNodeSequence irrnodes = mstMesher.getSingularNodes();
    numSingularitiesLineEdit->setText( QString::number(irrnodes.size()) );

    minElementSizeLineEdit->setText( QString::number(faceArea.front()));
    maxElementSizeLineEdit->setText( QString::number(faceArea.back()));

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 0;
    mrender->displayEntity[3] = 0;

    meshViewer->addObject(newQuadMesh);
    objectNameLineEdit->setText(QString(newQuadMesh->getName().c_str()));
}

///////////////////////////////////////////////////////////////////////////////
void JGmsh2DDialog :: rejectMesh()
{
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 1;
    mrender->displayEntity[1] = 1;
    mrender->displayEntity[2] = 2;
    mrender->displayEntity[3] = 0;

    if( newQuadMesh ) meshViewer->removeObject(newQuadMesh);
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JGmsh2DDialog :: closeDialog()
{
   this->close();
   parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JGmsh2DDialog :: makeConnections()
{
    connect( reparamBoundaryPushButton, SIGNAL( clicked() ), this, SLOT( reparamBoundary() ));
    connect( smoothBoundaryPushButton, SIGNAL( clicked() ), this, SLOT( smoothBoundary() ));
    connect( generatePushButton, SIGNAL( clicked() ), this, SLOT( generate() ));
    connect( rejectPushButton, SIGNAL( clicked() ), this, SLOT( rejectMesh() ));

    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}

///////////////////////////////////////////////////////////////////////////////



