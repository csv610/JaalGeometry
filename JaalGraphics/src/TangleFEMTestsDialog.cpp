#include "TangleFEMTestsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JTangleFEMTestsDialog :: JTangleFEMTestsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;

    randomTangle = 0;
    checkTangle  = 0;
    meshType     = 3;
    elemOrder    = 1;
    youngModulus = 2.0E+11;
    poissonRatio = 0.25;
    numGauss1    = 3;   // (1,3,7)
    numGauss2    = 1;   // (1,3,7)
    absJacobian  = 0;
    fieldType    = 1;
    shapeFamily  = 0;
    checkConvex  = 1;

    youngModulusLineEdit->setText( QString::number(youngModulus) );
    poissonRatioLineEdit->setText( QString::number(poissonRatio) );
}

///////////////////////////////////////////////////////////////////////////////

JTangleFEMTestsDialog :: ~JTangleFEMTestsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JTangleFEMTestsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    if( c == nullptr) c = JMeshViewer::registerComponent(viewManager);
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    setMesh();
}

///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: analyticalSol( const Point3D &p, double &u, double &v)
{
    double x = p[0];
    double y = p[1];

    if( fieldType == 1 ) {
        u =  0.1*x - 0.2*y + 1.0E-10;
        v = -0.3*x + 0.5*y + 1.0E-10;
    } else {
        double a1  = 10;
        double a2  = 10;
        double nu = 0.3;
        u = a1*x*x - 4.0*a2*x*y/(1+nu);
        v = -4.0*a1*x*y/(1+nu) + a2*y*y;
    }
}
/////////////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: setParams()
{
    QString qstr;
    checkTangle    = inversionCheckBox->isChecked();
    absJacobian    = absJacobianCheckBox->isChecked();
    checkConvex    = convexityCheckBox->isChecked();

    qstr = youngModulusLineEdit->text();
    youngModulus = qstr.toDouble();

    qstr = poissonRatioLineEdit->text();
    poissonRatio = qstr.toDouble();

    if( classicalShapeFuncRadioButton->isChecked() )   shapeFamily = 0;
    if( mixedShapeFuncRadioButton->isChecked() )       shapeFamily = 1;
    if( barycentricShapeFuncRadioButton->isChecked() ) shapeFamily = 2;

    if( linearShapeFuncRadioButton->isChecked() )    elemOrder = 1;
    if( quadraticShapeFuncRadioButton->isChecked() ) elemOrder = 2;

    if( linearFieldRadioButton->isChecked() )    fieldType = 1;
    if( quadraticFieldRadioButton->isChecked() ) fieldType = 2;

    /*
        elastic.setYoungModulus(youngModulus);
        elastic.setPoissonRatio(poissonRatio);
        elastic.checkTangle(checkTangle);
    //  elastic.checkConvexity(checkConvex);
    //  elastic.setShapeFamily(shapeFamily);
    //  elastic.setNumPrimaryFaceGaussPoints(numGauss1);
    //  elastic.setNumSecondaryFaceGaussPoints(numGauss2);

        bool val = absJacobianCheckBox->isChecked();
        elastic.useAbsoluteJacobian(val);
    */

}

///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: setMesh()
{
    /*
        if( meshViewer == nullptr) return;
        if( mesh ) meshViewer->remove(mesh);

        if( trimeshRadioButton->isChecked() )  meshType = 3;
        if( quadmeshRadioButton->isChecked() ) meshType = 4;

        problemID   = testNameComboBox->currentIndex();
        if( problemID == 0) genSingleNodeMesh();
        if( problemID == 1) genSimpleTangleMesh();
        if( problemID == 2) genRandomMesh();
        if( problemID == 3) genRotateQuadMesh();
        if( problemID == 4) genBendingBeamMesh();
        if( problemID == 5) genSquareCircle2Mesh();

        meshViewer->refreshDisplay();
        elastic.setMesh(mesh);
    */
}

///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: patchTest()
{
    /*
        if( mesh == nullptr) return;

        JNodeSequence boundnodes;
        mesh->getTopology()->getBoundary(boundnodes);

        assert( !boundnodes.empty() ) ;
        double uexact, vexact;
        for( int i = 0; i < boundnodes.size(); i++) {
            const Point3D &xyz =  boundnodes[i]->getXYZCoords();
            analyticalSol(xyz, uexact,vexact);
            elastic.fixX( boundnodes[i], uexact);
            elastic.fixY( boundnodes[i], vexact);
        }

        elastic.solve();

        vector<double> usol, vsol;
        elastic.getDisplacements(usol, vsol);

        size_t numnodes = mesh->getSize(0);
        double maxu = 0.0;
        double maxv = 0.0;
        for( int i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( !vtx->isBoundary() )  {
                const Point3D &xyz =  mesh->getNodeAt(i)->getXYZCoords();
                analyticalSol(xyz, uexact,vexact);
                double du = usol[i] - uexact;
                double dv = vsol[i] - vexact;
                if( fabs(uexact) > 1.0E-10)
                    maxu = max(maxu, fabs(du)/(fabs(uexact)));
                else
                    maxu = max(maxu, fabs(du));

                if( fabs(vexact) > 1.E-10)
                    maxv = max(maxv, fabs(dv)/(fabs(vexact)));
                else
                    maxv = max(maxv, fabs(dv));
            }
        }
        cout << "MaxError "  << maxu << " " << maxv << endl;
        uerrorLineEdit->setText( QString::number(maxu) );
        verrorLineEdit->setText( QString::number(maxv) );
        QDialog::update();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: solve()
{
    if( mesh == nullptr) return;

    setParams();

    int id  = testNameComboBox->currentIndex();

    if( id == 0 ) solveSingleNode();
    if( id == 1 ) solveSimpleTangle();
    if( id == 2 ) solveRandomTangle();
    if( id == 3 ) solveRotateQuads();
    if( id == 4 ) solveLinearElasticity();
}

///////////////////////////////////////////////////////////////////////////////
void JTangleFEMTestsDialog :: genSingleNodeMesh()
{
    mesh = JMesh::newObject();
    mesh->setName("TangleFEMModel");

    JNodeSequence nodes(5);
    nodes[0] = JNode::newObject();
    nodes[1] = JNode::newObject();
    nodes[2] = JNode::newObject();
    nodes[3] = JNode::newObject();
    nodes[4] = JNode::newObject();

    Point3D p3d;
    p3d[2] = 0.0;

    p3d[0] = 0.0;
    p3d[1] = 0.0;
    nodes[0]->setXYZCoords(p3d);

    p3d[0] = 1.0;
    p3d[1] = 0.0;
    nodes[1]->setXYZCoords(p3d);

    p3d[0] = 1.0;
    p3d[1] = 1.0;
    nodes[2]->setXYZCoords(p3d);

    p3d[0] = 0.0;
    p3d[1] = 1.0;
    nodes[3]->setXYZCoords(p3d);

    p3d[0] = 0.25;
    p3d[1] = 0.25;
    nodes[4]->setXYZCoords(p3d);

    mesh->addObjects(nodes);

    if( meshType == 3 ) {
        mesh->addObject( JTriangle::newObject(nodes[4], nodes[0], nodes[1] ));
        mesh->addObject( JTriangle::newObject(nodes[4], nodes[1], nodes[2] ));
        mesh->addObject( JTriangle::newObject(nodes[4], nodes[2], nodes[3] ));
        mesh->addObject( JTriangle::newObject(nodes[4], nodes[3], nodes[0] ));
    } else {
        mesh->addObject( JQuadrilateral::newObject(nodes[0], nodes[1], nodes[4], nodes[3] ));
        mesh->addObject( JQuadrilateral::newObject(nodes[1], nodes[2], nodes[3], nodes[4] ));
    }

    meshViewer->addObject(mesh);
    JNodeRenderPtr attrib;
    int err = nodes[4]->getAttribute("Render", attrib);
    if( !err) {
        attrib->glyph  = JNodeRender::NODE_AS_SPHERE;
        attrib->ballRadius = 0.020;
        attrib->color[0] = 1.0;
        attrib->color[1] = 0.0;
        attrib->color[2] = 0.0;
    }
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JTangleFEMTestsDialog :: genRotateQuadMesh()
{
    mesh = JMesh::newObject();
    mesh->setName("TangleFEMModel");

    JNodeSequence nodes(8);
    for( int i = 0; i < 8; i++) {
        nodes[i] = JNode::newObject();
        nodes[i]->setID(i);
    }
    Point3D xyz;
    xyz[2] = 0.0;

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    nodes[0]->setXYZCoords(xyz);

    xyz[0] = 1.0;
    xyz[1] = 0.0;
    nodes[1]->setXYZCoords(xyz);

    xyz[0] = 1.0;
    xyz[1] = 1.0;
    nodes[2]->setXYZCoords(xyz);

    xyz[0] = 0.0;
    xyz[1] = 1.0;
    nodes[3]->setXYZCoords(xyz);

    xyz[0] = 1.0/3.0;
    xyz[1] = 1.0/3.0;
    nodes[4]->setXYZCoords(xyz);

    xyz[0] = 2.0/3.0;
    xyz[1] = 1.0/3.0;
    nodes[5]->setXYZCoords(xyz);

    xyz[0] = 2.0/3.0;
    xyz[1] = 2.0/3.0;
    nodes[6]->setXYZCoords(xyz);

    xyz[0] = 1.0/3.0;
    xyz[1] = 2.0/3.0;
    nodes[7]->setXYZCoords(xyz);

    mesh->addObjects(nodes);

    JFaceSequence faces(5);
    faces[0] = JQuadrilateral::newObject( nodes[0], nodes[1], nodes[5], nodes[4] );
    faces[1] = JQuadrilateral::newObject( nodes[1], nodes[2], nodes[6], nodes[5] );
    faces[2] = JQuadrilateral::newObject( nodes[2], nodes[3], nodes[7], nodes[6] );
    faces[3] = JQuadrilateral::newObject( nodes[3], nodes[0], nodes[4], nodes[7] );
    faces[4] = JQuadrilateral::newObject( nodes[4], nodes[5], nodes[6], nodes[7] );
    mesh->addObjects( faces);

    if( meshType == 3)  {
        JMeshPtr trimesh = AllTriMeshGenerator::getFromQuadMesh(mesh,4);
        mesh = trimesh;
    }

    meshViewer->addObject(mesh);
    JNodeRenderPtr attrib;
    int err = nodes[4]->getAttribute("Render", attrib);
    if( !err) {
        attrib->glyph  = JNodeRender::NODE_AS_SPHERE;
        attrib->ballRadius = 0.020;
        attrib->color[0] = 1.0;
        attrib->color[1] = 0.0;
        attrib->color[2] = 0.0;
    }
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JTangleFEMTestsDialog :: genRandomMesh()
{
    JNodeSequence boundnodes;

    int    grid_dim[] = {20, 20};
    double length[]   = {1, 1};
    double origin[]   = {0.0, 0.0};

    mesh= AllQuadMeshGenerator::getStructuredMesh(grid_dim, length, origin);

    if( meshType == 3) {
        JMeshPtr trimesh = AllTriMeshGenerator::getFromQuadMesh(mesh,4);
        mesh = trimesh;
    }
    mesh->setName("TangleFEMModel");

    int numnodes = mesh->getSize(0);

    double xc = 0.5*length[0];
    double yc = 0.5*length[1];
    double rcut = length[0]/5.0;
    int nCount = 0;
    for( int i = 0; i < numnodes; i++) {
        JNodePtr v = mesh->getNodeAt(i);
        Point3D xyz = v->getXYZCoords();
        double dx = xyz[0] - xc;
        double dy = xyz[1] - yc;
        double r = sqrt(dx*dx + dy*dy);
        if( r < rcut)  {
            xyz[0] = xc + (-rcut + 2.0*rcut*drand48());
            xyz[1] = yc + (-rcut + 2.0*rcut*drand48());
            v->setXYZCoords(xyz);
            nCount++;
        }
    }
    meshViewer->addObject(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JTangleFEMTestsDialog :: genSimpleTangleMesh()
{
    mesh = JMesh::newObject();
    mesh->setName("TangleFEMModel");

    JNodeSequence nodes(8);
    for( int i = 0; i < 8; i++) {
        nodes[i] = JNode::newObject();
        nodes[i]->setID(i);
    }

    Point3D p3d;
    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;

    p3d[0] = 0.00;
    p3d[1] = 0.00;
    nodes[0]->setXYZCoords(p3d);
    p3d[0] = 1.00;
    p3d[1] = 0.00;
    nodes[1]->setXYZCoords(p3d);
    p3d[0] = 1.00;
    p3d[1] = 1.00;
    nodes[2]->setXYZCoords(p3d);
    p3d[0] = 0.00;
    p3d[1] = 1.00;
    nodes[3]->setXYZCoords(p3d);
    p3d[0] = 0.50;
    p3d[1] = 0.00;
    nodes[4]->setXYZCoords(p3d);
    p3d[0] = 0.50;
    p3d[1] = 1.00;
    nodes[5]->setXYZCoords(p3d);
    p3d[0] = 0.25;
    p3d[1] = 0.50;
    nodes[6]->setXYZCoords(p3d);
    p3d[0] = 0.75;
    p3d[1] = 0.50;
    nodes[7]->setXYZCoords(p3d);
    mesh->addObjects( nodes);
    mesh->addObject( JTriangle::newObject( nodes[0], nodes[4], nodes[6]) );
    mesh->addObject( JTriangle::newObject( nodes[1], nodes[7], nodes[4]) );
    mesh->addObject( JTriangle::newObject( nodes[2], nodes[5], nodes[7]) );
    mesh->addObject( JTriangle::newObject( nodes[3], nodes[6], nodes[5]) );
    mesh->addObject( JTriangle::newObject( nodes[1], nodes[2], nodes[7]) );
    mesh->addObject( JTriangle::newObject( nodes[0], nodes[6], nodes[3]) );
    mesh->addObject( JTriangle::newObject( nodes[4], nodes[7], nodes[6]) );
    mesh->addObject( JTriangle::newObject( nodes[5], nodes[6], nodes[7]) );

    if( meshType == 4) {
        AllQuadMeshGenerator qmesher;
        qmesher.setMesh(mesh);
        JMeshPtr qmesh = qmesher.getSimpleTris2Quads();
        mesh = qmesh;
    }

    meshViewer->addObject(mesh);
    JNodeRenderPtr attrib;
    int err = nodes[6]->getAttribute("Render", attrib);
    if( !err) {
        attrib->glyph  = JNodeRender::NODE_AS_SPHERE;
        attrib->ballRadius = 0.020;
        attrib->color[0] = 1.0;
        attrib->color[1] = 0.0;
        attrib->color[2] = 0.0;
    }

    err = nodes[7]->getAttribute("Render", attrib);
    if( !err) {
        attrib->glyph  = JNodeRender::NODE_AS_SPHERE;
        attrib->ballRadius = 0.020;
        attrib->color[0] = 0.0;
        attrib->color[1] = 1.0;
        attrib->color[2] = 1.0;
    }
    meshViewer->updateBuffers(mesh);

}

//////////////////////////////////////////////////////////////////////////
void JTangleFEMTestsDialog :: genBendingBeamMesh()
{
    int grid_dim[]   = {5, 4};    //  Mediaum
    double length[]  = {1, 0.2};
    double origin[]  = {0.0, 0.0};
    origin[0] = -0.5*length[0];
    origin[1] = -0.5*length[1];

    int gridresol;
    grid_dim[0] = 40;
    grid_dim[1] = 10;

    mesh = AllQuadMeshGenerator::getStructuredMesh(grid_dim, length, origin);

    if( meshType == 3)
        mesh  = AllTriMeshGenerator::getFromQuadMesh(mesh,4);

    mesh->setName("TangleFEMModel");

    /*
    elastic.setMesh(mesh);

    JEdgeSequence edges;
    int val;
    val  = 2;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
    elastic.setYForce(edges, 100.0);

    val  = 4;
    mesh->getEntities("Boundary", val, edges);
    assert( !edges.empty() ) ;
    elastic.fixEdges(edges);

        elastic.computeStress();
        vector<double> usol, vsol;
        elastic.getDisplacements(usol, vsol);
        cout << setprecision(10) << endl;
        int id;
        double sval;
        id = grid_dim[0]*grid_dim[1] - grid_dim[0];
        elastic.getValueAt( 'S', id, sval);
        cout << id << "  Stress Value at point D:  " <<  sval << endl;

        id = grid_dim[0]*ceil(0.5*grid_dim[1]) - 1;
        double uval, vval;
        elastic.getValueAt( 'U', id, uval);
        elastic.getValueAt( 'V', id, vval);
        cout << id << " Displacment values at point C:  " <<  uval << "  " << vval << endl;

        cout << "New Data stored in Udata.vtk" << endl;
        elastic.saveAs("Udata.vtk", "U");

        cout << "New Data stored in Vdata.vtk" << endl;
        elastic.saveAs("Vdata.vtk", "V");
        mesh->deleteAll();
        exit(0);
    */
    meshViewer->addObject(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: genSquareCircle1Mesh()
{
    JMeshNonlinearOptimization mopt;
    mopt.setBoundaryPreserve(1);

    double  xlength  = 100.0;
    double  ylength  = 40.0;
    double  radius   = 10.0;

    int  npoints = 100;
    double dtheta = 2*M_PI/( double)npoints;

    ofstream ofile( "test.poly", ios::out);

    ofile << npoints+4 << " 2  0  0" << endl;

    for( int i = 0; i < npoints; i++)  {
        double x = radius*cos(i*dtheta);
        double y = radius*sin(i*dtheta);
        ofile << i  << " " << x << "  " << y << endl;
    }

    ofile << npoints   << "  " << -0.5*xlength << " " << -0.5*ylength << endl;
    ofile << npoints+1 << "  " <<  0.5*xlength << " " << -0.5*ylength << endl;
    ofile << npoints+2 << "  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
    ofile << npoints+3 << "  " << -0.5*xlength << " " <<  0.5*ylength << endl;

    ofile << npoints + 4 <<  "  1 " << endl;
    for( int i = 0; i < npoints; i++)
        ofile << i << "  " << i << " " << (i+1)%npoints << " 5 " << endl;

    ofile << npoints   << "  " << npoints   << "  " << npoints+1 <<  "  1 " << endl;
    ofile << npoints+1 << "  " << npoints+1 << "  " << npoints+2 <<  "  2 " << endl;
    ofile << npoints+2 << "  " << npoints+2 << "  " << npoints+3 <<  "  3 " << endl;
    ofile << npoints+3 << "  " << npoints+3 << "  " << npoints   <<  "  4 " << endl;

    int landmark2 = npoints + 2;

    ofile << " 1 " << endl;
    ofile << " 0  0.0 0.0 " << endl;
    ofile.close();
    system( "triangle -peq30a5.0 test.poly");

    mesh = JMeshIO::readFile( "test.1.ele");
}
///////////////////////////////////////////////////////////////////////////////
void JTangleFEMTestsDialog :: genSquareCircle2Mesh()
{
    int    grid_dim[] = {4, 4};
    double length[]   = {10, 4};
    double origin[]   = {0.0, 0.0};

    mesh= AllQuadMeshGenerator::getStructuredMesh(grid_dim, length, origin);
    JEdgePtr edge;
    JFacePtr face;

    JEdgeSequence edges;
    JNodeSequence nodes;
    int   val;

//  Face-0
    face = mesh->getFaceAt(0);

    val = 1;
    edge= face->getEdgeAt(0);
    face->setAttribute("Boundary", val);

    val = 4;
    edge = face->getEdgeAt(3);
    edge->setAttribute("Boundary", val);

//  Face-1
    face = mesh->getFaceAt(1);
    val  = 1;
    edge = face->getEdgeAt(0);
    edge->setAttribute("Boundary", val);

    val  = 5;
    edge = face->getEdgeAt(2);
    edge->setAttribute("Boundary", val);

//  Face-2
    face = mesh->getFaceAt(2);

    val = 1;
    edge= face->getEdgeAt(0);
    face->setAttribute("Boundary", val);

    val = 2;
    edge = face->getEdgeAt(1);
    edge->setAttribute("Boundary", val);

//  Face-3
    face = mesh->getFaceAt(3);

    val = 5;
    edge= face->getEdgeAt(1);
    face->setAttribute("Boundary", val);

    val = 4;
    edge = face->getEdgeAt(3);
    edge->setAttribute("Boundary", val);

//   Face-5
    face = mesh->getFaceAt(3);

    val = 2;
    edge= face->getEdgeAt(1);
    face->setAttribute("Boundary", val);

    val = 5;
    edge = face->getEdgeAt(3);
    edge->setAttribute("Boundary", val);

//   Face-6
    face = mesh->getFaceAt(6);

    val = 3;
    edge= face->getEdgeAt(2);
    face->setAttribute("Boundary", val);

    val = 4;
    edge = face->getEdgeAt(3);
    edge->setAttribute("Boundary", val);

//   Face-7
    face = mesh->getFaceAt(7);

    val = 5;
    edge= face->getEdgeAt(0);
    face->setAttribute("Boundary", val);

    val = 3;
    edge = face->getEdgeAt(2);
    edge->setAttribute("Boundary", val);

//   Face-7
    face = mesh->getFaceAt(8);

    val = 2;
    edge= face->getEdgeAt(1);
    face->setAttribute("Boundary", val);

    val = 3;
    edge = face->getEdgeAt(2);
    edge->setAttribute("Boundary", val);

    mesh->getFaceAt(4)->setStatus(JMeshEntity::REMOVE);

    JQuadRefiner refiner;
    refiner.setMesh(mesh);
    JNodeSequence newnodes;
    JFaceSequence newfaces;

    int dim[2];
    dim[0] = 10;
    dim[1] = 10;
    refiner.refineAll(dim);
#ifdef CSV

    val  = 5;
    mesh->getEntities("Boundary", val, edges);
    JMeshTopology::getEntitySet( edges, nodes);

    double xc = 0.0;
    double yc = 0.0;
    int nsize = nodes.size();
    for( int i = 0; i < nsize; i++) {
        JNodePtr v = nodes[i];
        Point3D xyz = v->getXYZCoords();
        xc += xyz[0];
        yc += xyz[1];
    }

    xc = xc/(double)nsize;
    yc = yc/(double)nsize;
    double radius = 10.0;

    circleCenter[0] = xc;
    circleCenter[1] = yc;
    circleRadius    = radius;

//   double start_angle = 2.0*M_PI/(double)nodes.size();
    double start_angle = rotAngle*M_PI/180.0;

    int    pos = 0;
    double mindist = std::numeric_limits<double>::max();
    double t, dx, dy, dl;
    Point3D xyz;

    for( int i = 0; i < nodes.size(); i++) {
        JNodePtr v = nodes[i];
        xyz = v->getXYZCoords();
        dx = xyz[0] - xc;
        dy = xyz[1] - yc;

        t = atan2(dy,dx);
        xyz[0] = xc + radius*cos(t + start_angle);
        xyz[1] = yc + radius*sin(t + start_angle);
        xyz[2] = 0.0;
        v->setAttribute("TargetPos", xyz);

        dx = xyz[0] - xc ;
        dy = xyz[1] - yc - radius;
        dl = sqrt(dx*dx + dy*dy);
        if( dl < mindist) {
            pos = i;
            mindist = dl;
        }
    }
    int landmark = nodes[pos]->getID();
    xyz =  nodes[pos]->getXYZCoords();

    localmap.setMaxIterations(10000);
    localmap.solve();

    JMeshNonlinearOptimization mopt;
    mopt.setBoundaryPreservation(1);
    mopt.setNumIterations(100);
    mopt.setMesh(mesh);
    mopt.improveShapes();

    for( size_t i = 0; i < nodes.size(); i++) {
        JNodePtr v = nodes[i];
        Point3D xyz = v->getXYZCoords();
        double dx = xyz[0] - xc;
        double dy = xyz[1] - yc;
        double r  = sqrt(dx*dx + dy*dy);
        assert( fabs(r - radius) < 1.0E-06);
    }

    circleRadius = radius;
#endif

#ifdef CSV
    if( randomTangle) {
        double rcut = 10.0;
        int nCount = 0;
        int numnodes = mesh->getSize(0);
        for( int i = 0; i < numnodes; i++) {
            JNodePtr v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            double xt = 15.0;
            double yt = 20.0;
            double dx = xyz[0] - xt;
            double dy = xyz[1] - yt;
            double r = sqrt(dx*dx + dy*dy);
            if( r < rcut)  {
                xyz[0] = xt + (-rcut + 2.0*rcut*drand48());
                xyz[1] = yt + (-rcut + 2.0*rcut*drand48());
                v->setXYZCoords(xyz);
                nCount++;
            }
//         if( nCount == 10) break;
        }
    }


    if( randomTangle ) {
        val  = 5;
        mesh->getEntities("Boundary", val, edges);

        QuadChord qChord;
        qChord.setMesh(mesh);
        qChord.setSeed( edges[0] );
        JFaceSequence qFaces = qChord.getMeshFaces();

        qChord.setSeed( qFaces[3]->getEdgeAt(1) );
        qFaces = qChord.getMeshFaces();

        EdgeSet eset;
        for( int i = 0; i < qFaces.size(); i++) {
            eset.insert( qFaces[i]->getEdgeAt(1));
            eset.insert( qFaces[i]->getEdgeAt(3));
        }
        map<JNodePtr, Point3D> orgPoints;

        foreach_( JEdgePtr edge, eset) {
            JNodePtr v0 = edge->getNodeAt(0);
            JNodePtr v1 = edge->getNodeAt(1);
            orgPoints[v0] = v0->getXYZCoords();
            orgPoints[v1] = v1->getXYZCoords();
        }

        cout << "#Secondary Gauss Points (1,3,7)" << endl;
        cin >> numGauss2;

        double dr = 0.043;
        Point3D p0, p1;

        for( int i = 0; i < 25; i++) {
            double r = i*dr;
            outfile << r << " ";
            foreach_(JEdgePtr edge, eset) {
                JNodePtr v0 = edge->getNodeAt(0);
                JNodePtr v1 = edge->getNodeAt(1);

                p0[0]   = (1.0-r)*orgPoints[v0][0] + r*orgPoints[v1][0];
                p0[1]   = (1.0-r)*orgPoints[v0][1] + r*orgPoints[v1][1];
                p0[2]   = 0.0;

                p1[0]   = (1.0-r)*orgPoints[v1][0] + r*orgPoints[v0][0];
                p1[1]   = (1.0-r)*orgPoints[v1][1] + r*orgPoints[v0][1];
                p1[2]   = 0.0;

                v0->setXYZCoords(p0);
                v1->setXYZCoords(p1);
            }
//          patchTest(mesh);
            LinearElasticitySquareCircle2( mesh);
            ostringstream oss;
            oss << "anim";
            if( i < 10 ) oss << "0";
            oss << i << ".vtk";
            mesh->saveAs(oss.str());
        }
        exit(0);
    }
    LinearElasticitySquareCircle2( mesh);
//  patchTest(mesh);
#endif
}
///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: solveSingleNode()
{
    if( problemID != 0 ) return;
    Point3D xyz;
    xyz[1]  = 0.5;
    xyz[2]  = 0.0;

    double dx = 0.1;

    JNodePtr v4 = mesh->getNodeAt(4);

    for( int j = 1; j < 10; j++)
        for( int i = 1; i < 10; i++) {
            xyz[0] = i*dx;
            xyz[1] = j*dx;
            v4->setXYZCoords(xyz);
            meshViewer->updateBuffers(mesh);
            glFlush();
            patchTest();
            sleep(1.0);
        }
}
///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: solveSimpleTangle()
{
    if( problemID != 1 ) return;
    Point3D xyz;
    xyz[1]  = 0.5;
    xyz[2]  = 0.0;

    JNodePtr v6 = mesh->getNodeAt(6);
    JNodePtr v7 = mesh->getNodeAt(7);

    double dx = 0.01;

    for( int i = 0; i < 95; i++) {
        xyz[0] = 0.01 + i*dx;
        v6->setXYZCoords(xyz);

        xyz[0] = 0.99 - i*dx;
        v7->setXYZCoords(xyz);
        meshViewer->updateBuffers(mesh);
        glFlush();
        patchTest();
        sleep(0.5);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: solveRandomTangle()
{
    if( problemID != 2 ) return;
    patchTest();
}

///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: solveRotateQuads()
{
    if( problemID != 3 )  return;

    double angle, radius = 1.0/3.0;
    double dtheta = 5.0*M_PI/180.0;

    Point3D xyz;
    for( int i = 0; i < 72; i++) {
        patchTest();
        for( int j = 0; j < 4; j++) {
            const JNodePtr &vtx =  mesh->getNodeAt(4+j);
            xyz = vtx->getXYZCoords();
            angle = atan2(xyz[1] - 0.5, xyz[0]-0.5) + dtheta;
            xyz[0] = 0.5 + radius*cos(angle);
            xyz[1] = 0.5 + radius*sin(angle);
            vtx->setXYZCoords( xyz );
        }
        meshViewer->updateBuffers(mesh);
        sleep(0.5);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: solveLinearElasticity()
{
    /*
        if( problemID != 4) return;
        elastic.solve();
    */
}

///////////////////////////////////////////////////////////////////////////////
void JTangleFEMTestsDialog :: setShapeOrder()
{
    /*
        int order;
        if ( linearShapeFuncRadioButton->isChecked() )    order = 1;
        if ( quadraticShapeFuncRadioButton->isChecked() ) order = 2;
        elastic.setOrder(order);
        meshViewer->updateBuffers(mesh);
        displayQuadraticNodes();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JTangleFEMTestsDialog :: displayQuadraticNodes()
{
    /*
        JNodeSequence nodes;
        elastic.getQuadraticNodes(nodes);

        bool val = displayQuadraticNodesCheckBox->isChecked();
        JNodeRenderPtr attrib;

        size_t numnodes = nodes.size();
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = nodes[i];
            int err = vtx->getAttribute("Render", attrib);
            if( !err ) attrib->display = val;
        }
        meshViewer->updateBuffers(mesh);
    */
}
///////////////////////////////////////////////////////////////////////////////
void JTangleFEMTestsDialog :: makeConnections()
{
    CheckBox( displayQuadraticNodesCheckBox, [=] {displayQuadraticNodes();});

    ComboBox( testNameComboBox, [=] {setMesh(); });

    RadioButton( trimeshRadioButton, [=] {setMesh();});
    RadioButton( quadmeshRadioButton,[=] {setMesh();});
    RadioButton( quadraticShapeFuncRadioButton, [=] {setShapeOrder();});
    RadioButton( linearShapeFuncRadioButton,  [=] {setShapeOrder(); });

    PushButton( applyPushButton, [=] {solve();});

    PushButton( closePushButton, [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
