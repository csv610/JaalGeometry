#include "SurfaceParameterization.hpp"

JSurfaceParameterization :: JSurfaceParameterization()
{
    boundarytype = 0;
    weighttype   = 0;
    numIters     = 2000;
    gammaP       = 1.0;
    smooth       = 1;
    intrinsiclambda=0.5;
}
///////////////////////////////////////////////////////////////////////////////
int JSurfaceParameterization :: setMesh( const JMeshPtr &m)
{
    mesh = m;
}

///////////////////////////////////////////////////////////////////////////////

void JSurfaceParameterization :: genCharts()
{
    if( mesh == nullptr) return;
    invertedUVMeshes.clear();

    if( meshPart == nullptr) meshPart.reset( new JMetisPartitioner);

    meshPart->setMesh(mesh);
    meshPart->getTopologicalDisks();
    int numParts = meshPart->getNumPartitions();

    if( numParts ) {
        submeshes.resize(numParts);
#pragma omp parallel for
        for( int i = 0; i < numParts; i++) {
            submeshes[i] = meshPart->getSubMesh(i);
            setUVMesh( submeshes[i] );
        }
    } else {
        vector<JEdgeSequence> edges;
        mesh->getTopology()->getBoundary(edges);
        if( edges.size() == 1) {
            submeshes.resize(1);
            submeshes[0] = mesh;
            setUVMesh( mesh  );
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void JSurfaceParameterization :: normalize( const JMeshPtr &uvmesh)
{
    if( uvmesh == nullptr) return;

    JMeshAffineTransform affine;
    affine.setMesh( uvmesh );
    affine.toCenter();
    affine.normalize();
}

///////////////////////////////////////////////////////////////////////////////
double JSurfaceParameterization :: getIrregularity( const JMeshPtr &uvmesh)
{
   // Ref:Hierarchical Face Clustering on Polygonal Surfaces
   double A  = uvmesh->getGeometry()->getSurfaceArea();
   JEdgeSequence edges;
   uvmesh->getTopology()->getBoundary(edges);
   double l = JEdgeGeometry::getLength(edges);
   return l*l/(4.0*M_PI*A);
}
///////////////////////////////////////////////////////////////////////////////

int JSurfaceParameterization :: writePly2File(const JMeshPtr &submesh, const string &filename)
{
    if( submesh  == nullptr) return 1;

    size_t numnodes = submesh->getSize(0);
    assert( numnodes);
    size_t numfaces = submesh->getSize(2);
    assert( numfaces);

    int numTriangles = 0;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = submesh->getFaceAt(i);
        if( face->isActive() ) {
            int nn = face->getSize(0);
            if( nn == 3) numTriangles += 1;
            if( nn == 4) numTriangles += 2;
        }
    }

    ofstream ofile( filename.c_str(), ios::out);
    if( ofile.fail()) return 2;

    ofile << numnodes  << endl;
    ofile << numTriangles << endl;

    map<JNodePtr,size_t> localID;

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = submesh->getNodeAt(i);
        const Point3D &xyz = vtx->getXYZCoords();
        ofile << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
        localID[vtx] = i;
    }

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = submesh->getFaceAt(i);
        if( face->isActive() ) {
            int nn = face->getSize(0);
            if( nn == 3) {
                size_t  id0 = localID[face->getNodeAt(0)];
                size_t  id1 = localID[face->getNodeAt(1)];
                size_t  id2 = localID[face->getNodeAt(2)];
                ofile << "3 " << id0 << " " << id1 << " " << id2 << endl;
            }
            if( nn == 4) {
                size_t  id0 = localID[face->getNodeAt(0)];
                size_t  id1 = localID[face->getNodeAt(1)];
                size_t  id2 = localID[face->getNodeAt(2)];
                size_t  id3 = localID[face->getNodeAt(3)];
                ofile << "3 " << id0 << " " << id1 << " " << id2 << endl;
                ofile << "3 " << id0 << " " << id2 << " " << id3 << endl;
            }
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JSurfaceParameterization :: readPly2File( const  JMeshPtr &submesh, const string &fname)
{
    ifstream ifile( fname.c_str(), ios::in);
    if( ifile.fail()) return nullptr;

    JMeshPtr paramMesh = JMesh::newObject();

    size_t numnodes, numfaces;
    ifile >> numnodes;
    ifile >> numfaces;

    Point3D xyz;
    map<JNodePtr,JNodePtr>  nodeMap;
    for( size_t i = 0; i < numnodes; i++) {
        ifile >> xyz[0] >> xyz[1] >> xyz[2];
        JNodePtr newnode = JNode::newObject();
        newnode->setXYZCoords(xyz);
        paramMesh->addObject(newnode);
        const JNodePtr &vtx = submesh->getNodeAt(i);
        nodeMap[vtx] = newnode;
    }

    // We do not read the faces from the output. It comes from the original mesh ....
    numfaces = submesh->getSize(2);
    JNodeSequence nodes;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &oldface = submesh->getFaceAt(i);
        int nn = oldface->getSize(0);
        nodes.resize(nn);
        for( int j = 0; j < nn; j++) {
            nodes[j] = nodeMap[oldface->getNodeAt(j)];
        }
        JFacePtr newface = oldface->getClone();
        newface->setNodes( nodes );
        paramMesh->addObject(newface);
    }
    normalize(paramMesh);
    return paramMesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JSurfaceParameterization :: getTriMesh(const JMeshPtr &submesh)
{
    JMeshPtr triMesh;
    int elemType = submesh->getTopology()->getElementsType(2);
    if( elemType == JFace::TRIANGLE) {
        triMesh = submesh->deepCopy();
        assert( triMesh->getSize(0) == submesh->getSize(0) );
        assert( triMesh->getSize(2) == submesh->getSize(2) );
        triMesh->enumerate(0);
        triMesh->enumerate(2);
    }

    if( elemType == JFace::QUADRILATERAL) {
        JMeshPtr qmesh = submesh->deepCopy();
        AllTriMeshGenerator alltri;
        triMesh = alltri.getFromQuadMesh(qmesh);
        assert( triMesh->getSize(0) == submesh->getSize(0) );
        assert( triMesh->getSize(2) == 2.0*submesh->getSize(2) );
        triMesh->enumerate(0);
        triMesh->enumerate(2);
    }
    return triMesh;
}
///////////////////////////////////////////////////////////////////////////////

void JSurfaceParameterization :: LeastSquareConformalMap( const JMeshPtr &smesh)
{
#ifdef USE_IGL
    JMeshPtr triMesh = getTriMesh(smesh); // Since things work only for the triangle mesh..
    if( triMesh == nullptr) return;

    JMeshEigenMatrix mat;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    mat.setMesh(triMesh);
    V = mat.getNodeMatrix();
    F = mat.getFaceMatrix();

    // Fix two points on the boundary
    Eigen::VectorXi bnd, b(2,1);
    igl::boundary_loop(F, bnd);

    b(0) = bnd(0);
    b(1) = bnd(round(bnd.size()/2));

    Eigen::MatrixXd bc(2,2);
    bc<< 0,0,1,0;

    Eigen::MatrixXd uv;
    // LSCM parametrization
    igl::lscm(V, F, b, bc, uv);

    JMeshPtr uvMesh = smesh->deepCopy();
    JNodeSequence uvnodes = uvMesh->getNodes();
    Point3D xyz;
    int index = 0;
    for( const JNodePtr &vtx : uvnodes) {
        xyz[0] = uv.coeff(index,0);
        xyz[1] = uv.coeff(index,1);
        xyz[2] = 0.0;
        vtx->setXYZCoords(xyz);
        index++;
    }

//  rescaleUVMesh(smesh, uvMesh, NATURAL_BOUNDARY);
    uvMesh->setName("UVLSCM");
    normalize(uvMesh);
    smesh->setAttribute("UVMesh", uvMesh);
#endif
}

///////////////////////////////////////////////////////////////////////////////
void JSurfaceParameterization :: AsRigidAsPossible( const JMeshPtr &smesh)
{

#ifdef USE_IGL
    JMeshPtr triMesh = getTriMesh(smesh); // Since things work only for the triangle mesh..
    if( triMesh == nullptr) return;

    JMeshEigenMatrix mat;

    mat.setMesh(triMesh);
    Eigen::MatrixXd V = mat.getNodeMatrix();
    Eigen::MatrixXi F = mat.getFaceMatrix();

    Eigen::VectorXi bnd;
    igl::boundary_loop(F,bnd);
    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V,bnd,bnd_uv);

    Eigen::MatrixXd initial_guess;
    igl::harmonic(V,F,bnd,bnd_uv,1,initial_guess);

    // Add dynamic regularization to avoid to specify boundary conditions
    igl::ARAPData arap_data;
    arap_data.with_dynamics = true;
    Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
    Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);

    // Initialize ARAP
    arap_data.max_iter = 100;
    // 2 means that we're going to *solve* in 2d
    arap_precomputation(V,F,2,b,arap_data);

    // Solve arap using the harmonic map as initial guess
    Eigen::MatrixXd uv = initial_guess;

    arap_solve(bc,arap_data,uv);

    JMeshPtr uvMesh = smesh->deepCopy();
    JNodeSequence uvnodes = uvMesh->getNodes();
    Point3D xyz;
    int index = 0;
    for( const JNodePtr &vtx : uvnodes) {
        xyz[0] = uv.coeff(index,0);
        xyz[1] = uv.coeff(index,1);
        xyz[2] = 0.0;
        vtx->setXYZCoords(xyz);
        index++;
    }

//  rescaleUVMesh(smesh, uvMesh, NATURAL_BOUNDARY);
    uvMesh->setName("UVARAP");
    normalize(uvMesh);
    smesh->setAttribute("UVMesh", uvMesh);
#endif
}

///////////////////////////////////////////////////////////////////////////////
void JSurfaceParameterization :: rescaleUVMesh(const JMeshPtr &xyzMesh, const JMeshPtr &uvMesh, int boundarytype)
{
    if( rescaleByArea == 0) return;

    double xyzArea = fabs(xyzMesh->getGeometry()->getSurfaceArea());
    double uvArea  = fabs(uvMesh->getGeometry()->getSurfaceArea());

    JMeshAffineTransform affine;
    affine.setMesh(uvMesh);
    affine.toCenter();

    switch( boundarytype)
    {
    case SQUARE:
    {
        double sc = sqrt(xyzArea)/sqrt(uvArea);
        affine.scale(sc, sc, 1);
    }
    break;
    case CIRCLE:
    {
        double rx = sqrt( xyzArea/M_PI);
        double ru = sqrt( uvArea/M_PI);
        double sc =  rx/ru;
        affine.scale(sc, sc, 1);
    }
    break;
    default:
    {
        JFacePtr  quad = JMeshGeometry::getMinRectangle(uvMesh->getNodes());
        Point3D p1 = quad->getNodeAt(1)->getXYZCoords();
        Point3D p0 = quad->getNodeAt(0)->getXYZCoords();
        double  dx = p1[0] - p0[0];
        double  dy = p1[1] - p0[1];
        double  dt = atan2(dy,dx);
        affine.rotate( -dt, 2);
    }
    }
//    double uvAfter = fabs(uvMesh->getGeometry()->getSurfaceArea());
//    cout << "XYZ Area " << xyzArea << " UV before " << uvArea << " UV After " << uvAfter << endl;

}
///////////////////////////////////////////////////////////////////////////////

void JSurfaceParameterization :: setUVMesh( const JMeshPtr &submesh)
{
    if( submesh == nullptr) return;

    if( submesh->getTopology()->getEulerCharacteristic() != 1) {
        cout << "Warning: Submesh is not a topological disk " << endl;
        return;
    }

    JMeshPtr uvMesh;
    switch( weighttype )
    {
    case LEAST_SQUARE_CONFORMAL:
        LeastSquareConformalMap(submesh);
        break;
    case AS_RIGID_AS_POSSIBLE:
        AsRigidAsPossible(submesh);
        break;
    default:
        string infile  = "in.ply2";
        string outfile = "out.ply2";

        numIters = submesh->getSize(0);

        writePly2File( submesh, infile );

        ostringstream oss;
        oss << "/home/csverma/Disk/Software/Mesh/SurfParameterization/bin/surfparam ";
        oss << "-b " << boundarytype << " ";
        oss << "-g " << gammaP << " ";
        oss << "-l " << intrinsiclambda << " ";
        oss << "-n " << numIters << " ";
        oss << "-p " << paramtype;
        oss << "-s " << smooth << " ";
        oss << "-w " << weighttype << " ";
        oss << "-i " << infile  << " ";
        oss << "-o " << outfile;

        string cmd = oss.str();

        int stat = system(cmd.c_str() );
        if( stat != -1 ) {
            uvMesh = readPly2File( submesh, outfile );
            rescaleUVMesh( submesh, uvMesh, boundarytype);
            submesh->setAttribute("UVMesh", uvMesh);
        }
        break;
    }

    int err = submesh->getAttribute("UVMesh", uvMesh);
    if( err ) return;

    int nCount = 0;
    for( size_t i = 0; i < uvMesh->getSize(2); i++) {
        const JFacePtr &face = uvMesh->getFaceAt(i);
        double area = JFaceGeometry::getSignedArea(face);
        if( area < 0.0) nCount++;
    }

    nCount = 0;
    if( nCount > uvMesh->getSize(2)/2 ) {
        uvMesh->getTopology()->reverseAll();
        for( size_t i = 0; i < uvMesh->getSize(2); i++) {
            const JFacePtr &face = uvMesh->getFaceAt(i);
            double area = JFaceGeometry::getSignedArea(face);
            if( area < 0.0) nCount++;
        }
    }

    if( nCount ) {
        cout << "Warning: Out of total " << uvMesh->getSize(2) << " there are " << nCount << " Inverted faces in the UV mesh " << endl;
        invertedUVMeshes.push_back(uvMesh);
    }

    cout << "UV Patch Irregularity " << getIrregularity(uvMesh) << endl;

}

///////////////////////////////////////////////////////////////////////////////
