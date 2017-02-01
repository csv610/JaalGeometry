#include "MeshLaplacian.hpp"

////////////////////////////////////////////////////////////////////////////////

JLaplaceMeshSmoother ::  JLaplaceMeshSmoother()
{
    numIterations = 1;
    Kpb    = 0.1;
    lambda = 0.6307;
    mu     = -0.6732;
    setEdgesWeight( COMBINATORIAL );
    tolerance = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

double TutteWeight::getWeight(const JEdgePtr &edge, int dir)
{
    int ndeg = 0;

    JNodePtr vtx = nullptr;

    switch(dir) {
    case 1:
        vtx = edge->getNodeAt(0);
        ndeg = vtx->getNumRelations(0);
        assert( ndeg );
        break;
    case -1:
        vtx = edge->getNodeAt(1);
        ndeg = vtx->getNumRelations(0);
        assert( ndeg );
        break;
    }
    return 1.0/(double)ndeg;
}
////////////////////////////////////////////////////////////////////////////////
double NormalizedWeight::getWeight(const JEdgePtr &edge, int)
{
    JNodePtr vtx = nullptr;
    vtx = edge->getNodeAt(0);
    int di = vtx->getNumRelations(0);

    vtx = edge->getNodeAt(1);
    int dj = vtx->getNumRelations(0);
    return 1.0/sqrt(di*dj);
}
////////////////////////////////////////////////////////////////////////////////

double CotangentWeight:: getCotangent(const Point3D &org, const Point3D &p1, const Point3D &p2) const
{
    Vec3D u, v;
    JMath::make_vector(p1, org, u);
    JMath::make_vector(p2, org, v);

    double udotv = JMath::dot_product(u, v);
    double udotu = JMath::dot_product(u, u);
    double vdotv = JMath::dot_product(v, v);
    double denom = udotu * vdotv - udotv*udotv;

    if (denom < 1.0E-15) return 0;
    return udotv / sqrt(denom);
}
///////////////////////////////////////////////////////////////////////////////
double CotangentWeight::getWeight(const JEdgePtr &edge, int )
{
    static JFaceSequence faceNeighs;

    JEdge::getRelations(edge, faceNeighs);
    if( faceNeighs.empty() || faceNeighs.size() > 2) return 1;

    JNodePtr v0 = edge->getNodeAt(0);
    JNodePtr v1 = edge->getNodeAt(1);

    JNodePtr ov = nullptr;
    ov = JTriangle::getOppositeNode(faceNeighs[0], v0, v1);
    double cota = getCotangent(ov->getXYZCoords(), v0->getXYZCoords(), v1->getXYZCoords());

    if( faceNeighs.size() == 1) return 0.5*cota;

    ov = JTriangle::getOppositeNode(faceNeighs[1], v0, v1);
    double cotb = getCotangent(ov->getXYZCoords(), v0->getXYZCoords(), v1->getXYZCoords());
    return 0.5*(cota + cotb);
}
///////////////////////////////////////////////////////////////////////////////

double FloaterMeanValueWeight::getWeight(const JEdgePtr &edge, int dir)
{
    static JFaceSequence faceNeighs;

    JEdge::getRelations(edge, faceNeighs);
    if( faceNeighs.empty() || faceNeighs.size() > 2) return 1;

    JNodePtr v0 = edge->getNodeAt(0);
    JNodePtr v1 = edge->getNodeAt(1);

    double alpha1 = 0.0;
    double alpha2 = 0.0;

    switch(dir) {
    case 1:
        alpha1 =  JFaceGeometry::getAngleAt( faceNeighs[0], v0, ANGLE_IN_RADIANS);
        alpha2 =  JFaceGeometry::getAngleAt( faceNeighs[1], v0, ANGLE_IN_RADIANS);
        break;
    case -1:
        alpha1 =  JFaceGeometry::getAngleAt( faceNeighs[0], v1, ANGLE_IN_RADIANS);
        alpha2 =  JFaceGeometry::getAngleAt( faceNeighs[1], v1, ANGLE_IN_RADIANS);
        break;
    }
    double len = JEdgeGeometry::getLength(edge);
    double w   =  (tan(0.5*alpha1) + tan(0.5*alpha2))/len;
    return w;
}

//////////////////////////////////////////////////////////////////////////////

int JLaplaceMeshSmoother ::setMesh(const JMeshPtr &m)
{
    mesh = m;

    if( mesh == nullptr ) return 1;

    mesh->pruneAll();

    mesh->buildRelations(0,0);

    weightMatrix.reset(new Eigen::SparseMatrix<double>);
    size_t numNodes = mesh->getSize(0);
    weightMatrix->resize(numNodes,numNodes);

    size_t numEdges = mesh->getSize(1);
    mesh->enumerate(0);
    weightMatrix->reserve(2*numEdges);

    topDim = mesh->getTopology()->getDimension();

    // If the topology dimension is 3, then all the boundary faces and vertices
    // are constrained. If the dimension is 2, then all the boundary edges and
    // their nodes are constrained. They can however be smoothed by other methods.
    if( topDim == 3 ) { }

    return 0;
}

//////////////////////////////////////////////////////////////////////////////
void JLaplaceMeshSmoother :: setEdgesWeight( int w, bool update_every_step)
{
    edgeWeightID = w;
    update_weight_every_step = update_every_step;

    switch( w ) {
    case COMBINATORIAL:
        edgeWeight.reset(new CombinatorialWeight());
        break;
    case COTANGENT:
        edgeWeight.reset(new CotangentWeight());
        break;
    case EDGE_LENGTH:
        edgeWeight.reset(new EdgeLengthWeight());
        break;
    case FLOATER_MEAN_VALUE:
        edgeWeight.reset(new FloaterMeanValueWeight());
        break;
        break;
    case NORMALIZED:
        edgeWeight.reset(new NormalizedWeight());
        break;
    case TUTTE:
        edgeWeight.reset(new TutteWeight());
        break;
    default:
        edgeWeight.reset( new CombinatorialWeight());
        break;
    }
}
//////////////////////////////////////////////////////////////////////////////

int JLaplaceMeshSmoother :: getWeightMatrix()
{
    if( mesh == nullptr) return 1;

    string name = edgeWeight->getName();

    size_t numNodes = mesh->getSize(0);
    size_t numEdges = mesh->getSize(1);

    vector<Triplet_t> matCoeff;
    matCoeff.reserve(2*numEdges);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            size_t i  = edge->getNodeAt(0)->getID();
            size_t j  = edge->getNodeAt(1)->getID();
            matCoeff.push_back(Triplet_t(i, j, 1.0));
            matCoeff.push_back(Triplet_t(j, i, 1.0));
        }
    }
    weightMatrix->setFromTriplets(matCoeff.begin(), matCoeff.end() );

    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            double val = edgeWeight->getWeight(edge, 1);
            size_t i   = edge->getNodeAt(0)->getID();
            size_t j   = edge->getNodeAt(1)->getID();
            weightMatrix->coeffRef(i,j) = val;
            weightMatrix->coeffRef(j,i) = val;
        }
    }

    JNodeSequence neighs;
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            int iv = vertex->getID();
            JNode::getRelations(vertex, neighs);
            double sum = 0;
            for( size_t j = 0; j < neighs.size(); j++) {
                int jv = neighs[j]->getID();
                sum += weightMatrix->coeff(iv,jv);
            }
            weightMatrix->insert(i,i) = -sum;
        }
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

void JLaplaceMeshSmoother :: update(const JNodePtr &vtx)
{
    size_t id = vtx->getID();

    // No check; Just update, you are your own risk...
    if( !check_element_inversion ) {
        vtx->setXYZCoords(newPos[id] );
        return;
    }

    // If the new position results in positive face area, just updaate it.
    if( newNodeQuality[id] > 0.0) {
        vtx->setXYZCoords(newPos[id] );
        return;
    }

    // If the new position results in improving the bad quality (less negative
    // area, update it..
    if( newNodeQuality[id] > oldNodeQuality[id] ) {
        vtx->setXYZCoords(newPos[id] );
        return;
    }

    // You came here means that the face area is becoming more negative, therefore
    // Keep the position as it is ...

}
//////////////////////////////////////////////////////////////////////////////

int
JLaplaceMeshSmoother ::atomicOp( const JNodePtr &vertex, const JNodeSequence &neighs, double coeff)
{
    int numNeighs = neighs.size();
    assert(numNeighs);

    const Point3D &papex = vertex->getXYZCoords();

    if( numNeighs < 1) return 1;
    int iv = vertex->getID();

    Point3D xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    double wsum = 0.0;
    for(int j = 0; j < numNeighs; j++) {
        const Point3D &p3d = neighs[j]->getXYZCoords();
        int jv = neighs[j]->getID();
        double wij = weightMatrix->coeff(iv,jv);
        xyz[0] += wij*p3d[0];
        xyz[1] += wij*p3d[1];
        xyz[2] += wij*p3d[2];
        wsum   += wij;
    }

    if( fabs(wsum) < 1.0E-10) return 2;

    xyz[0] /= wsum;
    xyz[1] /= wsum;
    xyz[2] /= wsum;

    double dx = xyz[0] - papex[0];
    double dy = xyz[1] - papex[1];
    double dz = xyz[2] - papex[2];

    maxDiff = max(maxDiff, dx*dx + dy*dy + dz*dz );

    xyz[0] = papex[0] + coeff*dx;
    xyz[1] = papex[1] + coeff*dy;
    xyz[2] = papex[2] + coeff*dz;

    newPos[iv] = xyz;
    return 0;
}

//////////////////////////////////////////////////////////////////////////////
int JLaplaceMeshSmoother :: smooth( const JNodeSequence &freenodes)
{
    getWeightMatrix();
    if( weightMatrix == nullptr) {
        cout << "Warning: Edge-Weight matrix is null: No smoothing done " << endl;
        return 1;
    }

    size_t numnodes = mesh->getSize(0);
    newPos.resize(numnodes);
    for( size_t i = 0; i < numnodes; i++)
        newPos[i] = mesh->getNodeAt(i)->getXYZCoords();
    // If inverion is checked, get the minimum face area at every vertex ...
    if( check_element_inversion ) {
        mesh->buildRelations(0,topDim);
        getQuality();
    }

    JNodeSequence nodeneighs;
    for( int j = 0; j < numIterations; j++) {
        if( update_weight_every_step ) getWeightMatrix();

        size_t ncount = 0;
        maxDiff = 0.0;
        for( const JNodePtr &vertex: freenodes) {
            JNode::getRelations(vertex, nodeneighs );
            int err = atomicOp( vertex, nodeneighs, lambda);
            if( err == 0) ncount++;
        }

        // New positions available, get the new quality at every node.
        if( check_element_inversion) {
            oldNodeQuality = newNodeQuality;
            getQuality();
        }

        // Update the nodes positions ..
        for( const JNodePtr &vertex: freenodes) update(vertex);

        if( applyTaubinSteps ) {
            for( const JNodePtr &vertex: freenodes) {
                JNode::getRelations( vertex, nodeneighs );
                int err = atomicOp( vertex, nodeneighs, mu);
                if( err == 0) ncount++;
            }
            for( const JNodePtr &vertex: freenodes) update(vertex);
        }

        if( ncount == 0) break;

        maxDiff = sqrt(maxDiff);
        if( maxDiff < tolerance) break;

    }
}
/////////////////////////////////////////////////////////////////////////////////////

int JLaplaceMeshSmoother ::smooth_with_primal_nodes()
{
    int err = 1;
    // Constrainted nodes will not move ..
    size_t numnodes = mesh->getSize(0);
    JNodeSequence freenodes = mesh->getNodes();

/*
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            if( !vertex->hasAttribute("Constraint")) {
                freenodes.push_back(vertex);
            }
        }
    }
    cout << "Smooth nodes " << freenodes.size() << endl;
*/

    smooth( freenodes );
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

int JLaplaceMeshSmoother ::smooth_with_dual_nodes()
{
    if( mesh == nullptr) return 1;

    JNodeSequence freenodes;
    size_t numnodes = mesh->getSize(0);

    if( mesh->getAdjTable(0,2) == 0)
        mesh->buildRelations(0,2);

    JFaceSequence faceneighs;
    JNodeSequence dualneighs;

    JNodePtr dnode = nullptr;
    Point3D xyz;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            if( !vertex->hasAttribute("Constraint")) {
                freenodes.push_back(vertex);
                JNode::getRelations(vertex, faceneighs);
                dualneighs.clear();
                for( size_t j = 0; j < faceneighs.size(); j++) {
                    int err = faceneighs[j]->getAttribute("DualNode", dnode);
                    if( err ) {
                        dnode = JNode::newObject();
                        faceneighs[j]->setAttribute("DualNode", dnode);
                    }
                    JFaceGeometry::getDualPosition(faceneighs[j], xyz);
                    dnode->setXYZCoords(xyz);
                    dualneighs.push_back(dnode);
                }
                vertex->setAttribute("DualNeighs", dualneighs);
            }
        }
    }

    numnodes = freenodes.size();

    for( int j = 0; j < numIterations; j++) {
        size_t ncount = 0;
        maxDiff = 0.0;
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vertex = freenodes[i];
            vertex->getAttribute("DualNeighs", dualneighs);
            int stat = atomicOp(vertex, dualneighs, lambda);
            if( stat == 0) ncount++;
        }
        for( size_t i = 0; i < numnodes; i++) update( freenodes[i] );

        if( applyTaubinSteps ) {
            for( size_t i = 0; i < numnodes; i++) {
                const JNodePtr &vertex = freenodes[i];
                vertex->getAttribute("DualNeighs", dualneighs);
                int stat = atomicOp( vertex, dualneighs, mu);
                if( stat == 0) ncount++;
            }
            for( size_t i = 0; i < numnodes; i++) update( freenodes[i] );
        }

        maxDiff = sqrt(maxDiff);
        if( ncount == 0) return 1;
    }

    mesh->deleteFaceAttribute("DualNode");
    mesh->deleteNodeAttribute("DualNeighs");

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////
void JLaplaceMeshSmoother :: getQuality()
{
    if( mesh == nullptr) return;

    vector<double>  x, y;
    size_t ncount = 0;

    if( topDim == 2) {
        size_t numfaces = mesh->getSize(2);
        faceQuality.resize(numfaces);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            int np = face->getSize(0);
            x.resize(np);
            y.resize(np);
            for( int j = 0; j < np; j++) {
                size_t id  = face->getNodeAt(j)->getID();
                x[j] = newPos[id][0];
                y[j] = newPos[id][1];
            }
            faceQuality[i]  = JGeometry::getSignedArea( &x[0], &y[0], np);
            if( faceQuality[i] < 0.0) ncount++;
        }

        JFaceSequence faceneighs;
        size_t numnodes = mesh->getSize(0);
        newNodeQuality.resize(numnodes);

        vector<double> minQuality;
        for( size_t i = 0; i <  numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            JNode::getRelations(vtx, faceneighs);
            int nfaces = faceneighs.size();
            assert( nfaces );
            minQuality.resize( nfaces );
            for( int  j = 0; j < nfaces; j++) {
                size_t fid = faceneighs[j]->getID();
                minQuality[j]  = faceQuality[fid];
            }
            newNodeQuality[i] = *boost::min_element(minQuality);
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////

int JLaplaceMeshSmoother :: smoothAll()
{
    if( mesh == nullptr) return 1;

    addBoundaryConstraints(mesh);

    int err;
    if( numericalMethod == EXPLICIT_METHOD  ) {
        err = smooth_with_primal_nodes();
    }

//   useGPtoolBox();

    return err;
}
/////////////////////////////////////////////////////////////////////////////////
int JLaplaceMeshSmoother :: useGPtoolBox()
{
#ifdef USE_IGL
    JMeshEigenMatrix mat;
    mat.setMesh(mesh);
    Eigen::MatrixXd V = mat.getNodeMatrix();
    Eigen::MatrixXi T = mat.getElementMatrix();

    ostringstream oss;
    oss << "[U,~] = laplacian_smooth(V,T,b,";
    if( numericalMethod == EXPLICIT_METHOD)
        oss << "'explicit',";
    else
        oss << "'implicit',";
    if( edgeWeightID  == COMBINATORIAL)
        oss << "'uniform'";
    if( edgeWeightID == COTANGENT)
        oss << "'cotan'";

    oss << ")";

    vector<int> vb;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            if( vtx->hasAttribute("Constraint"))
                vb.push_back(vtx->getID());
        }
    }

    Eigen::MatrixXi b;
    int nrows = vb.size();

    b.resize(nrows, 1);
    for( size_t i = 0; i < vb.size(); i++) {
        b.coeffRef(i,0) = vb[i];
    }

    string cmd = oss.str();

    Engine *engine;
    igl::matlab::mlinit(&engine);
    igl::matlab::mlsetmatrix(&engine,"V",V);
    igl::matlab::mlsetmatrix(&engine,"T",T);
    igl::matlab::mlsetmatrix(&engine,"b",b);
    igl::matlab::mleval(&engine, cmd);

    Eigen::MatrixXd U;
    igl::matlab::mlgetmatrix(&engine,"U",U);
    igl::matlab::mlclose(&engine);

    size_t index = 0;
    Point3D xyz;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            xyz[0] = U.coeff(index,0);
            xyz[1] = U.coeff(index,1);
            xyz[2] = U.coeff(index,2);
            cout << V.coeff(index,0) << " " << V.coeff(index,1) << V.coeff(index,2) << endl;
            cout << U.coeff(index,0) << " " << U.coeff(index,1) << U.coeff(index,2) << endl;
            cout << "**************" << endl;
            vtx->setXYZCoords(xyz);
            index++;
        }
    }
    return 0;
#endif
}

/////////////////////////////////////////////////////////////////////////////////
