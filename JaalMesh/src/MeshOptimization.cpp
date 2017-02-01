#include "MeshOptimization.hpp"
#include "AllTriMeshGenerator.hpp"
#include "AllTetMeshGenerator.hpp"

#ifdef USE_IGL
#include <igl/matlab/matlabinterface.h>
#endif

using namespace Mesquite2;

///////////////////////////////////////////////////////////////////////////////////
JLogger* JMeshOptimization::logger = JLogger::getInstance();

void JMeshOptimization::setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    meshDim = mesh->getTopology()->getDimension();
    geomDim = mesh->getGeometry()->getDimension();
}

///////////////////////////////////////////////////////////////////////////////////

void JMeshGeometricOptimization::addBoundaryConstraints( const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;

    JNodeSequence nodes;
    mesh->getTopology()->getBoundary(nodes);

    int gid = 1;
    if( preserve_boundary ) {
        for( const JNodePtr &vtx : nodes) vtx->setAttribute("Constraint", gid);
    } else {
        for( const JNodePtr &vtx : nodes) vtx->deleteAttribute("Constraint");
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void JMeshGeometricOptimization::addConstraints( const JNodeSequence &nodes, int groupID)
{
    if( mesh == nullptr) return;
    for( const JNodePtr &vtx : nodes) vtx->setAttribute("Constraint", groupID);
}
/////////////////////////////////////////////////////////////////////////////////////////////////

void JMeshGeometricOptimization::removeAllConstraints()
{
    if( mesh == nullptr) return;
    mesh->deleteNodeAttribute("Constraint");
}
/////////////////////////////////////////////////////////////////////////////////////////////////


int JMeshNonlinearOptimization::untangle()
{
    if( mesh == nullptr) return 1;

    int dim = mesh->getTopology()->getDimension();
    if( dim == 3) {
        cout << "Warning: Currently untangling only for 2D meshes " << endl;
        return 1;
    }

    JMeshPtr jmesh = mesh->deepCopy();

    int topo = mesh->getTopology()->getElementsType(2);
    if( topo == JFace::QUADRILATERAL )
        jmesh = AllTriMeshGenerator::getFromQuadMesh(jmesh,4);

    addBoundaryConstraints(jmesh);

    Mesquite2::ArrayMesh *optmesh = nullptr;
    optmesh = JMeshNonlinearOptimization::jaal_to_mesquite(jmesh);

    if( optmesh == nullptr) return 1;

    int err = untangle( optmesh );
    if( !err) mesh->getGeometry()->setCoordsArray(vCoords, l2g);

    delete optmesh;

    vCoords.clear();
    l2g.clear();

    return err;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

int JMeshNonlinearOptimization::shapeOptimize()
{
    addBoundaryConstraints(mesh);

    MsqError err;

    Mesquite2::ArrayMesh *optmesh = nullptr;
    optmesh = JMeshNonlinearOptimization::jaal_to_mesquite( mesh);
    if( optmesh == nullptr ) return 1;

    MeshDomainAssoc *mesh_and_domain = nullptr;

    PlanarDomain mesh_plane(PlanarDomain::XY);

    if( meshDim == 2 && geomDim == 2)
        mesh_and_domain = new MeshDomainAssoc( optmesh, &mesh_plane);
    else
        mesh_and_domain = new MeshDomainAssoc( optmesh, nullptr);

    IdealWeightInverseMeanRatio  *q1 = new IdealWeightInverseMeanRatio;

    ShapeImprover smoother;
    smoother.quality_assessor().add_quality_assessment(q1);

    smoother.run_instructions(mesh_and_domain, err);

    if (!err) mesh->getGeometry()->setCoordsArray(vCoords, l2g);

    if( mesh_and_domain) delete mesh_and_domain;
    if( optmesh ) delete optmesh;

    return 0;
}


int JMeshNonlinearOptimization::execute(const JMeshPtr &mesh)
{
    addBoundaryConstraints(mesh);

    MsqError err;

    Mesquite2::ArrayMesh *optmesh = nullptr;
    optmesh = JMeshNonlinearOptimization::jaal_to_mesquite( mesh);
    if( optmesh == nullptr ) return 1;

    MeshDomainAssoc *mesh_and_domain = nullptr;

    PlanarDomain mesh_plane(PlanarDomain::XY);
    if( meshDim == 2 && geomDim == 2) {
        mesh_and_domain = new MeshDomainAssoc( optmesh, &mesh_plane);
    } else
        mesh_and_domain = new MeshDomainAssoc( optmesh, nullptr);

    boost::scoped_ptr<QualityMetric> qMetric;
    boost::scoped_ptr<QualityAssessor> qAssessor;

    IdealWeightInverseMeanRatio  *q1 = nullptr;
    IdealWeightMeanRatio         *q2 = nullptr;
    ConditionNumberQualityMetric *q3 = nullptr;
    EdgeLengthQualityMetric      *q4 = nullptr;

    switch( qualMetric) {
    case INVERSE_MEAN_RATIO:
        q1 = new IdealWeightInverseMeanRatio;
        q1->set_averaging_method(average_method, err);
        qMetric.reset( q1);
        qAssessor.reset( new QualityAssessor(q1));
        break;
    case MEAN_RATIO:
        q2 = new IdealWeightMeanRatio;
        q2->set_averaging_method(average_method, err);
        qMetric.reset( q2 );
        qAssessor.reset( new QualityAssessor(q2));
        break;
    case CONDITION_NUMBER:
        q3 = new ConditionNumberQualityMetric;
        q3->set_averaging_method(average_method, err);
        qMetric.reset( q3 );
        qAssessor.reset( new QualityAssessor(q3));
        break;
    case EDGE_LENGTH:
        q4 = new EdgeLengthQualityMetric;
        q4->set_averaging_method(average_method, err);
        qMetric.reset( q4 );
        qAssessor.reset( new QualityAssessor(q4));
        break;
    }

    /*
        IdealShapeTarget target;
        TShapeB1  mu;
        TQualityMetric metric_0( &target, &mu) ;
        ElementPMeanP metric( 1.0, &metric_0);
        qAssessor->add_quality_assessment( &metric );
    */

    // sets the default objective function template
    boost::scoped_ptr<ObjectiveFunctionTemplate>  obj_func;

    if( norm == LP_NORM)
        obj_func.reset( new LPtoPTemplate( qMetric.get(), normVal, err));

    if( norm == LMEAN_NORM)
        obj_func.reset( new PMeanPTemplate( 2.0, qMetric.get()));

    if( norm == LINF_NORM)
        obj_func.reset( new LInfTemplate( qMetric.get()));

//  LInfTemplate obj_func(&untangle);

    TerminationCriterion tc_outer;

    TerminationCriterion tc_inner;
    tc_inner.add_absolute_gradient_L2_norm(1e-4);
    tc_inner.add_iteration_limit(numIter);

//  tc_inner.write_iterations("opt.dat", err);

    boost::scoped_ptr<QuasiNewton>       qn;
    boost::scoped_ptr<SteepestDescent>   sp;
    boost::scoped_ptr<TrustRegion>       tr;
    boost::scoped_ptr<FeasibleNewton>    fn;
//  boost::scoped_ptr<ConjugateGradient> cg;
    boost::scoped_ptr<SmartLaplacianSmoother> lp;

    Mesquite2::QualityImprover *improver = nullptr;

    switch (algorithm) {
    case STEEPEST_DESCENT:
        sp.reset(new SteepestDescent( obj_func.get() ));
        if( patch_type  == 0)
            sp->use_element_on_vertex_patch();
        else
            sp->use_global_patch();
        improver = sp.get();
        tc_inner.add_iteration_limit(numIter);
        improver->set_inner_termination_criterion(&tc_inner);
        break;
    case QUASI_NEWTON:
        qn.reset(new QuasiNewton( obj_func.get() ));
        if( patch_type  == 0)
            qn->use_element_on_vertex_patch();
        else
            qn->use_global_patch();
        improver = qn.get();
        tc_inner.add_iteration_limit(numIter);
        improver->set_inner_termination_criterion(&tc_inner);
        break;
    case TRUST_REGION:
        tr.reset(new TrustRegion( obj_func.get() ));
        if( patch_type  == 0)
            tr->use_element_on_vertex_patch();
        else
            tr->use_global_patch();
        improver = tr.get();
        tc_inner.add_iteration_limit(numIter);
        improver->set_inner_termination_criterion(&tc_inner);
        break;
    case FEASIBLE_NEWTON:
        fn.reset(new FeasibleNewton( obj_func.get() ));
        if( patch_type  == 0)
            fn->use_element_on_vertex_patch();
        else
            fn->use_global_patch();

        improver = fn.get();
        tc_inner.add_iteration_limit(numIter);
        improver->set_inner_termination_criterion(&tc_inner);
        break;
    case CONJUGATE_GRADIENT:
        /*
                cg.reset(new ConjugateGradient(obj_func.get() ));
                if( patch_type  == 0)
                    cg->use_element_on_vertex_patch();
                else
                    cg->use_global_patch();

                improver = cg.get();
                tc_inner.add_iteration_limit(numIter);
                improver->set_inner_termination_criterion(&tc_inner);
        i*/
        break;
    case SMART_LAPLACIAN:
        lp.reset(new SmartLaplacianSmoother( obj_func.get() ));
        improver = lp.get();
        tc_outer.add_iteration_limit(numIter);
        improver->set_outer_termination_criterion(&tc_outer);
        break;
    default:
        logger->setWarn("Invalid optimization algorithm selected");
    }


    if( improver ) {
        // creates an instruction queue
        InstructionQueue queue;

        queue.add_quality_assessor(qAssessor.get(), err);
        queue.set_master_quality_improver(improver, err);
        queue.add_quality_assessor(qAssessor.get(), err);

        queue.run_instructions(mesh_and_domain, err);

        if (!err) mesh->getGeometry()->setCoordsArray(vCoords, l2g);
    }

    if( mesh_and_domain) delete mesh_and_domain;
    if( optmesh ) delete optmesh;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////

int JMeshNonlinearOptimization::qualityOptimize2D()
{
    int topo = mesh->getTopology()->getElementsType(2);

    if( topo == 0 ) {
        logger->setInfo("Optimizing polygonal mesh");
        JMeshPtr dummyMesh = mesh->deepCopy();
        JMeshPtr trimesh = AllTriMeshGenerator::getFromPolyMesh(dummyMesh);
        trimesh->getTopology()->searchBoundary();
        execute(trimesh);
        size_t numNodes = mesh->getSize(0);
        size_t index = 0;
        for( size_t i = 0; i < numNodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->isActive() ) {
                const JNodePtr &vtri = trimesh->getNodeAt(index++);
                const Point3D &p3d = vtri->getXYZCoords();
                vtx->setXYZCoords(p3d);
            }
        }
        trimesh->deleteAll();
        return 0;
    }

    if( topo == JFace::TRIANGLE ) {
        execute(mesh);
        return 0;
    }

    if( topo == JFace::QUADRILATERAL ) {
        if( !use_simplicial ) {
            execute(mesh);
            return 0;
        }
        logger->setInfo("Optimizing triangle mesh of the quadrilateral mesh: ");
        JMeshPtr dummyMesh = mesh->deepCopy();
        JMeshPtr trimesh = AllTriMeshGenerator::getFromQuadMesh(dummyMesh, 4);
        execute(trimesh);
        size_t numNodes = mesh->getSize(0);
        size_t index = 0;
        for( size_t i = 0; i < numNodes; i++) {
            const JNodePtr &vquad = mesh->getNodeAt(i);
            if( vquad->isActive() ) {
                JNodePtr vtri = trimesh->getNodeAt(index++);
                const Point3D &p3d = vtri->getXYZCoords();
                vquad->setXYZCoords(p3d);
            }
        }
        trimesh->deleteAll();
//      delete trimesh;
//      delete dummyMesh;
        return 0;
    }
    return 1;
}

//////////////////////////////////////////////////////////////////////////////////

int JMeshNonlinearOptimization::qualityOptimize3D()
{
    int topo = mesh->getTopology()->getElementsType(3);

    if( topo == JCell::TETRAHEDRON ) {
        execute(mesh);
        return 0;
    }

    if( topo == JCell::HEXAHEDRON ) {
        if(!use_simplicial) {
            execute(mesh);
            return 0;
        }

        AllTetMeshGenerator alltets;
        JMeshPtr tetmesh = alltets.fromHexMesh(mesh);
        execute(tetmesh);
        tetmesh->deleteCells();
        size_t numNodes = mesh->getSize(0);
        size_t index = 0;
        for( size_t i = 0; i < numNodes; i++) {
            const JNodePtr &vhex = mesh->getNodeAt(i);
            if( vhex->isActive() ) {
                const JNodePtr vtet = tetmesh->getNodeAt(index++);
                const Point3D &p3d = vtet->getXYZCoords();
                vhex->setXYZCoords(p3d);
            }
        }
        return 0;
    }
    return 1;
}
//////////////////////////////////////////////////////////////////////////////////

int JMeshNonlinearOptimization::improveQuality(bool simplices)
{
    assert( mesh );
    use_simplicial = simplices;
    logger->setInfo("Mesh Optimization starts");

    int meshDim = mesh->getTopology()->getDimension();

    if( meshDim == 2) {
        qualityOptimize2D();
        return 0;
    }

    if( meshDim == 3) {
        qualityOptimize3D();
        return 0;
    }

    logger->setWarn("Optimization for mixed elements currently not supported yet");

    return 1;
}


///////////////////////////////////////////////////////////////////////////////
Mesquite2::ArrayMesh* JMeshNonlinearOptimization::jaal_to_mesquite( const JMeshPtr &mesh)
{
    MsqError err;

    unsigned long int numnodes = mesh->getSize(0);
    unsigned long int numfaces = mesh->getSize(2);
    unsigned long int numcells = mesh->getSize(3);

    //
    //////////////////////////////////////////////////////////////////////////////
    // Mesquite works only when the faces are convex, and our mesh may contain some
    // truly bad elements. Therefore, we will skip those faces and mark them fix.
    // Hopefully, later by some other operations ( smoothing, decimation, refinement
    // etc.) those bad elements may get removed.

    // By ignoring some of the faces, the mesh may get disconnected. ( I need to test
    // such cases ).
    //
    // The mesh must be doublets free. ( I found out this after some pains ).
    //////////////////////////////////////////////////////////////////////////////

    size_t index = 0;
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) vertex->setID( index++ );
    }
    size_t noptnodes = index;

    vfixed.resize(noptnodes);

//  #pragma omp parallel for
    size_t ncount = 0;
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        int lid = vertex->getID();
        if( vertex->hasAttribute("Constraint") ) {
            vfixed[lid] = 1;
            ncount++;
        } else {
            vfixed[lid] = 0;
        }
    }

    if( ncount == numnodes ) {
        cout << "Warning: all nodes are constrained; nothing to move  " << endl;
        return nullptr;
    }

    JMeshQuality  mq;
    mq.setMesh(mesh);
    size_t nCount = 0;
    for (size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            double val = mq.getJacobian(face);
            if( val < 0.0) {
                nCount++;
                for( int j = 0; j < face->getSize(0); j++) {
                    int lid = face->getNodeAt(j)->getID();
                    vfixed[lid] = 1;
                }

            }
        }
    }

    if( nCount > numfaces/2) {
        mesh->getTopology()->reverseAll();
    }

    mesh->getGeometry()->getCoordsArray( vCoords, l2g );
    mesh->getTopology()->getNodesArray( vNodes, etopo );

    Mesquite::ArrayMesh *optmesh = nullptr;

    int topo = mesh->getTopology()->getElementsType(meshDim);

    assert( noptnodes == vCoords.size()/3 );
    assert( noptnodes == vfixed.size());

    if( meshDim == 2) {
        size_t noptfaces = 0;
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() ) noptfaces++;
        }
        if (topo == JFace::QUADRILATERAL) {
            optmesh = new Mesquite::ArrayMesh(3, noptnodes, &vCoords[0], &vfixed[0],
                                              noptfaces, Mesquite::QUADRILATERAL, &vNodes[0]);
        }
        if (topo == JFace::TRIANGLE) {
            optmesh = new Mesquite::ArrayMesh(3, noptnodes, &vCoords[0], &vfixed[0],
                                              noptfaces, Mesquite::TRIANGLE, &vNodes[0]);
        }
    }

    if( meshDim == 3) {
        size_t noptcells = 0;
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &c = mesh->getCellAt(i);
            if( c->isActive() ) noptcells++;
        }
        if (topo == JCell::TETRAHEDRON)
            optmesh = new Mesquite::ArrayMesh(3, noptnodes, &vCoords[0], &vfixed[0],
                                              noptcells, Mesquite::TETRAHEDRON, &vNodes[0]);
        if (topo == JCell::HEXAHEDRON) {
            optmesh = new Mesquite::ArrayMesh(3, noptnodes, &vCoords[0], &vfixed[0],
                                              noptcells, Mesquite::HEXAHEDRON, &vNodes[0]);
        }
    }
    return optmesh;
}

///////////////////////////////////////////////////////////////////////////////
int  JCyclicQuadMeshOptimizer :: atomicOp( const JNodePtr &vtx)
{
    if( vtx->isBoundary() ) return 2;

    JNodeSequence vneighs;
    JNode::getRelations(vtx,vneighs);

    if( vneighs.size() != 4 ) return 1;

    bool nearboundary = 0;
    for( const JNodePtr &v : vneighs) {
        if( v->isBoundary() ) {
            nearboundary = 1;
            break;
        }
    }
    if( !nearboundary) return 1;

    // Ref : A New Smoothing Algorithm for Quadrilateral and Hexahedral Meshes
    // Sanjay Kumar Khatri1

    Point3D newpos, xyz, apex;
    apex = vtx->getXYZCoords();

    JFaceSequence faceneighs;
    JNode::getRelations(vtx,faceneighs);

    double xlapCorr = 0.0;
    double ylapCorr = 0.0;
    double zlapCorr = 0.0;
    for( int i = 0; i < 4; i++) {
        JQuadrilateralPtr face = JQuadrilateral::down_cast( faceneighs[i] );
        assert( face );
        const JNodePtr &v =  JQuadrilateral::getDiagonalNode(face, vtx);
        xyz = v->getXYZCoords();
        xlapCorr  += xyz[0] - apex[0];
        ylapCorr  += xyz[1] - apex[1];
        zlapCorr  += xyz[2] - apex[2];
    }

    double xlapStd = 0.0;
    double ylapStd = 0.0;
    double zlapStd = 0.0;
    for( int i = 0; i < 4; i++) {
        const JNodePtr &v =  vneighs[i];
        xyz = v->getXYZCoords();
        xlapStd   += xyz[0] - apex[0];
        ylapStd   += xyz[1] - apex[1];
        zlapStd   += xyz[2] - apex[2];
    }

    newpos    =  vtx->getXYZCoords();
    newpos[0] += +0.25*xlapCorr - 0.50*xlapStd;
    newpos[1] += +0.25*ylapCorr - 0.50*ylapStd;
    newpos[2] += +0.25*zlapCorr - 0.50*zlapStd;

    vtx->setXYZCoords(newpos);
    vtx->setAttribute("TargetPos", newpos);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////
int JCyclicQuadMeshOptimizer :: smoothAll()
{
    if( mesh == nullptr) {
        cout << "Warning: A null mesh object passed in -JCyclicQuadMeshOptimizer" << endl;
        return 1;
    }

    mesh->pruneNodes();
    mesh->buildRelations(0,0);
    mesh->buildRelations(0,2);

    size_t numnodes = mesh->getSize(0);
    if( preserve_boundary) {
        for(size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->isBoundary() ) {
                const Point3D &p = vtx->getXYZCoords();
                vtx->setAttribute("TargetPos", p);
            }
        }
    } else {
        mesh->deleteAttribute("TargetPos");
    }

    for(size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() && !vtx->isBoundary() )
            atomicOp(vtx);
    }

    /*
        JLocallyInjectiveMap lim;
        lim.setMesh(mesh);
        lim.solve();
    */
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JCyclicQuadMeshOptimizer :: smooth( const JNodeSequence &nodes)
{
    if( mesh == nullptr) return 1;

    mesh->pruneNodes();

    mesh->buildRelations(0,0);
    mesh->buildRelations(0,2);

    size_t numnodes = mesh->getSize(0);
    newCoords.resize(numnodes);

    int index = 0;
    for( const JNodePtr &v : nodes)
        newCoords[index++] = v->getXYZCoords();

    for( int iter = 0; iter < numIter; iter++) {
        for( const JNodePtr &vtx : nodes)
            atomicOp(vtx);
        index = 0;
        for( const JNodePtr &vtx : nodes)
            vtx->setXYZCoords( newCoords[index++] );
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JBoundedDistortionOptimizer :: optimize( const JMeshPtr &mesh, double K)
{

#ifdef USE_IGL
    size_t numnodes = mesh->getSize(0);
    size_t numcells = mesh->getSize(3);

    Eigen::MatrixXd V;    // Soring Vertices ...
    Eigen::MatrixXi T;    // Soring Tets

    V.resize(numnodes,3);
    T.resize(numcells,4);

    size_t index = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->setID(index);
            const Point3D  &xyz = vtx->getXYZCoords();
            V.coeffRef(index,0) = xyz[0];
            V.coeffRef(index,1) = xyz[1];
            V.coeffRef(index,2) = xyz[2];
            index++;
        }
    }
    vector<int> conn(4);

    // Notice: Matlab index starts from 1;
    // Lipman's Code expect reverse orientation ...
    index = 0;
    for( size_t i = 0; i < numcells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            for( int j = 0; j < 4; j++)
                conn[j] = cell->getNodeAt(j)->getID();
            T.coeffRef(index,0) = conn[0];
            T.coeffRef(index,1) = conn[2];
            T.coeffRef(index,2) = conn[1];
            T.coeffRef(index,3) = conn[3];
            index++;
        }
    }

    Engine* engine;
    igl::matlab::mlinit(&engine);

    // Send Laplacian matrix to matlab
    igl::matlab::mlsetmatrix(&engine,"V",V);
    igl::matlab::mlsetmatrix(&engine,"T",T);

    // Calling Yaron Lipman's Matlab code ....
    ostringstream oss;
    oss << "[X,~] = improve_mesh(V,T," << K << ")";
    string cmd = oss.str();
    cout << "solve using matlab enginer ... " << endl;
    igl::matlab::mleval(&engine, cmd);

    Eigen::MatrixXd X;
    igl::matlab::mlgetmatrix(&engine,"X", X);

    Point3D xyz;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        xyz[0] = X.coeff(i,0);
        xyz[1] = X.coeff(i,1);
        xyz[2] = X.coeff(i,2);
        if( vtx->isActive() ) vtx->setXYZCoords(xyz);
    }

    igl::matlab::mlclose(&engine);
    return 0;
#endif
}

///////////////////////////////////////////////////////////////////////////////


int JMeshNonlinearOptimization::untangle(Mesquite2::ArrayMesh *optmesh)
{
    if( optmesh == nullptr) return 1;

    /*
        MsqError err;
        MeshDomainAssoc *mesh_and_domain = nullptr;

        PlanarDomain mesh_plane(PlanarDomain::XY);
        if( meshDim == 2 && geomDim == 2) {
            mesh_and_domain = new MeshDomainAssoc( optmesh, &mesh_plane);
        } else {
            mesh_and_domain = new MeshDomainAssoc( optmesh, nullptr);
        }

        InstructionQueue queue1;
        ConditionNumberQualityMetric shape_metric;
        UntangleBetaQualityMetric  untangle(2);
        LInfTemplate obj_func(&untangle);

        ConjugateGradient cg(&obj_func, err);
        if( err ) return 1;
        QualityAssessor  stop_qa  = QualityAssessor( &shape_metric);
        QualityAssessor  stop_qa2 = QualityAssessor( &shape_metric);
        stop_qa.add_quality_assessment(&untangle);
        queue1.add_quality_assessor( &stop_qa, err);
        if( err ) return 1;
        queue1.set_master_quality_improver( &cg, err);
        if( err ) return 1;
        queue1.add_quality_assessor( &stop_qa2, err);
        if( err ) return 1;

        TerminationCriterion tc_inner;
        tc_inner.add_absolute_gradient_L2_norm(1e-4);
        tc_inner.add_iteration_limit(numIter);
        cg.set_inner_termination_criterion(&tc_inner);

        queue1.run_instructions( mesh_and_domain, err);
        if( err ) return 1;
        return 0;

            ShapeImprover smoother;
            IdealWeightInverseMeanRatio extra_metric;
            smoother.quality_assessor().add_quality_assessment(&extra_metric);
            smoother.run_instructions( mesh_and_domain, err );

            return err;
        */

    return 1;

}

