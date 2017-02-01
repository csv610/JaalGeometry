#include "ShapeOpt.hpp"

////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    size_t numnodes = mesh->getActiveSize(0);
    pCloud.resize(3,numnodes);

    numnodes = mesh->getSize(0);
    size_t index = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            const Point3D  &xyz = vtx->getXYZCoords();
            pCloud(0,index) = xyz[0];
            pCloud(1,index) = xyz[1];
            pCloud(2,index) = xyz[2];
            index++;
        }
    }
    solver.setPoints(pCloud);
}

////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addAreaConstraint( const JFaceSequence &faces, double weight)
{
    double avgArea = 0.0;
    for( const JFacePtr &face : faces)
        avgArea += JFaceGeometry::getSignedArea(face);
    avgArea /= ( double) faces.size();

    vector<int> vid;
    for( const JFacePtr &face : faces) {
        double area = JFaceGeometry::getSignedArea(face);
        double factor = avgArea;
        double minval = 0.99*factor;
        double maxval = 1.01*factor;
        vid.clear();
        int nn = face->getSize(0);
        for( int i = 0; i < nn; i++)
            vid.push_back( face->getNodeAt(i)->getID() );
        auto c = std::make_shared<ShapeOp::AreaConstraint>(vid, weight, solver.getPoints(), minval, maxval);
        cdb["Area"].push_back(c);
    }
}
////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addPlaneConstraint( const JNodeSequence &nodes, double weight)
{
    vector<int> vid;
    for( const JNodePtr & vtx : nodes) {
        vid.push_back(vtx->getID() );
    }
    auto c = std::make_shared<ShapeOp::PlaneConstraint>(vid, weight, solver.getPoints() );
    cdb["Plane"].push_back(c);
}
////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addCloseConstraint( const JNodeSequence &nodes, double weight)
{
    vector<int> vid;
    ShapeOp::Vector3  pos;
    for( const JNodePtr &vtx : nodes) {
        vid.clear();
        vid.push_back(vtx->getID() );
        auto c = std::make_shared<ShapeOp::ClosenessConstraint>(vid, weight, solver.getPoints() );
        const Point3D &xyz = vtx->getXYZCoords();
        pos[0] = xyz[0];
        pos[1] = xyz[1];
        pos[2] = xyz[2];
        c->setPosition(pos);
        cdb["Closeness"].push_back(c);
    }
}
////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addCircleConstraint( const JFaceSequence &faces, double weight)
{
    vector<int> vid;
    for( const JFacePtr &face : faces) {
        int nn = face->getSize(0);
        vid.clear();
        for( int i = 0; i < nn; i++)
            vid.push_back( face->getNodeAt(i)->getID() );
        auto c = std::make_shared<ShapeOp::CircleConstraint>(vid, weight, solver.getPoints() );
        cdb["Circle"].push_back(c);
    }
}

////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addEdgeStrainConstraint( const JEdgeSequence &edges, double len, double weight)
{
    double infy = std::numeric_limits<double>::max();
    vector<int> vid(2);
    ShapeOp::EdgeStrainConstraint *c;
    for( const JEdgePtr &edge : edges) {
        vid[0] = edge->getNodeAt(0)->getID();
        vid[1] = edge->getNodeAt(1)->getID();
        auto c = std::make_shared<ShapeOp::EdgeStrainConstraint>(vid, weight, solver.getPoints() );
        c->setEdgeLength(len);
        cdb["EdgeStrain"].push_back(c);
    }
}

////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addRectangleConstraints(double weight)
{
    vector<int> vid;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int nn = face->getSize(0);
            if(nn == 4 ) {
                vid.clear();
                for( int j = 0; j < 4; j++)
                    vid.push_back( face->getNodeAt(j)->getID() );
                auto c = std::make_shared<ShapeOp::RectangleConstraint>(vid, weight, solver.getPoints() );
                cdb["Rectangle"].push_back(c);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addParallelogramConstraints(double weight)
{
    vector<int> vid;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int nn = face->getSize(0);
            if(nn == 4 ) {
                vid.clear();
                for( int j = 0; j < 4; j++)
                    vid.push_back( face->getNodeAt(i)->getID() );
                auto c = std::make_shared<ShapeOp::ParallelogramConstraint>(vid, weight, solver.getPoints() );
                cdb["Parallelogram"].push_back(c);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
/*
void JShapeOptimizer :: addLengthConstraint()
{
    double elen = mesh->getGeometry()->getMeanBoundaryEdgeLength();
    JEdgeSequence edges;
    mesh->getTopology()->getInternal(edges);
    addLengthConstraint(edges, elen);
}

////////////////////////////////////////////////////////////////////////////////
*/

void JShapeOptimizer :: addBoundaryConstraints(double weight)
{
    if( mesh == nullptr) return;
    JNodeSequence bdnodes;
    mesh->getTopology()->searchBoundary();
    mesh->getTopology()->getBoundary( bdnodes);
    addCloseConstraint( bdnodes);
}
////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addCircleConstraints(double weight)
{
    JFaceSequence faces = mesh->getFaces();
    addCircleConstraint(faces);
}

////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addAreaConstraints(double weight)
{
    JFaceSequence faces = mesh->getFaces();
    addAreaConstraint(faces);

    JNodeSequence bdnodes;
    mesh->getTopology()->searchBoundary();
    mesh->getTopology()->getBoundary( bdnodes);

    addCloseConstraint( bdnodes);
}

////////////////////////////////////////////////////////////////////////////////
void JShapeOptimizer :: addPlaneConstraints(double weight)
{
    JNodeSequence nodes = mesh->getNodes();
    addPlaneConstraint(nodes);
}

////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: addLaplaceConstraints(double weight)
{
    mesh->buildRelations(0,0);
    JNodeSequence nodes = mesh->getNodes();
    JNodeSequence vneighs;

    vector<int> vid;
    for( const JNodePtr &vi : nodes) {
        vid.clear();
        if (!vi->isBoundary() ) {
            vid.push_back(vi->getID() );
            JNode::getRelations(vi, vneighs);
            for( const JNodePtr &vj : vneighs) {
                vid.push_back(vj->getID() );
            }
            auto c = std::make_shared<ShapeOp::UniformLaplacianConstraint>(vid, weight, solver.getPoints(), 1);
            cdb["Laplace"].push_back(c);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

void JShapeOptimizer :: solve()
{
    solver.initialize(1);
    solver.setTimeStep(0.1);
    solver.solve(1000);

    pCloud = solver.getPoints();

    size_t numnodes = mesh->getSize(0);

    size_t index = 0;
    Point3D xyz;

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            xyz[0] = pCloud(0,index);
            xyz[1] = pCloud(1,index);
            xyz[2] = pCloud(2,index);
            index++;
            vtx->setXYZCoords(xyz);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
