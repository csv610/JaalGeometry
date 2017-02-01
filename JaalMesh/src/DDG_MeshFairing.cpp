#include "DDG_MeshFairing.hpp"

using namespace DDG;

void DDG::JMeshFairing :: setMesh( const JMeshPtr &m)
{
    jmesh = m;
    if( jmesh == nullptr) return;

    JMeshOBJExporter mexp;
    mexp.writeFile( jmesh, "tmp.obj");

    mesh.read("tmp.obj");
    x = DenseMatrix<Real>( mesh.vertices.size(), 3 );
    timeStep = 0.0001;
}

////////////////////////////////////////////////////////////////////
void DDG::JMeshFairing :: nextStep()
{
    DDG::SparseMatrix<Real> star0;
    DDG::SparseMatrix<Real> star1;
    DDG::SparseMatrix<Real> d0;

    HodgeStar0Form<Real>::build( mesh, star0 );
    HodgeStar1Form<Real>::build( mesh, star1 );
    ExteriorDerivative0Form<Real>::build( mesh, d0 );

    cout << "Current Time Step " << timeStep << endl;
    SparseMatrix<Real> L = d0.transpose() * star1 * d0;
    SparseMatrix<Real> A = star0 + Real(timeStep) * L;

    getPositions();
    DenseMatrix<Real> rhs = star0 * x;

    solvePositiveDefinite(A, x, rhs);

    setPositions();
    updateMesh();
}

////////////////////////////////////////////////////////////////////

void DDG::JMeshFairing :: getPositions()
{
    for ( VertexCIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++)
    {
        for( int i = 0; i < 3; ++i)
            x(v->index, i) = v->position[i];
    }
}

////////////////////////////////////////////////////////////////////
void DDG::JMeshFairing :: setPositions()
{
    for ( VertexIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++)
    {
        v->position = Vector(x(v->index, 0),
                             x(v->index, 1),
                             x(v->index, 2));
    }
}
////////////////////////////////////////////////////////////////////
void DDG::JMeshFairing :: updateMesh()
{
    if( jmesh == nullptr) return;

    cout << "Update mesh " << endl;

    Point3D xyz;
    for ( VertexCIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++)
    {
        xyz[0] = x(v->index,0);
        xyz[1] = x(v->index,1);
        xyz[2] = x(v->index,2);
        const JNodePtr &vtx = jmesh->getNodeAt(v->index);
        vtx->setXYZCoords(xyz);
    }
}
////////////////////////////////////////////////////////////////////
