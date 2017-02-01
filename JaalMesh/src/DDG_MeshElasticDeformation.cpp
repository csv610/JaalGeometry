#include "DDG_MeshElasticDeformation.hpp"

void DDG::JMeshElasticDeformation :: setSource( const JMeshPtr &m)
{
    srcMesh = m;
    if( srcMesh == nullptr) return;

    JMeshOBJExporter mexp;
    mexp.writeFile( srcMesh, "tmp.obj");

    source.read("tmp.obj");
    mesh.read("tmp.obj");
    computeLaplacian(source, src_L);

    interpolatedMesh = srcMesh->deepCopy();
}

////////////////////////////////////////////////////////////////////////////////
void DDG::JMeshElasticDeformation :: setTarget( const JMeshPtr &m)
{
    dstMesh = m;
    if( dstMesh == nullptr) return;

    JMeshOBJExporter mexp;
    mexp.writeFile( dstMesh, "tmp.obj");

    target.read("tmp.obj");
    computeLaplacian(target, tgt_L);
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr DDG::JMeshElasticDeformation :: getInterpolatedMesh(double t)
{
    assert( t >= 0.0 && t <= 1.0);

    double tolerance = 1.0E-08;

    SparseMatrix<Complex> L  = Complex(1.-t)*src_L   + Complex(t)*tgt_L;
    SparseFactor<Complex> LL;
    LL.build( L );

    DenseMatrix<Complex> src_angle, tgt_angle;
    DenseMatrix<Complex> src_rhs, tgt_rhs;
    DenseMatrix<Complex> rhs, x;

    // Init ....
    VertexCIter sv = source.vertices.begin();
    VertexCIter tv = target.vertices.begin();
    VertexIter  v  = mesh.vertices.begin();
    for( ; v != mesh.vertices.end(); v++, tv++, sv++)
    {
        v->position = (1.-t)*sv->position + t*tv->position;
    }

    for( int iter = 0; iter < numIters; iter++ )
    {
        // extract rotation from deformation gradient
        computeRotation(source, mesh, src_angle);
        computeRotation(target, mesh, tgt_angle);

        // optimize vertices' position
        computeDivergence(source, src_angle, src_rhs);
        computeDivergence(target, tgt_angle, tgt_rhs);

        rhs = Complex(1.-t)*src_rhs + Complex(t)*tgt_rhs;

        get2DPositions(mesh, x);

        double res = residual(L, x, rhs);
        if (res < tolerance) break;

        backsolvePositiveDefinite(LL, x, rhs);
        x.removeMean();
        assign2DPositions(x, mesh);
    }

    return interpolatedMesh;
}

////////////////////////////////////////////////////////////////////////////////

void DDG::JMeshElasticDeformation :: computeRotation(const Mesh& meshA,
        const Mesh& meshB, DenseMatrix<Complex>& angle) const
{
    angle = DenseMatrix<Complex>(meshA.faces.size(),1);
    FaceCIter fA = meshA.faces.begin();
    FaceCIter fB = meshB.faces.begin();

    DenseMatrix<Real> A0(2,2);
    DenseMatrix<Real> A1(2,2);

    for( ; fA != meshA.faces.end(); fA++, fB++)
    {
        // edge vectors from A
        VertexIter v0 = fA->he->vertex;
        VertexIter v1 = fA->he->next->vertex;
        VertexIter v2 = fA->he->next->next->vertex;

        Vector e01 = v1->position - v0->position;
        Vector e12 = v2->position - v1->position;

        A0(0,0) = e01.x;
        A0(0,1) = e12.x;
        A0(1,0) = e01.y;
        A0(1,1) = e12.y;

        // edge vectors from B
        v0 = fB->he->vertex;
        v1 = fB->he->next->vertex;
        v2 = fB->he->next->next->vertex;

        e01 = v1->position - v0->position;
        e12 = v2->position - v1->position;

        A1(0,0) = e01.x;
        A1(0,1) = e12.x;
        A1(1,0) = e01.y;
        A1(1,1) = e12.y;

        // ortho(A1*inv(A0))
        DenseMatrix<Real> M = A1 * PolarDecomposition2x2::invert(A0);
        Complex z = PolarDecomposition2x2::extractOrthogonalPart(M);
        angle(fA->index,0) = z;
    }
}

////////////////////////////////////////////////////////////////////////////////

void DDG::JMeshElasticDeformation :: computeLaplacian(const Mesh& mesh,
        SparseMatrix<Complex>& L) const
{
    SparseMatrix<Complex> star0;
    SparseMatrix<Complex> star1;
    HodgeStar0Form<Complex>::build( mesh, star0 );
    HodgeStar1Form<Complex>::build( mesh, star1 );

    SparseMatrix<Complex> d0;
    ExteriorDerivative0Form<Complex>::build( mesh, d0 );

    L = d0.transpose() * star1 * d0;
    L += Complex(1e-8)*star0;
}

////////////////////////////////////////////////////////////////////////////////
void DDG::JMeshElasticDeformation :: computeDivergence(const Mesh& mesh,
        const DenseMatrix<Complex>& angle,
        DenseMatrix<Complex>& rhs) const
{
    rhs = DenseMatrix<Complex>(mesh.vertices.size(),1);
    for( VertexCIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++ )
    {
        Complex sum;
        Vector pi = v->position;
        HalfEdgeIter he = v->he;
        do
        {
            if( not he->onBoundary )
            {
                Vector area_grad = 0.5*he->next->rotatedEdge();
                Complex sum_ijk( area_grad.x, area_grad.y );
                Complex Rijk = angle(he->face->index, 0);
                sum += Rijk * sum_ijk;
            }
            he = he->flip->next;
        }
        while( he != v->he );
        rhs(v->index,0) = sum;
    }
}
////////////////////////////////////////////////////////////////////////////////

void DDG::JMeshElasticDeformation :: assign2DPositions(const DenseMatrix<Complex>& x, Mesh& mesh)
{
    Point3D xyz;
    for( VertexIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++ )
    {
        Complex xi = x(v->index,0);
        v->position = Vector(xi.re, xi.im, 0.);
        const JNodePtr &vtx = interpolatedMesh->getNodeAt(v->index);
        xyz[0] = xi.re;
        xyz[1] = xi.im;
        xyz[2] = 0.0;
        vtx->setXYZCoords(xyz);
    }

}
////////////////////////////////////////////////////////////////////////////////

void DDG::JMeshElasticDeformation :: get2DPositions(const Mesh& mesh, DenseMatrix<Complex>& x) const
{
    x = DenseMatrix<Complex>(mesh.vertices.size(),1);
    for( VertexCIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++ )
    {
        Complex z(v->position.x, v->position.y);
        x(v->index,0) = z;
    }
}
////////////////////////////////////////////////////////////////////////////////
