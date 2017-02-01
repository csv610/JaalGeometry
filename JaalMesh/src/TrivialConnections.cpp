#include "TrivialConnections.hpp"

using namespace DDG;

#ifdef CSV

////////////////////////////////////////////////////////////////////////////////
int JTrivialConnections :: solve()
{
    bool ok = checkGaussBonnet();
    if( not ok )
    {
        std::cout << "Gauss-Bonnet thm does not hold" << std::endl;
        return false;
    }

    int t0 = clock();
    solveForTrivialHolonomy();
    int t1 = clock();
    cout << "[trivial] time: " << seconds( t0, t1 ) << "s" << "\n";

    t0 = clock();
    solveForNonTrivialHolonomy();
    t1 = clock();
    cout << "[nontrivial] time: " << seconds( t0, t1 ) << "s" << "\n";

    return true;
}
////////////////////////////////////////////////////////////////////////////////

bool JTrivialConnections :: checkGaussBonnet()
{
    // vertex singularity
    int k = 0;
    for(VertexCIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v++)
    {
        k += v->singularity;
    }

    cout << "#Num Singular " << k << endl;

    // generator singularity
    // TODO: include singularity for all generators
    k += firstGeneratorIndex;

    return ( mesh.getEulerCharacteristicNumber() == k );
}

////////////////////////////////////////////////////////////////////////////////

void JTrivialConnections :: solveForTrivialHolonomy()
{
    // Neumann boundary condition => prescribing geodesic curvature
    // For now, keeping original geodesic curvature
    DenseMatrix<Real> b( mesh.vertices.size() );
    for(VertexIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v++)
    {
        double value = 0.0;
        if( not v->onBoundary() )
        {
            value -= ( 2. * M_PI - v->theta() );
            value += 2. * M_PI * v->singularity;
        }
        b( v->index ) = value;
    }

    DenseMatrix<Real> u( mesh.vertices.size() );
    if( b.norm() > 1.0e-8 ) backsolvePositiveDefinite( L, u, b );

    for(VertexIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v++)
    {
        v->potential = u( v->index );
    }
}
////////////////////////////////////////////////////////////////////////////////


void JTrivialConnections :: solveForNonTrivialHolonomy()
{
    unsigned nb = numberHarmonicBases();
    if( nb == 0 ) return;
    harmonicCoefs = std::vector<double>(nb, 0.0);

    DenseMatrix<Real>  b(nb);
    SparseMatrix<Real> H(nb,nb);

    int row = 0;
    bool skipBoundaryLoop = true;
    for(unsigned i = 0; i < generators.size(); ++i)
    {
        const Generator& cycle = generators[i];
        if( skipBoundaryLoop and isBoundaryGenerator(cycle) )
        {
            skipBoundaryLoop = false;
            continue;
        }

        for(unsigned j = 0; j < cycle.size(); ++j)
        {
            HalfEdgeIter he = cycle[j];
            for(unsigned col = 0; col < nb; ++col)
            {
                H(row,col) += he->harmonicBases[col];
            }
        }

        double value = - generatorHolonomy( cycle );
        if( row == 0 )
        {
            value += 2.0 * M_PI * firstGeneratorIndex ;
        }
        b(row) = value;
        row++;
    }

    /*
                 DenseMatrix<Real> x(nb);
               if( b.norm() > 1.0e-8 ) solve(H, x, b);

                 for(unsigned i = 0; i < nb; ++i)
                    harmonicCoefs[i] = x(i);
    */
}


bool JTrivialConnections :: isBoundaryGenerator(const Generator& cycle) const
{
    if( cycle.size() == 0 ) return false;
    return ( cycle[0]->vertex->onBoundary() or
             cycle[0]->flip->vertex->onBoundary() );
}
//////////////////////////////////////////////////////////////////////////////

unsigned JTrivialConnections :: numberHarmonicBases() const
{

    unsigned nb = 0;
    bool skipBoundaryLoop = true;
    for(unsigned i = 0; i < generators.size(); ++i)
    {
        const Generator& cycle = generators[i];
        if( skipBoundaryLoop and isBoundaryGenerator(cycle) )
        {
            skipBoundaryLoop = false;
            continue;
        }
        nb++;
    }
    return nb;
}
//////////////////////////////////////////////////////////////////////////////

void JTrivialConnections :: faceFrame(HalfEdgeIter h, Vector& a, Vector& b) const
{
    if( h->onBoundary ) return;

    VertexIter v0 = h->vertex;
    VertexIter v1 = h->next->vertex;

    a = ( v1->position - v0->position ).unit();
    Vector n = h->face->normal();
    b = cross( n, a );
}
//////////////////////////////////////////////////////////////////////////////

double JTrivialConnections :: parallelTransport(HalfEdgeIter h) const
{
    if( h->onBoundary or h->flip->onBoundary ) return 0.0;

    VertexIter v0 = h->vertex;
    VertexIter v1 = h->next->vertex;
    Vector e = v1->position - v0->position;

    Vector aL, bL;
    faceFrame( h->face->he, aL, bL );

    Vector aR, bR;
    faceFrame( h->flip->face->he, aR, bR );

    double deltaL = atan2( dot(e,bL), dot(e,aL) );
    double deltaR = atan2( dot(e,bR), dot(e,aR) );
    return (deltaL - deltaR);
}
//////////////////////////////////////////////////////////////////////////////

double JTrivialConnections :: connectionOneForm(HalfEdgeIter h) const
{
    double angle = 0.0;

    // coclosed term
    double star1 = 0.5 * ( h->cotan() + h->flip->cotan() );
    double u0 = h->flip->vertex->potential;
    double u1 = h->vertex->potential;
    angle += star1*(u1 - u0);

    // harmonic term
    for(unsigned k = 0; k < h->harmonicBases.size(); ++k)
        angle += this->harmonicCoefs[k] * h->harmonicBases[k];

    return angle;
}
//////////////////////////////////////////////////////////////////////////////

double JTrivialConnections :: vertexHolonomy(VertexIter vertex) const
{
    double sum = 2.0*M_PI;
    HalfEdgeIter h = vertex->he;
    do
    {
        sum += parallelTransport(h);
        sum += connectionOneForm(h);
        h = h->flip->next;
    }
    while( h != vertex->he );
    return sum;
}
//////////////////////////////////////////////////////////////////////////////

double JTrivialConnections :: generatorHolonomy(const Generator& cycle)
{
    double sum = 0.0;
    if( cycle.empty() ) return sum;

    for(unsigned k = 0; k < cycle.size(); ++k)
    {
        HalfEdgeIter h = cycle[k];
        sum += parallelTransport(h);
        sum += connectionOneForm(h);
    }

    while( sum <  0.0      ) sum += 2.0*M_PI;
    while( sum >= 2.0*M_PI ) sum -= 2.0*M_PI;
    return sum;
}
//////////////////////////////////////////////////////////////////////////////

void JTrivialConnections :: init()
{
    // Laplacian with Neumann boundary condition
    SparseMatrix<Real> star0, star1, d0, Delta;
    HodgeStar0Form<Real>::build( mesh, star0 );
    HodgeStar1Form<Real>::build( mesh, star1 );
    ExteriorDerivative0Form<Real>::build( mesh, d0 );
    Delta = d0.transpose() * star1 * d0;

    // make L positive-definite
    Delta += Real(1.0e-8)*star0;

    // pre-factorize
    L.build(Delta);

    // generators
    int t0 = clock();
    TreeCotree tct;
    tct.build( mesh );
    int t1 = clock();
    cout << "[generators] time: " << seconds( t0, t1 ) << "s" << "\n";

    unsigned nb = this->numberHarmonicBases();
    if( nb == 0 ) return;

    // harmonic coefs
    harmonicCoefs = std::vector<double>( nb, 0.0 );

    // harmonic bases
    t0 = clock();
    HarmonicBases bases;
    bases.compute( mesh );
    t1 = clock();
    cout << "[harmonic] time: " << seconds( t0, t1 ) << "s" << "\n";
}
#endif
