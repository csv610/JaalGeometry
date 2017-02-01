#pragma once

#include "Mesh.hpp"

#include "DDG_Mesh.hpp"
#include "DDG_SparseMatrix.hpp"
#include "DDG_DenseMatrix.hpp"
#include "DDG_Direction.hpp"
#include "DDG_DiscreteExteriorCalculus.hpp"
#include "DDG_HarmonicBases.hpp"
#include "DDG_Quaternion.hpp"
#include "DDG_Real.hpp"
#include "DDG_TreeCotree.hpp"
#include "DDG_Utility.hpp"

using namespace DDG;

class JTrivialConnections
{
    typedef std::vector<HalfEdgeIter> Generator;
public:
    void setMesh(const JMeshPtr &m);
    int  solve();
private:
    DDG::Mesh mesh;
    DDG::DirectionField field;
    DDG::SparseFactor<Real> L;
    // pre-factorization of Laplacian

    std::vector<Generator> generators;
    // non-contractible loops

    std::vector<double> harmonicCoefs;
    // coefficients for linear combination of harmonic bases

    double firstGeneratorIndex;

    void init();
    bool checkGaussBonnet();
    void solveForTrivialHolonomy();
    void solveForNonTrivialHolonomy();

    bool isBoundaryGenerator(const Generator& cycle) const;
    // returns true if generator is a boundary loop

    unsigned numberHarmonicBases() const;
    // returns 2g + (m-1) where g = genus and m = number of boundary loops

    double connectionOneForm(HalfEdgeIter h) const;
    // returns rotation angle by crossing h

    double parallelTransport(HalfEdgeIter h) const;
    // returns rotation around h->vertex

    void faceFrame(HalfEdgeIter h, Vector& a, Vector& b) const;
    // returns unit edge vector parallel to h and
    // its rotation by \pi/2 around the face's normal

    double vertexHolonomy(VertexIter vertex) const;
    // returns defect angle by rotating around vertex

    double generatorHolonomy(const Generator& cycle);
    // returns defect angle by rotating around cycle
};
