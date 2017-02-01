#pragma once

#include <Eigen/Sparse>
#include "MeshOptimization.hpp"
#include "MeshQuality.hpp"
#include "MeshMatrix.hpp"
#include "SparseMatrix.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

struct EdgeWeight {
    virtual ~EdgeWeight() {}
    virtual string getName() const = 0;
    virtual bool isSymmetric() const = 0;
    virtual double getWeight(const JEdgePtr &edge, int) = 0;
};

struct CombinatorialWeight : public EdgeWeight {
    ~CombinatorialWeight() {}
    string getName() const
    {
        return "Combinatorial";
    }
    bool isSymmetric() const
    {
        return 1;
    }

    double getWeight(const JEdgePtr &, int )
    {
        return 1.0;
    }
};

struct UntangleWeight : public EdgeWeight {
    ~UntangleWeight() {}
    string getName() const
    {
        return "UntangleCombinatorial";
    }
    bool isSymmetric() const
    {
        return 1;
    }

    double getWeight(const JEdgePtr &edge, int )
    {
        if( isTangled(edge) )
            return 1.0;
        else
            return val;
    }
    double val;
private:
    JMeshQuality mq;
    bool isTangled( const JNodePtr &vtx)
    {
        JFaceSequence neighs;
        JNode::getRelations(vtx, neighs);
        double val;
        for( const JFacePtr &face : neighs) {
            val = mq.getJacobian(face);
            if( val < 0.0) return 1;
        }
        return 0;
    }


    bool isTangled( const JEdgePtr &edge)
    {
        if( isTangled( edge->getNodeAt(0) )) return 1;
        if( isTangled( edge->getNodeAt(1) )) return 1;
        return 0;
    }
};

struct EdgeLengthWeight : public EdgeWeight {
    ~EdgeLengthWeight() {}
    string getName() const
    {
        return "Length";
    }
    bool isSymmetric() const
    {
        return 1;
    }
    double getWeight(const JEdgePtr &edge, int)
    {
        return JEdgeGeometry::getLength( edge );
    }
};

struct TutteWeight : public EdgeWeight {
    ~TutteWeight() {}
    string getName() const
    {
        return "Tutte";
    }
    bool isSymmetric() const
    {
        return 0;
    }
    double getWeight( const JEdgePtr &edge, int dir = 1);
};

struct CotangentWeight : public EdgeWeight {
    ~CotangentWeight() {}
    string getName() const
    {
        return "Cotangent";
    }
    bool   isSymmetric() const
    {
        return 1;
    }
    double getWeight( const JEdgePtr &edge, int dir = 1);
    double getCotangent(const Point3D &org, const Point3D &p1, const Point3D &p2) const;
};

struct FloaterMeanValueWeight : public EdgeWeight {
    ~FloaterMeanValueWeight() {}
    string getName() const
    {
        return "FloaterMeanValue";
    }

    bool  isSymmetric() const
    {
        return 0;
    }

    double getWeight( const JEdgePtr &, int dir = 1);
};

/*
struct LeastSquareWeight : public EdgeWeight {
    ~LeastSquareWeight() {}
    string getName() const { return "Normalized"; }

    bool isSymmetric() const {
        return 1;
    }
    double getWeight(const Edge *, int dir = 1){
        return 1.0;
    }
};
*/

struct NormalizedWeight : public EdgeWeight {
    ~NormalizedWeight() {}
    string getName() const
    {
        return "Normalized";
    }
    bool isSymmetric() const
    {
        return 1;
    }

    double getWeight(const JEdgePtr &edge, int dir = 1);
};

///////////////////////////////////////////////////////////////////////////////
class JMeshSmoother: public JMeshGeometricOptimization
{
};

class JLaplaceMeshSmoother: public JMeshSmoother
{
//    typedef Math::GeneralSparseMatrix<double> GMatrix;
//    typedef Math::SymmetricSparseMatrix<double> SMatrix;
    typedef Eigen::Triplet<double> Triplet_t;
public:

    static const int  COMBINATORIAL      = 0;
    static const int  COTANGENT          = 1;
    static const int  EDGE_LENGTH        = 3;
    static const int  FLOATER_MEAN_VALUE = 4;
    static const int  LEAST_SQUARE       = 5;
    static const int  NORMALIZED         = 6;
    static const int  TUTTE              = 7;

    static const int  EXPLICIT_METHOD    = 0;
    static const int  IMPLICIT_METHOD    = 1;

    JLaplaceMeshSmoother();

    ~JLaplaceMeshSmoother() {}

    bool isSymmetric() const;

    void setTaubinSteps( bool a)
    {
        applyTaubinSteps = a;
    }

    void setNumericalMethod(int t)
    {
        numericalMethod = t;
    }

    void setBandPassFreq( double f)
    {
        Kpb = f;
    }

    void setLambda(double val)
    {
        lambda = std::max(0.001, val);
    }

    void setMu(double m)
    {
        mu = m;
    }

    // Set the mesh. When we set the mesh, boundary nodes are included in
    // constraint nodes sequence.. We can add more constraints, after
    // calling this function...
    int  setMesh( const JMeshPtr &m);

    void setEdgesWeight(int w, bool update_weight_every_step = 0);

    void setNumIterations( int n )
    {
        numIterations = std::max(1,n);
    }

    void checkInversion( bool b) {
        check_element_inversion = b;
    }

    int getNumIterations() const {
        return numIterations;
    }

    void setTolerance( double t ) {
        tolerance = t;
    }

    double getMaxResidual() const
    {
        return maxDiff;
    }
    void setDualSmooth( bool b) {
        dualSmooth = b;
    }

    int smoothAll();
    int smooth( const JNodeSequence &nodes);

    void clearAll()
    {
        if( weightMatrix ) weightMatrix.reset();
        shrink_to_zero( newPos);
        shrink_to_zero( fixednodes);
    }

//  int untangle();

protected:
    JMeshPtr  mesh;
    int     topDim;
    int     numIterations = 1;    // Number of iterations for explicit schemes.
    int     numericalMethod  = EXPLICIT_METHOD;
    int     edgeWeightID;
    int     applyTaubinSteps = 0;
    bool    dualSmooth = 0;
    bool    update_weight_every_step = 0;
    bool    check_element_inversion  = 0;
    double  lambda, mu, Kpb;    // Taubin's method parameters, see the paper for details.
    double  maxDiff, tolerance;

    boost::scoped_ptr<Eigen::SparseMatrix<double> > weightMatrix;
    boost::scoped_ptr<EdgeWeight> edgeWeight;
    vector<Point3D> newPos;
    vector<double>  faceQuality, oldNodeQuality, newNodeQuality;

    JNodeSequence fixednodes;
    JNodeSequence neighs;

    void update(const JNodePtr &vtx);

    void assignEdgeWeight();
    void getQuality();

    int atomicOp(const JNodePtr &v, const JNodeSequence &vneigh, double coeff);
    int atomicOp(const JNodePtr &v, const JEdgeSequence &vneigh, double coeff);

    int smooth_with_primal_nodes();
    int smooth_with_dual_nodes();


    int getWeightMatrix();
    int useGPtoolBox();
};

