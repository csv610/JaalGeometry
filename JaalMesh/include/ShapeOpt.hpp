#pragma once

#include "MeshOptimization.hpp"
#include "Constraint.hpp"
#include "ShapeOpSolver.hpp"

class JShapeOptimizer : public JMeshGeometricOptimization
{
public:
    void setMesh( const JMeshPtr &m);

    void addPlaneConstraint( const JNodeSequence &n, double w = 1.0);
    void addCloseConstraint( const JNodeSequence &n, double w = 1.0);
    void addLineConstraint( const JEdgeSequence &f, double w = 1.0);
    void addEdgeStrainConstraint( const JEdgeSequence &f, double len, double w = 1.0);
    void addAreaConstraint( const JFaceSequence &f, double w = 1.0);
    void addCircleConstraint( const JFaceSequence &f, double w = 1.0);
    void addPlaneConstraint( const JFaceSequence &f, double w = 1.0);
    void addCloseConstraint( const JFaceSequence &f, double w = 1.0);
    void addVolumeConstraint( const JCellSequence &f, double w = 1.0);
    void addSphereConstraint( const JCellSequence &f, double w = 1.0);

    void addAreaConstraints( double w = 1.0);
    void addEdgeStrainConstraints( double w = 1.0);
    void addTriangleStrainConstraints( double w = 1.0);
    void addCircleConstraints( double w = 1.0);
    void addPlaneConstraints( double w = 1.0);
    void addBoundaryConstraints( double w = 1.0);
    void addLaplaceConstraints( double w = 1.0);
    void addRectangleConstraints( double w = 1.0);
    void addParallelogramConstraints( double w = 1.0);
    void addSimilarityConstraints( double w = 1.0);
    void addRigidConstraintss( double w = 1.0);
    void addLengthConstraints( double w = 1.0);

    // Specific to 3D 
    void addVolumeConstraints( double w = 1.0);
    void addSphereConstraints( double w = 1.0);

    void removeAreaConstraints() { cdb.erase("Area"); }
    void removeEdgeStrainConstraints() { cdb.erase("EdgeStrain"); }
    void removeTriangleStrainConstraints() { cdb.erase("TriangleStrain"); }
    void removeCircleConstraints() { cdb.erase("Circle"); }
    void removePlaneConstraints() { cdb.erase("Plane"); }
    void removeLaplaceConstraints() { cdb.erase("Laplace"); }
    void removeRectangleConstraints() { cdb.erase("Rectangle"); }
    void removeParallelogramConstraints() { cdb.erase("Parallelogram"); }
    void removeSimilarityConstraints() { cdb.erase("Similarity"); }
    void removeRigidConstraints() { cdb.erase("Rigid"); }
    void removeVolumeConstraints() { cdb.erase("Volume"); }
    void removeSphereConstraints() { cdb.erase("Sphere"); }
    void removeCrossFieldConstraints() { cdb.erase("CrossField"); }

    vector<Point3D> getTarget(const string &s);

    void solve();
private:
    std::map<string, vector<std::shared_ptr<ShapeOp::Constraint>>> cdb;
    ShapeOp::Matrix3X pCloud;
    ShapeOp::Solver   solver;
};

