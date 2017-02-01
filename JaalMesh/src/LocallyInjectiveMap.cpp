// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "LocallyInjectiveMap.hpp"

using namespace std;
using namespace Eigen;

///////////////////////////////////////////////////////////////////////////////

JLocallyInjectiveMap::JLocallyInjectiveMap()
{
    initParams();
};

///////////////////////////////////////////////////////////////////////////////

JLocallyInjectiveMap::~JLocallyInjectiveMap()
{
};

///////////////////////////////////////////////////////////////////////////////
int JLocallyInjectiveMap :: getEnergyType( const string &str)
{
    if( str == "Dirichlet") return DIRICHLET_ENERGY;
    if( str == "Laplacian") return LAPLACIAN_ENERGY;
    if( str == "ARAP")      return AS_RIGID_AS_POSSIBLE_ENERGY;
    if( str == "Conformal") return LEAST_SQUARE_CONFORMAL_ENERGY;
    return -1;
}
///////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMap:: initParams()
{
    dim   = 0;
    energyType = LAPLACIAN_ENERGY;
    enableOutput = false;
    enableBarriers = true;
    enableAlphaUpdate = true;
    tolerance  = 1.0E-05;

    beta  = 0.01;
    gamma = 1.0;

    findLocalMinima = true;
    enableSubstepping = true;
    enableLogBarriers = false;
    enableNeoHookeanBarriers = false;
    enableBarrierCompensation = false;

    positionalConstraintError = 0;

    /*
        alpha = 1.0E+08;
        alphaRatio = 1.0E+03;
    */

    smallestArea = -1;       // Default: 1.0E-05*smallest triangle...
    maxIterations = 1000;
    /*
        limData = nullptr;
    */
}

///////////////////////////////////////////////////////////////////////////////
void JLocallyInjectiveMap:: setMesh( const JMeshPtr &m )
{
    inMesh = m;
    derivedMesh = 0;
    initMesh();

}
///////////////////////////////////////////////////////////////////////////////
void JLocallyInjectiveMap::createSimplicialMesh()
{
    simplicialMesh.reset();
    if( inMesh == nullptr) return;

    inMesh->pruneAll();
    inMesh->enumerate(0);
    dim = inMesh->getTopology()->getDimension();
    if( dim < 2 || dim > 3) {
        cout << "Warning: Locallly injective mapping only for 2D and 3D simplicial elements" << endl;
        return;
    }

    int eType =  inMesh->getTopology()->getElementsType(dim);

    int cid;
    Point3D xyz;
    simplicialMesh = inMesh;
    if( dim == 2 && eType == JFace::QUADRILATERAL) {
        AllTriMeshGenerator alltri;
        JMeshPtr tmpmesh = inMesh->deepCopy();
        simplicialMesh = alltri.getFromQuadMesh(tmpmesh, 4);
        int nSize = inMesh->getSize(0);
        for( size_t i = 0; i < nSize; i++) {
            const JNodePtr &v0 = inMesh->getNodeAt(i);
            int err = v0->getAttribute("TargetPos", xyz);
            if( !err ) {
                const JNodePtr &v1 = simplicialMesh->getNodeAt(i);
                v1->setAttribute("TargetPos", xyz);
            }
            err = v0->getAttribute("Constraint", cid);
            if( !err) {
                const JNodePtr &v1 = simplicialMesh->getNodeAt(i);
                v1->setAttribute("Constraint", cid);
            }

        }

        derivedMesh = 1;
        return;
    }

    /*
            if( dim == 3 && elemType != JCell::TETRAHEDRON) {
                cout << "Warning: Locallly injective mapping only for 3D tetrahedra elements" << endl;
                return;
            }
    */
}
///////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMap::initMesh()
{
    if( inMesh == nullptr ) return;

    createSimplicialMesh();

    if( simplicialMesh == nullptr) return;

    elemType =  simplicialMesh->getTopology()->getElementsType(dim);

    int numNodes = simplicialMesh->getSize(0);

    initialVertices.resize(numNodes,3);
    deformedVertices.resize(numNodes,3);

//  #pragma omp parallel for
    for(int i=0; i< numNodes; i++) {
        const JNodePtr &vtx = simplicialMesh->getNodeAt(i);
        const Point3D  &xyz = vtx->getXYZCoords();
        for(int j=0; j< 3; j++) {
            initialVertices.coeffRef(i,j)  = xyz[j];
            deformedVertices.coeffRef(i,j) = xyz[j];
        }
    }

    int nSize  = numNodes*dim;
    constraintMatrix.resize(nSize, nSize);
    constraintTargets.resize(nSize);
    constraintTargets.setZero();

    if( elemType == JFace::TRIANGLE ) {
        int numTriangles = simplicialMesh->getSize(2);
        elements.resize(numTriangles,3);
//      #pragma omp parallel for
        for(int i=0; i < numTriangles; i++) {
            const JFacePtr &face = simplicialMesh->getFaceAt(i);
            for(int j= 0; j< 3; j++)
                elements.coeffRef(i,j) = face->getNodeAt(j)->getID();
        }
    }

    if( elemType == JCell::TETRAHEDRON ) {
        int numCells = simplicialMesh->getSize(3);
        elements.resize(numCells,4);
//        #pragma omp parallel for
        for(int i= 0; i < numCells; i++) {
            const JCellPtr &cell = simplicialMesh->getCellAt(i);
            for(int j= 0; j< 4; j++)
                elements.coeffRef(i,j) = cell->getNodeAt(j)->getID();
        }
    }

    vector<Eigen::Triplet<double> > triplets;
    for(int i=0; i< nSize; i++)
        triplets.push_back(Triplet<double>(i,i,1));
    constraintMatrix.setFromTriplets(triplets.begin(),triplets.end());
}
///////////////////////////////////////////////////////////////////////////////
double JLocallyInjectiveMap::getMaxDistance() const
{
    if( simplicialMesh == nullptr ) return 0.0;

    Point3D pdst, pcurr;

    double maxDiff = 0;
    int numNodes   =  simplicialMesh->getSize(0);

    #pragma omp parallel for
    for( int i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = simplicialMesh->getNodeAt(i);
        int err = vtx->getAttribute("TargetPos", pdst);
        if( !err) {
            pcurr     = vtx->getXYZCoords();
            double dx =  pdst[0] - pcurr[0];
            double dy =  pdst[1] - pcurr[1];
            double dz =  pdst[2] - pcurr[2];
            maxDiff   =  std::max(maxDiff, dx*dx + dy*dy + dz*dz);
        }
    }
    return maxDiff;
}

///////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMap::setTargets()
{
    if( simplicialMesh == nullptr ) return;

    int numNodes = simplicialMesh->getSize(0);

    for(int i=0; i< numNodes*dim; i++)
        constraintMatrix.coeffRef(i,i) = 0;

    constraintTargets.setZero();

    Point3D xyz;
//  #pragma omp parallel for  private(xyz)
    for( int i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = simplicialMesh->getNodeAt(i);
        if( vtx->isActive() ) {
            if( vtx->hasAttribute("TargetPos") ) {
                int vid = vtx->getID();
                vtx->getAttribute("TargetPos", xyz);
                for(int j=0; j< dim; j++) {
                    constraintMatrix.coeffRef(vid*dim+j,vid*dim+j) = 1;
                    constraintTargets.coeffRef(vid*dim+j) = xyz[j];
                }
            }
            if( vtx->hasAttribute("Constraint")) {
                int vid = vtx->getID();
                xyz = vtx->getXYZCoords();
                for(int j=0; j< dim; j++) {
                    constraintMatrix.coeffRef(vid*dim+j,vid*dim+j) = 1;
                    constraintTargets.coeffRef(vid*dim+j) = xyz[j];
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMap::setEnergyType(int energy)
{
    energyType = energy;
}
///////////////////////////////////////////////////////////////////////////////

int JLocallyInjectiveMap::getEnergy() const
{
    return energyType;
}
///////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMap::updateMesh()
{
    if( simplicialMesh == nullptr) return;

    Point3D xyz;

    size_t numNodes  = simplicialMesh->getSize(0);
    #pragma omp parallel for  private(xyz)
    for(size_t i=0; i< numNodes; i++) {
        for(int j=0; j< 3; j++)
            xyz[j] = deformedVertices.coeff(i,j);
        const JNodePtr &vtx = simplicialMesh->getNodeAt(i);
        vtx->setXYZCoords(xyz);
    }

    if( derivedMesh ) {
        numNodes  = inMesh->getSize(0);
        #pragma omp parallel for  private(xyz)
        for(size_t i=0; i< numNodes; i++) {
            for(int j=0; j< 3; j++)
                xyz[j] = deformedVertices.coeff(i,j);
            const JNodePtr &vtx = inMesh->getNodeAt(i);
            vtx->setXYZCoords(xyz);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void JLocallyInjectiveMap::setOriginalMesh() const
{
    if( simplicialMesh == nullptr) return;

    Point3D xyz;
    size_t numNodes  = simplicialMesh->getSize(0);
    #pragma omp parallel for  private(xyz)
    for(int i=0; i< numNodes; i++) {
        for(int j=0; j< 3; j++)
            xyz[j] = initialVertices.coeff(i,j);
        const JNodePtr &vtx = simplicialMesh->getNodeAt(i);
        vtx->setXYZCoords(xyz);
    }

    if( derivedMesh) {
        numNodes  = inMesh->getSize(0);
        #pragma omp parallel for  private(xyz)
        for(size_t i=0; i< numNodes; i++) {
            for(int j=0; j< 3; j++)
                xyz[j] = initialVertices.coeff(i,j);
            const JNodePtr &vtx = inMesh->getNodeAt(i);
            vtx->setXYZCoords(xyz);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMap::setDeformedMesh() const
{
    if( simplicialMesh == nullptr ) return;
    size_t numNodes  = simplicialMesh->getSize(0);

    Point3D xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    #pragma omp parallel for  private(xyz)
    for(size_t i=0; i< numNodes; i++) {
        const JNodePtr &vtx = simplicialMesh->getNodeAt(i);
        for(int j=0; j< 3; j++)
            xyz[j] = deformedVertices.coeff(i,j);
        vtx->setXYZCoords(xyz);
    }
}

///////////////////////////////////////////////////////////////////////////////////
int JLocallyInjectiveMap::initSolve()
{
    setTargets();
    /*
        limData = InitLIM( deformedVertices, initialVertices, elements,
                           borderVertices, gradients, constraintMatrix,
                           constraintTargets, energyType, enableOutput,
                           enableBarriers, enableAlphaUpdate);
    */
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////

int JLocallyInjectiveMap::stepSolve()
{
    /*
        if( limData->solver->CurrentStepSize < 1e-15 ||
          ( limData->solver->CurrentPositionalEnergy <= tolerance &&
          ( findLocalMinima == false || limData->solver->CurrentStepSize < 1e-15)))
            return 1; // termination criteria fulfilled

    //  setTargets();

    int stat;
    limData->iteration = 0;
    for( int i = 0; i < maxIterations; i++) {
        if(limData->iteration >= maxIterations) return -1;
        stat = ComputeLIM_Step(limData, deformedVertices);
    }

    if( stat == -2) updateMesh();

    return stat;
    */
}

///////////////////////////////////////////////////////////////////////////////////////

int JLocallyInjectiveMap::solve()
{
    setTargets();

    EnergyType entype = static_cast<EnergyType>(energyType);

    int stat = ComputeLIM( deformedVertices, initialVertices, elements,
                           borderVertices, gradients, constraintMatrix,
                           constraintTargets, entype, tolerance,
                           maxIterations, findLocalMinima, enableOutput,
                           enableAlphaUpdate);
    if( stat != -2) updateMesh();
    return stat;


    /*
      if( limData == nullptr)  {
      limData = InitLIM(deformedVertices, initialVertices, elements, borderVertices, gradients, constraintMatrix, constraintTargets, energyType, enableOutput, enableBarriers, enableAlphaUpdate);

      }

      limData->iteration = 0;

      int result = 0;
      while(result == 0)
      {
        if( limData->solver->CurrentStepSize < 1e-15 || (limData->solver->CurrentPositionalEnergy <= tolerance &&
            (findLocalMinima == false || limData->solver->CurrentStepSize < 1e-15)))
          result = 1; // termination criteria fulfilled

        if(limData->iteration >= maxIterations)
          result = -1; // max iteration reached

        if(result == 0)
        {
          if(limData->solver->Solve() == -1)
            result = -2; // state not feasible -> inverted elements
          else
          {
            // swap vertex buffers
            Eigen::Matrix<double,Eigen::Dynamic,3>* temp = limData->mesh->DeformedVertices;
            limData->mesh->DeformedVertices = limData->mesh->PredictedVertices;
            limData->mesh->PredictedVertices = temp;

            limData->iteration++;
          }
        }
      }

       if( result != -2) updateMesh();

      // assign resulting vertices
       deformedVertices = *limData->mesh->DeformedVertices;
       return result;
    */

}

////////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMap:: setConstraintsPosition(const vector<int> &constraintVertices, const Matrix<double,Dynamic,3>& positions)
{
    if( limData == nullptr) return;

    for(int i=0; i<constraintVertices.size(); i++)
    {
        int idx = constraintVertices[i];
        for(int j=0; j< dim; j++)
        {
            limData->mesh->ConstraintTargets->coeffRef(idx*dim+j) = positions.coeff(i,j);
        }
    }
    limData->solver->Restart();
}

////////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMap ::setConstraints(const vector<int>& constraintVertices)
{
    if( limData == nullptr) return;

    // free all constraint vertices
    int nRows = limData->mesh->InitalVertices->rows()*dim;

    for(int i=0; i< limData->mesh->InitalVertices->rows()*dim; i++)
        limData->mesh->ConstraintMatrix->coeffRef(i,i) = 0;
    limData->mesh->ConstraintTargets->setZero();

    // set new constraint vertices
    for(int i=0; i<constraintVertices.size(); i++)
    {
        int idx = constraintVertices[i];
        for(int c=0; c< dim; c++)
        {
            limData->mesh->ConstraintMatrix->coeffRef(idx*dim+c,idx*dim+c) = 1;
            limData->mesh->ConstraintTargets->coeffRef(idx*dim+c) = limData->mesh->DeformedVertices->coeff(idx,c);
        }
    }
    limData->solver->UpdatePositionalConstraintMatrix();
    limData->solver->Restart();
}
////////////////////////////////////////////////////////////////////////////////


