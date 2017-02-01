#include "ConjugateField.hpp"

using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////

void JConjugateField :: genRandomConstraints(int nRandom)
{
    int num = 2;

    b.resize(nRandom);

    int minVal = 0;
    int maxVal = mesh->getSize(2)-1;

    for( int i = 0; i < nRandom; i++)
        b <<  JMath::random_value(minVal, maxVal);

    bc.resize(b.size(), num*3);
    for (unsigned i=0; i<b.size(); ++i)
    {
        VectorXd t = genRandomVector(B1.row(b(i)), B2.row(b(i)),num);
        bc.row(i) = t;
    }
}

////////////////////////////////////////////////////////////////////////////////

int JConjugateField :: genField()
{
/*
        // Interpolate to get a smooth field
        igl::n_polyvector(V, F, b, bc, smooth_pvf);

        // Initialize conjugate field with smooth field
        csdata = new igl::ConjugateFFSolverData<Eigen::MatrixXd,Eigen::MatrixXi>(V,F);
        conjugate_pvf = smooth_pvf;

        // Optimize the field
        int    conjIter = 20;
        bool   doHardConstraints = true;
        double lambdaOrtho = .1;
        double lambdaInit = 100;
        double lambdaMultFactor = 1.01;
        double lambdaOut;

        VectorXi isConstrained = VectorXi::Constant(F.rows(),0);
        for (unsigned i=0; i<b.size(); ++i)
            isConstrained(b(i)) = 1;

        igl::conjugate_frame_fields(*csdata, isConstrained, conjugate_pvf, conjugate_pvf,
                                    conjIter, lambdaOrtho, lambdaInit, lambdaMultFactor,
                                    doHardConstraints, &lambdaOut);

        vecField = JMesh::newObject();
        addEdges(vecField, B - global_scale*conjugate_pvf.block(0,0,F.rows(),3),
                 B + global_scale*conjugate_pvf.block(0,0,F.rows(),3));
        addEdges(vecField, B - global_scale*conjugate_pvf.block(0,3,F.rows(),3),
                 B + global_scale*conjugate_pvf.block(0,3,F.rows(),3));

        delete csdata;
*/


    /*
        // local representations of field vectors
        Eigen::Matrix<double, Eigen::Dynamic, 2> pvU, pvV;
        pvU.resize(F.rows(),2);
        pvV.resize(F.rows(),2);
        //smooth
        const Eigen::MatrixXd &Us = smooth_pvf.leftCols(3);
        const Eigen::MatrixXd &Vs = smooth_pvf.rightCols(3);
        pvU << igl::dot_row(Us,B1), igl::dot_row(Us,B2);
        pvV << igl::dot_row(Vs,B1), igl::dot_row(Vs,B2);
        csdata->evaluateConjugacy(pvU, pvV, conjugacy_s);

        //conjugate
        const Eigen::MatrixXd &Uc = conjugate_pvf.leftCols(3);
        const Eigen::MatrixXd &Vc = conjugate_pvf.rightCols(3);
        pvU << igl::dot_row(Uc,B1), igl::dot_row(Uc,B2);
        pvV << igl::dot_row(Vc,B1), igl::dot_row(Vc,B2);
        csdata->evaluateConjugacy(pvU, pvV, conjugacy_c);

        // Highlight in red the constrained faces
        MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
          for (unsigned i=0; i<b.size();++i)
              C.row(b(i)) << 1, 0, 0;
        double maxC = std::max(conjugacy_c.maxCoeff(), conjugacy_s.maxCoeff());
        double minC = std::min(conjugacy_c.minCoeff(), conjugacy_s.minCoeff());

        Eigen::VectorXd valS = conjugacy_s;
        // Eigen::VectorXd valS = (valS.array() - minC)/(maxC-minC);
        // valS = 1 - valS.array();
        Eigen::VectorXd valC = conjugacy_c;
        // Eigen::VectorXd valC = (valC.array() - minC)/(maxC-minC);
        // valC = 1 - valC.array();
        MatrixXd CS, CC;
        igl::jet(valS, 0, 0.004, CS);
        igl::jet(valC, 0, 0.004, CC);

        if (key == 1)
        {
            // Frame field constraints
            MatrixXd F1_t = MatrixXd::Zero(F.rows(),3);
            MatrixXd F2_t = MatrixXd::Zero(F.rows(),3);

            for (unsigned i=0; i<b.size(); ++i)
            {
                F1_t.row(b(i)) = bc.block(i,0,1,3);
                F2_t.row(b(i)) = bc.block(i,3,1,3);
            }

            addEdges(B - global_scale*F1_t, B + global_scale*F1_t);
            addEdges(B - global_scale*F2_t, B + global_scale*F2_t);
        }


        if (key == 2)
        {
            // Interpolated result
            addEdges(B - global_scale*smooth_pvf.block(0,0,F.rows(),3),
                     B + global_scale*smooth_pvf.block(0,0,F.rows(),3)),
            addEdges(B - global_scale*smooth_pvf.block(0,3,F.rows(),3),
                     B + global_scale*smooth_pvf.block(0,3,F.rows(),3));
        }

        if (key == 3)
        {
            // Conjugate field
        }
    */
    return 0;
}

/////////////////////////////////////////////////////////////////////////

void JConjugateField :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

#ifdef USE_IGL

    JMeshEigenMatrix mat;
    mat.setMesh(mesh);
    V = mat.getNodeMatrix();
    F = mat.getFaceMatrix();

    igl::barycenter(V, F, B);

    // Local bases (needed for conjugacy)
    igl::local_basis(V, F, B1, B2, B3);

    // Compute scale for visualizing fields
    global_scale =  0.5*igl::avg_edge_length(V, F);
#endif
}
