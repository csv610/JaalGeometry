#include "HeatConduction3D.hpp"


HeatConduction3D:: HeatConduction3D()
{
    mesh = nullptr;
    numIterations = 100;
    attribname = "Temperature";
}

////////////////////////////////////////////////////////////////////////////////

void HeatConduction3D:: get_local_tet_mat(Vec4D &x, Vec4D &y, Vec4D &z, Mat4D &lmat)
{
    Mat4D g;

    double x23 = x[1] - x[2];
    double y23 = y[1] - y[2];
    double z23 = z[1] - z[2];
    double x34 = x[2] - x[3];
    double y34 = y[2] - y[3];
    double z34 = z[2] - z[3];
    double x41 = x[3] - x[0];
    double y41 = y[3] - y[0];
    double z41 = z[3] - z[0];
    double x12 = x[0] - x[1];
    double y12 = y[0] - y[1];
    double z12 = z[0] - z[1];

    double xy234 = x34*y23 - x23*y34;
    double xy341 = x34*y41 - x41*y34;
    double xy412 = x12*y41 - x41*y12;
    double xy123 = x12*y23 - x23*y12;
    double xz234 = x23*z34 - x34*z23;
    double xz341 = x41*z34 - x34*z41;
    double xz412 = x41*z12 - x12*z41;
    double xz123 = x23*z12 - x12*z23;
    double yz234 = y34*z23 - y23*z34;
    double yz341 = y34*z41 - y41*z34;
    double yz412 = y12*z41 - y41*z12;
    double yz123 = y12*z23 - y23*z12;

    g(0,0) = -(x[1]*yz234 + y[1]*xz234 + z[1]*xy234);
    g(1,0) = -(x[2]*yz341 + y[2]*xz341 + z[2]*xy341);
    g(2,0) = -(x[3]*yz412 + y[3]*xz412 + z[3]*xy412);
    g(3,0) = -(x[0]*yz123 + y[0]*xz123 + z[0]*xy123);

    g(0,1) = yz234;
    g(1,1) = yz341;
    g(2,1) = yz412;
    g(3,1) = yz123;

    g(0,2) = xz234;
    g(1,2) = xz341;
    g(2,2) = xz412;
    g(3,2) = xz123;

    g(0,3) = xy234;
    g(1,3) = xy341;
    g(2,3) = xy412;
    g(3,3) = xy123;
    //
    //    vol36 = 1/(36*volume) = 1/(6*determinant)
    //
    double det, vol36;
    det = fabs(x[0]*yz234 + x[1]*yz341 + x[2]*yz412 + x[3]*yz123);

    vol36 = 1.0/(6.0*det);
    if( vol36 < 1.0E-15) {
        cout << "Fatal Error: Volume of tet is: " << vol36 << endl;
        exit(0);
    }

    for( int j = 0; j < 4; j++) {
        for( int i = 0; i < 4; i++) {
            lmat(i,j) = vol36*( g(i,1)*g(j,1) +
                                g(i,2)*g(j,2) +
                                g(i,3)*g(j,3) );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void HeatConduction3D::build_global_matrix( )
{
    mesh->enumerate(0);

    Mat4D lmat;
    Vec4D xt, yt, zt;

    size_t numRows = mesh->getSize(0);
    M.resize(numRows,numRows);

    size_t numCells = mesh->getSize(3);
    for( size_t icell = 0; icell < numCells; icell++) {
        JCellPtr tet  = mesh->getCellAt(icell);
        if( tet->isActive() ) {
            for( int inode = 0; inode < 4; inode++) {
                JNodePtr vertex = tet->getNodeAt(inode);
                const Point3D &xyz = vertex->getXYZCoords();
                xt[inode] = xyz[0];
                yt[inode] = xyz[1];
                zt[inode] = xyz[2];
            }
            get_local_tet_mat(xt, yt, zt, lmat);
            for( int i = 0; i < 4; i++) {
                int ii  = tet->getNodeAt(i)->getID();
                for( int j = 0; j < 4; j++) {
                    int jj  = tet->getNodeAt(j)->getID();
                    M.coeffRef(ii,jj) += 0.5*lmat(i,j);
                }
            }
        }
    }
    matrix_ready = 1;
}

//////////////////////////////////////////////////////////////////////

void HeatConduction3D :: apply_boundary_conditions()
{
    A = M;

    int numRows = mesh->getSize(0);

    rhs.resize( numRows );
    for( int i = 0; i < numRows; i++) rhs[i] = 0.0;

    JNodeSequence vneighs;

    numBoundConditions = 0;
    double coeff, val = 0.0;
    for( size_t k = 0; k < numRows; k++) {
        JNodePtr vertex = mesh->getNodeAt(k);
        if( vertex->hasAttribute("Dirichlet") ) {
            numBoundConditions++;
            vertex->getAttribute("Dirichlet", val);
            JNode::getRelations(vertex, vneighs);

            int irow =  vertex->getID();
            for( size_t j = 0; j < vneighs.size(); j++) {
                int icol =  vneighs[j]->getID();
                coeff = M.coeffRef(irow, icol);
                rhs[irow] -= coeff*val;
                A.coeffRef(irow,icol) = 0.0;

                int jrow = icol;
                int jcol = irow;
                coeff = M.coeffRef(jrow, jcol);
                rhs[jrow] -= coeff*val;
                A.coeffRef(jrow,jcol) = 0.0;
            }
            A.coeffRef(irow,irow) = 1.0;
        }
    }
}
//////////////////////////////////////////////////////////////////////

void HeatConduction3D::solve_linear_system()
{
    ConjugateGradient<SparseMatrix<double> > cg;
    cg.compute(A);
    cg.setMaxIterations(numIterations);
    temperature = cg.solve(rhs);
}

//////////////////////////////////////////////////////////////////////

int HeatConduction3D :: solve()
{
    if( mesh == nullptr ) return 1;
    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 3 ) return 2;

    int elemType = mesh->getTopology()->getElementsType(3);
    if( elemType != JCell::TETRAHEDRON) return 2;

    if( !matrix_ready )  build_global_matrix();

    apply_boundary_conditions();

    if( numBoundConditions == 0) return 1;

    solve_linear_system();

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() )
            vtx->setAttribute(attribname, temperature[vtx->getID()] );
    }

    return 0;
}
//////////////////////////////////////////////////////////////////////

