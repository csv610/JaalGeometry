#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "basic_math.hpp"

#include "slu_ddefs.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

int LinearSystem::solve( SparseMatrix<double> &K, vector<double> &b, vector<double> &x)
{
    int nRows = K.rows();
    int nCols = K.cols();
    int nonZeros = K.nonZeros();
    assert( nRows == nCols );
    assert( nonZeros > 0);

    SuperMatrix A, L, U, B;

    K.makeCompressed();
    double *val     = K.valuePtr();
    assert( val != nullptr);
    int    *col     = K.innerIndexPtr();
    assert( col != nullptr);
    int    *offset  = K.outerIndexPtr();
    assert( offset != nullptr);

    dCreate_CompCol_Matrix(&A, nRows, nCols, nonZeros, val, col, offset, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&B, nRows, 1, &b[0], nRows, SLU_DN, SLU_D, SLU_GE);

    vector<int> permRow, permCol;
    permRow.resize(nRows);
    permCol.resize(nCols);

    superlu_options_t options;
    set_default_options(&options);

    options.ColPerm = NATURAL;

    SuperLUStat_t stat;
    StatInit(&stat);

    int info;
    dgssv(&options, &A, &permCol[0], &permRow[0], &L, &U, &B, &stat, &info);

    DNformat *AStore = (DNformat *)B.Store;
    x.resize(nRows);
    double *xtmp = (double *)AStore->nzval;
    for( int i = 0; i < nRows; i++) x[i] = xtmp[i];

#ifdef DEBUG
    for( int i = 0; i < nRows; i++) {
        if( permRow[i] != i ) {
            cout << "Warning: There was row permutation in the matrix " << endl;
            break;
        }
    }
#endif

    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);

    return info;
}

