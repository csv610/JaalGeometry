#include "KillingVectorsDeformation.hpp"

#ifdef CSV


#include "logspiral.h"
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
#include <complex>
#include <limits>

#define SQRT_2 1.41421356
#define POINT_SIZE_SCALE 2

/*
#define VF_SCALE 1
#define LINE_WIDTH 2
#undef min
#undef max
*/

using namespace std;

#define ALLOC_MEMORY(p,type,size) \
p = (type *) malloc ((((size) <= 0) ? 1 : (size)) * sizeof (type)) ; \
if (p == (type *) NULL) \
{ \
    exit(1) ; \
}

#define REALLOC_MEMORY(p,type,size) \
p = (type *) realloc (p,(((size) <= 0) ? 1 : (size)) * sizeof (type)) ; \
if (p == (type *) NULL) \
{ \
    exit(1) ; \
}

#define FREE_MEMORY(p,type) \
if (p != (type *) NULL) \
{ \
    free (p) ; \
    p = (type *) NULL ; \
}

///////////////////////////////////////////////////////////////////////////////

void error_handler(int status, const char *file, int line, const char *message)
{
    /*
        qWarning("CHOLMOD error status %d", status);
        qWarning("File: %s", file);
        qWarning("Line: %d", line);
        qWarning("Message: %s", message);
    */
}

///////////////////////////////////////////////////////////////////////////////

template<class T>
KillingVectorsDeformation<T>::KillingVectorsDeformation(const string &filename)
{
    if (filename.rfind("obj") != string::npos) {
        modelType = 2;
        numVertices = 0;
        numFaces = 0;
        ifstream infile(filename.c_str(), ios::in);
        bool is_vt = false;
        T x,y,z;
        string a,b,c;

        while (!infile.eof()) {
            // get line
            string curLine;
            getline(infile, curLine);

            // read type of the line
            istringstream issLine(curLine);
            string linetype;
            issLine >> linetype;

            if (linetype == "v") {
                numVertices++;
                issLine >> x >> y >> z;
                Point2D<T> p(x,y);
                vertices.push_back(p);
                continue;
            }
            if (linetype == "vt") {
                is_vt = true;
                issLine >> x >> y;
                Point2D<T> p(x,y);
                texCoords.push_back(p);
                continue;
            }
            if (linetype == "f") {
                numFaces++;
                issLine >> a >> b >> c;
                char* a_dup = strdup(a.c_str());
                char* b_dup = strdup(b.c_str());
                char* c_dup = strdup(c.c_str());
                char* aa = strtok(a_dup,"/");
                char* bb = strtok(b_dup,"/");
                char* cc = strtok(c_dup,"/");
                Face fa;
                fa.f[0] = atoi(aa)-1;
                fa.f[1] = atoi(bb)-1;
                fa.f[2] = atoi(cc)-1;
                faces.push_back(fa);
                continue;
            }
            if (linetype == "#") continue;
//           if (linetype == "mtllib") issLine >> mtlFile; //supported mtl: stores first line indicating needed mtl for future saved obj models
        }
        if (is_vt == false) {
            texCoords.resize(numVertices);
            for (int i = 0; i < numVertices; i++)
                texCoords[i] = vertices[i];
        }
    }

    if (filename.rfind("off") != string::npos) {
        modelType = 1;
        ifstream infile(filename.c_str(), ios::in);
        string temp;
        infile >> temp;

        infile >> numVertices >> numFaces >> temp;

        vertices.resize(numVertices);
        T z;
        for (int i = 0; i < numVertices; i++)
            infile >> vertices[i].x >> vertices[i].y >> z;

        texCoords.resize(numVertices);
        for (int i = 0; i < numVertices; i++)
            texCoords[i] = vertices[i];

        int three;
        faces.resize(numFaces);
        for (int i = 0; i < numFaces; i++)
            infile >> three >> faces[i][0] >> faces[i][1] >> faces[i][2];
    }

    //setNumExtraEigenvectors(1);
    initialize();
}
////////////////////////////////////////////////////////////////////////////////

template<class T>
KillingVectorsDeformation<T>::KillingVectorsDeformation(KillingVectorsDeformation<T> &m)
{
    numVertices = m.numVertices;
    numFaces = m.numFaces;
    vertices = m.vertices;
    texCoords = m.texCoords;
    faces = m.faces;
    wireframeTrans = m.wireframeTrans;
    drawVFMode = m.drawVFMode;
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
KillingVectorsDeformation<T>::~KillingVectorsDeformation()
{
    FREE_MEMORY(Ax,double);
    FREE_MEMORY(Parent, LDL_int);
    FREE_MEMORY(Flag, LDL_int);
    FREE_MEMORY(Lp, LDL_int);
    FREE_MEMORY(Lnz, LDL_int);
    FREE_MEMORY(D, double);
    FREE_MEMORY(Pattern, LDL_int);
    FREE_MEMORY(Y, double);
    FREE_MEMORY(X, double);
    FREE_MEMORY(Pfw, LDL_int);
    FREE_MEMORY(Pinv, LDL_int);
    FREE_MEMORY(Li, LDL_int);
    FREE_MEMORY(Lx, double);
    cholmod_free_factor(&L2, cm);
    cholmod_finish(cm);
}
////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::initialize()
{
    transCount = 0;
    pCount = 0;

    minX = numeric_limits<T>::max();
    maxX = numeric_limits<T>::min();
    minY = numeric_limits<T>::max();
    maxY = numeric_limits<T>::min();

    T sumX = 0.0;
    T sumY = 0.0;

    for (int i = 0; i < numVertices; i++) {
        minX = min(minX, vertices[i].x);
        minY = min(minY, vertices[i].y);
        maxX = max(maxX, vertices[i].x);
        maxY = max(maxY, vertices[i].y);
        sumX = sumX + vertices[i].x;
        sumY = sumY + vertices[i].y;
    }

    T avgX = sumX/numVertices;
    T avgY = sumY/numVertices;

    for (int i=0; i<numVertices; i++) {
        vertices[i].x = vertices[i].x - avgX;
        vertices[i].y = vertices[i].y - avgY;
    }

    undoIndex = 0; //save for undo
    undoVertices[0] = vertices; //save for undo

    updateCovariance(); // needed to know how many nonzeros

    int nz = covariance.getNumNonzero();
    int n = 2*numVertices;
    ALLOC_MEMORY(Ax,double,nz);
    ALLOC_MEMORY(Parent, LDL_int, n);
    ALLOC_MEMORY(Flag, LDL_int, n);
    ALLOC_MEMORY(Lp, LDL_int, n+1);
    ALLOC_MEMORY(Lnz, LDL_int, n);
    ALLOC_MEMORY(D, double, n);
    ALLOC_MEMORY(Pattern, LDL_int, n) ;
    ALLOC_MEMORY(Y, double, n);
    ALLOC_MEMORY(X, double, n);
    ALLOC_MEMORY(Pfw, LDL_int, n);
    ALLOC_MEMORY(Pinv, LDL_int, n);
    Li = NULL;
    Lx = NULL;

    Ai = covariance.getAi();
    Ap = covariance.getAp();

    cout << "Computing permutation..." << endl;
    clock_t t = clock();
    double Info[AMD_INFO];
    if (amd_order (n, Ap, Ai, Pfw, (double *) NULL, Info) < AMD_OK) {
        cout << "call to AMD failed " << endl;
        exit(1);
    }
    amd_control((double*)NULL);
    amd_info(Info);

    cout << "AMD time: " << (clock()-t)/(double)CLOCKS_PER_SEC << endl;
    cout << "Doing symbolic ldl..." << endl;

    t = clock();
    ldl_symbolic (2*numVertices, Ap, Ai, Lp, Parent, Lnz, Flag, Pfw, Pinv);
    cout << "Symbolic time: " << (clock()-t)/(double)CLOCKS_PER_SEC << endl;
    //recomputeEigenvectors();


    neighbors.resize(numVertices);
    map< int , map<int,int> > edgeCount;
    for (int i = 0; i < numFaces; i++) {
        int a = faces[i][0], b = faces[i][1], c = faces[i][2];
        neighbors[a].insert(b);
        neighbors[a].insert(c);
        neighbors[b].insert(a);
        neighbors[b].insert(c);
        neighbors[c].insert(a);
        neighbors[c].insert(b);
        edgeCount[a][b]++;
        edgeCount[b][a]++;
        edgeCount[a][c]++;
        edgeCount[c][a]++;
        edgeCount[c][b]++;
        edgeCount[b][c]++;
    }

    for (int i = 0; i < numFaces; i++) {
        int a = faces[i][0], b = faces[i][1], c = faces[i][2];
        if (edgeCount[a][b] == 1) {
            boundaryVertices.insert(a);
            boundaryVertices.insert(b);
        }
        if (edgeCount[b][c] == 1) {
            boundaryVertices.insert(b);
            boundaryVertices.insert(c);
        }
        if (edgeCount[a][c] == 1) {
            boundaryVertices.insert(a);
            boundaryVertices.insert(c);
        }
    }

    cm = &Common;
    cholmod_start(cm);
    cm->error_handler = error_handler;

    /////////////////////////////////////////////////////////
    // Prefactor

    cholmod_dense *boundaryRHS = cholmod_zeros(2*numVertices, 1, CHOLMOD_REAL, cm);
    double *boundaryX = (double*)boundaryRHS->x;

    for (int i = 0; i < 2*numVertices; i++)
        boundaryX[i] = 0;

    // for fun constrain boundary to (1,1)
    for (set<int>::iterator it = boundaryVertices.begin(); it != boundaryVertices.end(); ++it) {
        boundaryX[*it] = 1;
        boundaryX[*it + numVertices] = 1;
    }

    pCount = 0;
    getP(Pcopy);
    pCount = 0;

    double *rhsMove = (double*)malloc(Pcopy.numRows()*sizeof(double));

    Pcopy.multiply(boundaryX, rhsMove);

    for (int i = 0; i < Pcopy.numRows(); i++)
        rhsMove[i] *= -1;

    Pcopy.zeroOutColumns(boundaryVertices);
    Pcopy.zeroOutColumns(boundaryVertices, numVertices);

    cholmod_dense *B2 = cholmod_zeros(Pcopy.numCols(), 1, CHOLMOD_REAL, cm);
    double *Bx = (double*)B2->x;
    for (int i = 0; i < Pcopy.numCols(); i++) Bx[i] = 0;

    Pcopy.transposeMultiply(rhsMove,Bx);

    vector<int> constrained;
    for (set<int>::iterator it = boundaryVertices.begin(); it != boundaryVertices.end(); ++it) {
        int bv = *it;
        constrained.push_back(*it);
        constrained.push_back(*it+numVertices);
        Bx[bv] += 1;
        Bx[bv+numVertices] += 1;
    }

    Pcopy.addConstraint(constrained,1);

    nz = Pcopy.getNumNonzero();
    T *AxTemp = Pcopy.getAx();
    REALLOC_MEMORY(Ax, double, nz);
    for (int i = 0; i < nz; i++) Ax[i] = (double)AxTemp[i];
    Ai = Pcopy.getAi();
    Ap = Pcopy.getAp();

    cSparse.nrow = Pcopy.numCols(); // P is row-major, so cholmod thinks it's P*
    cSparse.ncol = Pcopy.numRows();
    cSparse.nzmax = Pcopy.getNumNonzero();
    cSparse.p = Ap;
    cSparse.i = Ai;
    cSparse.x = Ax;
    cSparse.stype = 0;
    cSparse.itype = CHOLMOD_INT;
    cSparse.xtype = CHOLMOD_REAL;
    cSparse.sorted = 1;
    cSparse.packed = 1;
    A = &cSparse;

    L2 = cholmod_analyze(A, cm);

    cholmod_free_dense(&B2, cm);
    cholmod_free_dense(&boundaryRHS, cm);
    free(rhsMove);
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::displaceMesh(vector<int> &indices, vector< Vector2D<T> > &displacements, T alpha)
{
    if (indices.size() == 1) { // when only one vertex is constrained, move parallel
        for (int i = 0; i < numVertices; i++) {
            vertices[i][0] += displacements[0].x;
            vertices[i][1] += displacements[0].y;
        }
        return;
    }

    time_t time1 = clock();
    getP(P);
    Pcopy.copy(P);
    double constructTime = (clock() - time1)/(double)CLOCKS_PER_SEC;
    cout << "Construct P time: " << constructTime << endl;

    time_t total_time = clock();
    alpha = alpha / (2*indices.size()) * P.infinityNorm();

    rhs.resize(2*numVertices);

    for (int i = 0; i < 2*numVertices; i++) rhs[i] = 0;

    for (unsigned int i = 0; i < indices.size(); i++) {
        rhs[ indices[i] ] = displacements[i][0]*alpha*alpha;
        rhs[ indices[i]+numVertices ] = displacements[i][1]*alpha*alpha;
    }

    vector<int> indices2;
    for (unsigned int i = 0; i < indices.size(); i++) {
        indices2.push_back(indices[i]);
        indices2.push_back(indices[i] + numVertices);
    }

    P.addConstraint(indices2, alpha);

    int nz = P.getNumNonzero();
    T *AxTemp = P.getAx();
    REALLOC_MEMORY(Ax, double, nz);
    for (int i = 0; i < nz; i++) Ax[i] = (double)AxTemp[i];
    Ai = P.getAi();
    Ap = P.getAp();

    cholmod_sparse cSparse;
    cSparse.nrow = P.numCols(); // P is row-major, so cholmod thinks it's P*
    cSparse.ncol = P.numRows();
    cSparse.nzmax = P.getNumNonzero();
    cSparse.p = Ap;
    cSparse.i = Ai;
    cSparse.x = Ax;
    cSparse.stype = 0;
    cSparse.itype = CHOLMOD_INT;
    cSparse.xtype = CHOLMOD_REAL;
    cSparse.sorted = 1;
    cSparse.packed = 1;
    cholmod_sparse *A = &cSparse;

    // Try cholmod!
    cholmod_dense *Xcholmod, *B;
    int n = cSparse.nrow;

    B = cholmod_zeros(n, 1, cSparse.xtype, cm);
    double* Bx = (double*)B->x;

    for (int i = 0; i < n; i++)
        Bx[i] = rhs[i];

    cholmod_factor *L = cholmod_analyze(A, cm);

    cholmod_factorize(A, L, cm);

    Xcholmod = cholmod_solve(CHOLMOD_A, L, B, cm);

    double* Xx = (double*)Xcholmod->x;
    if (drawVFMode) {
        vfOrig.resize(numVertices);
        for (int i = 0; i < numVertices; i++) vfOrig[i] = Vector2D<double>(Xx[i],Xx[i+numVertices]);
    }

    // NOW, DIRICHLET SOLVE ////////////////////////////////////////////////////

    clock_t dirichletTime = clock();

    cholmod_dense *boundaryRHS = cholmod_zeros(2*numVertices, 1, cSparse.xtype, cm);
    double *boundaryX = (double*)boundaryRHS->x;

    for (int i = 0; i < 2*numVertices; i++)
        boundaryX[i] = 0;

    for (set<int>::iterator it = boundaryVertices.begin(); it != boundaryVertices.end(); ++it) {
        boundaryX[*it] = Xx[*it];
        boundaryX[*it + numVertices] = Xx[*it + numVertices];
    }

    double *rhsMove = (double*)malloc(Pcopy.numRows()*sizeof(double));
    Pcopy.multiply(boundaryX, rhsMove);

    for (int i = 0; i < Pcopy.numRows(); i++)
        rhsMove[i] *= -1;

    Pcopy.zeroOutColumns(boundaryVertices);
    Pcopy.zeroOutColumns(boundaryVertices, numVertices);

    cholmod_dense *B2 = cholmod_zeros(Pcopy.numCols(), 1, cSparse.xtype, cm);
    Bx = (double*)B2->x;
    for (int i = 0; i < Pcopy.numCols(); i++) Bx[i] = 0;

    Pcopy.transposeMultiply(rhsMove,Bx);

    vector<int> constrained;
    for (set<int>::iterator it = boundaryVertices.begin(); it != boundaryVertices.end(); ++it) {
        int bv = *it;
        constrained.push_back(*it);
        constrained.push_back(*it+numVertices);
        Bx[bv] += Xx[bv];
        Bx[bv+numVertices] += Xx[bv+numVertices];
    }

    Pcopy.addConstraint(constrained,1);

    nz = Pcopy.getNumNonzero();
    AxTemp = Pcopy.getAx();
    REALLOC_MEMORY(Ax, double, nz);
    for (int i = 0; i < nz; i++) Ax[i] = (double)AxTemp[i];
    Ai = Pcopy.getAi();
    Ap = Pcopy.getAp();

    cSparse.nrow = Pcopy.numCols(); // P is row-major, so cholmod thinks it's P*
    cSparse.ncol = Pcopy.numRows();
    cSparse.nzmax = Pcopy.getNumNonzero();
    cSparse.p = Ap;
    cSparse.i = Ai;
    cSparse.x = Ax;
    cSparse.stype = 0;
    cSparse.itype = CHOLMOD_INT;
    cSparse.xtype = CHOLMOD_REAL;
    cSparse.sorted = 1;
    cSparse.packed = 1;
    A = &cSparse;

    //cholmod_factor *L2 = cholmod_analyze(A, cm);
    cholmod_factorize(A, L2, cm);
    cholmod_dense *Xcholmod2 = cholmod_solve(CHOLMOD_A, L2, B2, cm);
    Xx = (double*)Xcholmod2->x; // UNCOMMENT TO RESTORE NON-BAD SOLVE

    cout << "Dirichlet time: " <<  (clock()-dirichletTime)/(double)CLOCKS_PER_SEC << endl;

    ////////////////////////////////////////////////////////////////////////////

    // Try complex number trick /////

    // This log spiral implementation is slow and double-counts edges -- eventually we can make it
    // faster and parallelize

    newPoints.resize(numVertices);
    counts.resize(numVertices);

    for (int i = 0; i < numVertices; i++) {
        counts[i] = 0;
        newPoints[i] = Point2D<double>(0,0);
    }

    for (int i = 0; i < numFaces; i++)
        for (int j = 0; j < 3; j++) {
            int e1 = faces[i][j];
            int e2 = faces[i][(j+1)%3];
            int vtx = faces[i][(j+2)%3];

            complex<double> v1(Xx[e1], Xx[e1+numVertices]);
            complex<double> v2(Xx[e2], Xx[e2+numVertices]);
            complex<double> p1(vertices[e1][0], vertices[e1][1]);
            complex<double> p2(vertices[e2][0], vertices[e2][1]);

            complex<double> z = (v1-v2)/(p1-p2);
            complex<double> p0 = (p2*v1-p1*v2)/(v1-v2);

            double c = z.real();
            double alpha = z.imag();
            Point2D<double> p(p0.real(),p0.imag());
            Point2D<double> l1(vertices[e1][0], vertices[e1][1]);
            Point2D<double> l2(vertices[e2][0], vertices[e2][1]);

            LogSpiral<double> spiral;
            spiral.p0 = p;
            spiral.c = c;
            spiral.alpha = alpha;

            Point2D<double> result1 = spiral.evaluate(l1,1);//logSpiral(p,c,alpha,l1,1);
            Point2D<double> result2 = spiral.evaluate(l2,1);//logSpiral(p,c,alpha,l2,1);

            // compute cotangent weights
            Vector2D<T> d1 = vertices[e1] - vertices[vtx];
            Vector2D<T> d2 = vertices[e2] - vertices[vtx];
            double angle = fabs(rotationAngle(d1,d2));
            double cotangent = 1;// / tan(angle);

            counts[e1] += cotangent;
            counts[e2] += cotangent;

            newPoints[e1] += result1*cotangent;
            newPoints[e2] += result2*cotangent;
        }

    /////////////////////////////////

    vf.resize(numVertices);

    if (drawVFMode) {
        for (int i = 0; i < numVertices; i++) vf[i] = Vector2D<double>(Xx[i],Xx[i+numVertices]);
    } else {
        for (int i = 0; i < numVertices; i++) {
            //vertices[i][0] += Xx[i];
            //vertices[i][1] += Xx[i+numVertices];
            Point2D<double> answer = newPoints[i]/counts[i];
            vertices[i] = Point2D<T>(answer.x,answer.y);
        }
    }

    cholmod_free_factor(&L, cm);
    //cholmod_free_factor(&L2, cm);
    cholmod_free_dense(&Xcholmod, cm);
    cholmod_free_dense(&Xcholmod2, cm);
    cholmod_free_dense(&B, cm);
    cholmod_free_dense(&B2, cm);
    cholmod_free_dense(&boundaryRHS, cm);

    free(rhsMove);
    double totalSolveTime = (clock() - total_time)/(double)CLOCKS_PER_SEC;
    cout << "Total solve time: " <<  totalSolveTime << endl;

    double fullTime = constructTime + totalSolveTime;
    double fps = 1./fullTime;
    cout << "Total fps: " << fps;
}
////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::updateCovariance()
{
    //static SimpleSparseMatrix<T> P;
    getP(P);

    //static SimpleSparseMatrix<T> trans;

    // do a second time for timing with memory allocated
    clock_t t = clock();
    P.transpose(trans);
    cout << "Transpose time: " << (clock()-t)/(double)CLOCKS_PER_SEC << endl;

    t = clock();
    if (transCount == 0)
        trans.multiply(P, covariance);
    else
        trans.parallelMultiply(P, covariance);

    transCount++;
    cout << "Covariance product time: " <<  (clock()-t)/(double)CLOCKS_PER_SEC << endl;
}
////////////////////////////////////////////////////////////////////////////////

template<class T>
bool checkBad(T *array, int n)
{
    for (int i = 0; i < n; i++)
        if (array[i] == array[i]-1) {
            return true;
        }
    return false;
}
////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::recomputeEigenvectors()
{
    updateCovariance();

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2*numVertices; j++)
            eigenvectors[i][j] = 0;

    time_t t = clock();

    T invNorm = 1. / sqrt(numVertices);
    T dots[2] = {0,0};
    for (int i = 0; i < numVertices; i++) {
        eigenvectors[0][i] = invNorm;
        eigenvectors[1][i+numVertices] = invNorm;
        eigenvectors[2][i] = -vertices[i][1];
        eigenvectors[2][i+numVertices] = vertices[i][0];
        dots[0] += invNorm * -vertices[i][1];
        dots[1] += invNorm * vertices[i][0];
    }

    eigenvectors[2] -= eigenvectors[0].dot(eigenvectors[2])*eigenvectors[0] + eigenvectors[1].dot(eigenvectors[2])*eigenvectors[1];
    eigenvectors[2].normalize();

    cout << "Known vector computation time: " << (clock()-t)/(double)CLOCKS_PER_SEC << endl;


    // matrix is symmetric so we can be less than careful about getting this right...
    int nz = covariance.getNumNonzero();
    int n = 2*numVertices;

    // Convert Ax, just in case our matrix is made of floats
    T* AxTemp = covariance.getAx();
    for (int i = 0; i < nz; i++) Ax[i] = (double)AxTemp[i];

    Ai = covariance.getAi();
    Ap = covariance.getAp();

    cout << "Computing permutation..." << endl;
    t = clock();
    double Info[AMD_INFO];
    if (amd_order (n, Ap, Ai, Pfw, (double *) NULL, Info) < AMD_OK) {
        cout << "call to AMD failed" << endl;
        exit(1);
    }
    amd_control((double*)NULL);
    amd_info(Info);

    cout << "AMD time: " << (clock()-t)/(double)CLOCKS_PER_SEC << endl;
    cout << "Doing symbolic ldl..." << endl;;

    t = clock();
    ldl_symbolic (2*numVertices, Ap, Ai, Lp, Parent, Lnz, Flag, Pfw, Pinv);
    cout << "Symbolic time: " << (clock()-t)/(double)CLOCKS_PER_SEC << endl;

    int lnz = Lp[n];
    REALLOC_MEMORY(Li, LDL_int, lnz);
    REALLOC_MEMORY(Lx, double, lnz);
    cout << "lnz = " << lnz;

    cout << "Doing numeric ldl..." << endl;
    t = clock();
    ldl_numeric (2*numVertices, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Pattern, Flag, Pfw, Pinv);
    cout << "Numeric time: " << (clock()-t)/(double)CLOCKS_PER_SEC << endl;

    T minVal[4];
    int smallestIndices[4];

    t = clock();
    for (int i = 0; i < 4; i++) {
        minVal[i] = 1e10;
        smallestIndices[i] = -1;

        for (int j = 0; j < n; j++)
            if (minVal[i] > D[j]) {
                bool bad = false;
                for (int k = 0; k < i; k++)
                    if (smallestIndices[k] == j) bad = true;
                if (!bad) {
                    smallestIndices[i] = j;
                    minVal[i] = D[j];
                }
            }
    }

    for (int i = 0; i < 3; i++)
        D[ smallestIndices[i] ] = 0;
    cout << "Zero out smallest time: " << (clock()-t)/(double)CLOCKS_PER_SEC << endl;

    t = clock();
    Vector old;
    for (unsigned int i = 3; i < eigenvectors.size(); i++) {
        // compute eigenvector orthogonal to first i eigenvectors

        // use old eigenvector and orthogonalize to start
        for (unsigned int j = 0; j < i; j++)
            eigenvectors[i] -= eigenvectors[i].dot(eigenvectors[j])*eigenvectors[j];

        if (eigenvectors[i].norm() < 1e-2) {
            eigenvectors[i].Random();
            for (unsigned int j = 0; j < i; j++)
                eigenvectors[i] -= eigenvectors[i].dot(eigenvectors[j])*eigenvectors[j];
        }

        eigenvectors[i].normalize();

        old = eigenvectors[i];
        double error = 5000;
        int numIterations = 0;

        cout << "Starting norm: " << eigenvectors[i].norm() << endl;
        while (error > 1e-5) {
            numIterations++;
            // solve with result on rhs
            for (int j = 0; j < 2*numVertices; j++)
                X[j] = eigenvectors[i][j];

            // the factorization is LDL' = PAP'
            LDL_perm (n, Y, X, Pfw) ;			// y = Pb
            LDL_lsolve (n, Y, Lp, Li, Lx) ;		// y = L\y
            LDL_dsolve (n, Y, D) ;			// y = D\y
            LDL_ltsolve (n, Y, Lp, Li, Lx) ;		// y = L'\y
            LDL_permt (n, X, Y, Pfw) ;			// x = P'y

            old = eigenvectors[i];
            for (int j = 0; j < 2*numVertices; j++)
                eigenvectors[i][j] = X[j];

            for (unsigned int j = 0; j < i; j++)
                eigenvectors[i] -= eigenvectors[i].dot(eigenvectors[j])*eigenvectors[j];

            eigenvectors[i].normalize();
            error = (old-eigenvectors[i]).norm();
            old = eigenvectors[i];
        }
        if (eigenvectors[i][0] < 0) eigenvectors[i] *= -1;
        cout << i << " Eigenvector took iterations" << numIterations << endl;
    }

    double itTime = (clock()-t)/(double)CLOCKS_PER_SEC;
    cout << "Iteration time: " << itTime << endl;
}
////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::getP(SimpleSparseMatrix<T> &prod)
{
    P2.reshape(numFaces*3, numFaces*4, 4*numFaces);
    dx2.reshape(numFaces, numVertices, 3*numFaces);
    dy2.reshape(numFaces, numVertices, 3*numFaces);

    // examine matrix push_back when in right order
    P2.startMatrixFill();
    dx2.startMatrixFill();
    dy2.startMatrixFill();
    for (int f = 0; f < numFaces; f++) {
        int i = faces[f][0], j = faces[f][1], k = faces[f][2];

        int temp;
        if (i > j) {
            temp = i;
            i = j;
            j = temp;
        }
        if (i > k) {
            temp = i;
            i = k;
            k = temp;
        }
        if (j > k) {
            temp = j;
            j = k;
            k = temp;
        }

        Vector2D<T> d1 = vertices[i] - vertices[k],
                    d2 = vertices[j] - vertices[i];

        T area = fabs(d1[1]*d2[0] - d1[0]*d2[1]);
        Vector2D<T> c1(-d1[1]/area,d1[0]/area),
                 c2(-d2[1]/area,d2[0]/area);

        dx2.addElement(f,i,-c1[0] - c2[0]);
        dx2.addElement(f,j,c1[0]);
        dx2.addElement(f,k,c2[0]);

        dy2.addElement(f,i,-c1[1] - c2[1]);
        dy2.addElement(f,j,c1[1]);
        dy2.addElement(f,k,c2[1]);

        P2.addElement(3*f, f, 2);
        P2.addElement(3*f+1, f+numFaces, SQRT_2);
        P2.addElement(3*f+1, f+2*numFaces, SQRT_2);
        P2.addElement(3*f+2, f+3*numFaces, 2);
    }

    int colShift[4] = {0, 0, numVertices, numVertices};

    SimpleSparseMatrix<T> *list[4] = {&dx2, &dy2, &dx2, &dy2};

    stacked.stack(list, 4, colShift);

    if (pCount == 0) P2.multiply(stacked, prod);
    else P2.parallelMultiply(stacked,prod);

    pCount++;
}
////////////////////////////////////////////////////////////////////////////////
#ifdef CSV

template<class T>
void KillingVectorsDeformation<T>::replacePoints(const QString &filename)   //todo: handle obj file too
{
    ifstream infile(filename.toAscii());

    if (filename.endsWith("off")) {
        string temp;
        infile >> temp;
        infile >> numVertices >> numFaces >> temp;

        qWarning("Mesh:  %d vertices, %d faces", numVertices, numFaces);

        vertices.resize(numVertices);
        T z;
        for (int i = 0; i < numVertices; i++)
            infile >> vertices[i].x >> vertices[i].y >> z;
    }

    if (filename.endsWith("obj")) {
        numVertices = 0;
        numFaces = 0;
        vertices.clear();
        T x,y,z;

        while (!infile.eof()) {
            // get line
            string curLine;
            getline(infile, curLine);

            // read type of the line
            istringstream issLine(curLine);
            string linetype;
            issLine >> linetype;

            if (linetype == "v") {
                numVertices++;
                issLine >> x >> y >> z;
                Point2D<T> p(x,y);
                vertices.push_back(p);
                continue;
            }
            if (linetype == "f") {
                numFaces++;
                continue;
            }
        }
        qWarning("numVertices = %d, numFaces = %d", numVertices, numFaces);
    }
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::saveVertices(ofstream& outfile, const QString &filename)
{

    if (filename.endsWith("off")) {
        for (int i = 0; i < numVertices; i++)
            outfile << vertices[i][0] << ' ' << vertices[i][1] << " 0\n";
    }

    if (filename.endsWith("obj")) {
        if (modelType == 2 && mtlFile != "") outfile << "mtllib " << mtlFile << endl;
        for (int i = 0; i < numVertices; i++)
            outfile << "v " << vertices[i][0] << ' ' << vertices[i][1] << " 0\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::saveTextureUVs(ofstream& outfile, const QString &filename)
{

    if (filename.endsWith("obj")) {
        for (int i = 0; i < numVertices; i++)
            outfile << "vt " << texCoords[i][0] << ' ' << texCoords[i][1] << endl;
        if (modelType == 2 && mtlFile != "") outfile << "usemtl Explorer_Default" << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::saveFaces(ofstream& outfile, const QString &filename)
{

    if (filename.endsWith("off")) {
        for (int i = 0; i < numFaces; i++)
            outfile << "3 " << faces[i][0] << ' ' << faces[i][1] << ' ' << faces[i][2] << endl;
    }

    if (filename.endsWith("obj")) {
        for (int i = 0; i < numFaces; i++)
            outfile << "f " << faces[i][0]+1 << '/' << faces[i][0]+1 << ' ' << faces[i][1]+1 << '/' << faces[i][1]+1 << ' ' << faces[i][2]+1 << '/' << faces[i][2]+1 << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::renderVertex(T left, T bottom, T meshWidth, T width, T height, Point2D<T> p)
{
    glLineWidth(LINE_WIDTH);
    T right = left + meshWidth;
    T meshHeight = (maxY - minY)*meshWidth/(maxX-minX);
    T top = bottom + meshHeight;

    T wFrac = (right-left)/width;
    T totWidth = (maxX - minX)/wFrac;
    T lowX = minX - totWidth * left / width;

    T hFrac = (top-bottom)/height;
    T totHeight = (maxY - minY)/hFrac;
    T lowY = minY - totHeight * bottom / height;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(lowX, lowX+totWidth, lowY, lowY+totHeight, 0, 1);
    glMatrixMode(GL_MODELVIEW);

    float s=totWidth / 500 * POINT_SIZE_SCALE;
    glBegin(GL_QUADS);
    glVertex2f(p[0]-s,p[1]-s);
    glVertex2f(p[0]-s,p[1]+s);
    glVertex2f(p[0]+s,p[1]+s);
    glVertex2f(p[0]+s,p[1]-s);
    glEnd(/*GL_QUADS*/);
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::renderSelectedVertex(T left, T bottom, T meshWidth, T width, T height, int v)
{
    glLineWidth(LINE_WIDTH);
    T right = left + meshWidth;
    T meshHeight = (maxY - minY)*meshWidth/(maxX-minX);
    T top = bottom + meshHeight;

    T wFrac = (right-left)/width;
    T totWidth = (maxX - minX)/wFrac;
    T lowX = minX - totWidth * left / width;

    T hFrac = (top-bottom)/height;
    T totHeight = (maxY - minY)/hFrac;
    T lowY = minY - totHeight * bottom / height;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(lowX, lowX+totWidth, lowY, lowY+totHeight, 0, 1);
    glMatrixMode(GL_MODELVIEW);

    Point2D<T> p = vertices[v];
    float s=totWidth / 500 * POINT_SIZE_SCALE;
    glBegin(GL_QUADS);
    glVertex2f(p[0]-s,p[1]-s);
    glVertex2f(p[0]-s,p[1]+s);
    glVertex2f(p[0]+s,p[1]+s);
    glVertex2f(p[0]+s,p[1]-s);
    glEnd(/*GL_QUADS*/);

    if (drawVFMode && vf.size() == numVertices && vfOrig.size() == numVertices) {
        double totalNorm = 0;

        for (int i = 0; i < numVertices; i++)
            totalNorm += vfOrig[i].normSquared();

        totalNorm = sqrt(totalNorm);

        glColor3f(0,0,1);
        glBegin(GL_LINES);
        glVertex2f(vertices[v][0], vertices[v][1]);
        glVertex2f(vertices[v][0]+vfOrig[v][0]/totalNorm*VF_SCALE, vertices[v][1]+vfOrig[v][1]/totalNorm*VF_SCALE);
        glEnd();

        totalNorm = 0;

        for (int i = 0; i < numVertices; i++)
            totalNorm += vf[i].normSquared();

        totalNorm = sqrt(totalNorm);

        glColor3f(0,.5,0);
        glBegin(GL_LINES);
        glVertex2f(vertices[v][0], vertices[v][1]);
        glVertex2f(vertices[v][0]+vf[v][0]/totalNorm*VF_SCALE, vertices[v][1]+vf[v][1]/totalNorm*VF_SCALE);
        glEnd();

        glColor3f(1,0,0);
    }
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
void KillingVectorsDeformation<T>::render(T left,T bottom,  T meshWidth, T width, T height)
{
    glLineWidth(LINE_WIDTH);
    T right = left + meshWidth;
    T meshHeight = (maxY - minY)*meshWidth/(maxX-minX);
    T top = bottom + meshHeight;

    T wFrac = (right-left)/width;
    T totWidth = (maxX - minX)/wFrac;
    T lowX = minX - totWidth * left / width;

    T hFrac = (top-bottom)/height;
    T totHeight = (maxY - minY)/hFrac;
    T lowY = minY - totHeight * bottom / height;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(lowX, lowX+totWidth, lowY, lowY+totHeight, 0, 1);
    glMatrixMode(GL_MODELVIEW);


    glColor3f(1,1,1);
    glEnable(GL_TEXTURE_2D);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // GL_LINE
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < numFaces; i++)
        for (int j = 0; j < 3; j++) {
            glTexCoord2f(texCoords[ faces[i][j] ][0],texCoords[ faces[i][j] ][1]);
            glVertex2f(vertices[ faces[i][j] ][0], vertices[ faces[i][j] ][1]);
        }
    glEnd(/*GL_TRIANGLES*/);
    glDisable(GL_TEXTURE_2D);

    //wireframe overlay
    glColor4f(0,0,0,wireframeTrans);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < numFaces; i++)
        for (int j = 0; j < 3; j++) {
            glVertex2f(vertices[ faces[i][j] ][0], vertices[ faces[i][j] ][1]);
        }
    glEnd();

    if (drawVFMode && vf.size() == numVertices && vfOrig.size() == numVertices) {
        double totalNorm = 0;

        for (int i = 0; i < numVertices; i++)
            totalNorm += vfOrig[i].normSquared();

        totalNorm = sqrt(totalNorm);

        glColor3f(0,0,1);
        glBegin(GL_LINES);
        for (int i = 0; i < numVertices; i++) {
            glVertex2f(vertices[i][0], vertices[i][1]);
            glVertex2f(vertices[i][0]+vfOrig[i][0]/totalNorm*VF_SCALE, vertices[i][1]+vfOrig[i][1]/totalNorm*VF_SCALE);
        }
        glEnd();

        totalNorm = 0;

        for (int i = 0; i < numVertices; i++)
            totalNorm += vf[i].normSquared();

        totalNorm = sqrt(totalNorm);

        glColor3f(0,.5,0);
        glBegin(GL_LINES);
        for (int i = 0; i < numVertices; i++) {
            glVertex2f(vertices[i][0], vertices[i][1]);
            glVertex2f(vertices[i][0]+vf[i][0]/totalNorm*VF_SCALE, vertices[i][1]+vf[i][1]/totalNorm*VF_SCALE);
        }
        glEnd();
    }

    /*glColor4f(0,0,0,1);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < numFaces; i++)
        for (int j = 0; j < 3; j++) {
            glVertex2f(vertices[ faces[i][j] ][0], vertices[ faces[i][j] ][1]);
        }
    glEnd();
    */
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
int KillingVectorsDeformation<T>::getClosestVertex(Point2D<T> point, T dist)   // linear time -- could make faster...
{
    int closest = -1;
    T closestDistance = numeric_limits<T>::max();

    for (int i = 0; i < numVertices; i++) {
        T distance = vertices[i].distanceSquared(point);

        if (distance < closestDistance && distance < dist) {
            closestDistance = distance;
            closest = i;
        }
    }

    return closest;
}
////////////////////////////////////////////////////////////////////////////////

/*
template<class T>
void KillingVectorsDeformation<T>::changePinnedVertices(set<int> &pinned)
{
}
*/
////////////////////////////////////////////////////////////////////////////////

template class KillingVectorsDeformation<float>;
template class KillingVectorsDeformation<double>;
#endif

#endif
