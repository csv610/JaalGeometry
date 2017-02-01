#pragma once

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>

#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include "basic_math.hpp"

using namespace std;

template<class T>
class JSparseMatrixRep
{
    typedef T& reference;
    typedef const T& const_reference;
public:
    typedef typename std::map<size_t, T> SparseRow; // Seems to be the fastest.
    //  typedef typename __gnu_cxx::hash_map<size_t, T> SparseRow;
    //  typedef tr1::unordered_map<size_t, T> SparseRow;
    typedef typename SparseRow::const_iterator const_row_iterator;

    JSparseMatrixRep() {
        nrows = 0;
        ncols = 0;
        default_value = 0;
    }

    int setSize(size_t m, size_t n) {
        nrows = m;
        ncols = n;
        if (nrows) data.resize(nrows);
        return 0;
    }

    void clear() {
        data.clear();
        nrows = 0;
        ncols = 0;
    }

    SparseRow &getRow(size_t irow) {
        return data[irow];
    }

    const SparseRow &getRow(size_t irow) const {
        return data[irow];
    }

    void clearRow(size_t irow) {
        data[irow].clear();
    }

    /**
     * Set the value at (i,j) position. If the value is zero, it will not
     * be added in the matrix and if the exisiting value is again set to
     * zero, the entry will be removed.
     */
    void setValue(size_t i, size_t j, const T &val) {
        if ( val == default_value ) {
            erase(i,j);
            return;
        }

        if (i >= nrows || j >= ncols) {
            cout << "Warning: Row Index is out of range " << i << endl;
            return;
        }

        data[i][j] = val;

    }

    bool isZero(size_t i, size_t j) const {
        typename SparseRow::const_iterator it;
        it = data[i].find(j);
        if (it == data[i].end()) return 1;
        return 0;
    }

    size_t getDegree(size_t irow) const {
        size_t degree = data[irow].size();
        if (data[irow].find(irow) != data[irow].end()) degree -= 1;
        return degree;
    }

    vector<size_t> getDegrees() const {
        vector<size_t> degree;
        degree.resize(nrows);
        for (size_t i = 0; i < nrows; i++) degree[i] = getDegree(i);
        return degree;
    }

    size_t getNumNonZeros() {
        size_t sum = 0;
        for (size_t i = 0; i < nrows; i++) sum += data[i].size();
        return sum;
    }

    const_reference at(size_t i, size_t j) const {
        if (i >= nrows || j >= ncols) {
            cout << "Fatal Error: Matrix Index out of range " << endl;
            exit(0);
        }
        typename SparseRow::const_iterator it;
        it = data[i].find(j);
        if (it == data[i].end()) return default_value;
        return it->second;
    }

    reference at(size_t i, size_t j) {
        if (i >= nrows || j >= ncols) {
            cout << "Fatal Error: Matrix Index out of range " << endl;
            exit(0);
        }
        typename SparseRow::iterator it;
        it = data[i].find(j);
        if (it == data[i].end()) data[i][j] = default_value;
        return data[i][j];
    }

    void erase( size_t i, size_t j) {
        data[i].erase(j);
    }

    void prune() {
        for( size_t i = 0; i < nrows; i++) {
            while(1) {
                bool finished = 1;
                for( auto keyVal:data[i] ) {
                    if( keyVal.second == default_value) {
                        erase(i, keyVal.first);
                        finished = 0;
                        break;
                    }
                }
                if( finished ) break;
            }
        }
    }

    int saveAs(ofstream &ofile) const {
        typename SparseRow::const_iterator it;
        for (size_t ir = 0; ir < nrows; ir++) {
            for (it = data[ir].begin(); it != data[ir].end(); ++it) {
                ofile << ir << " " << it->first << " " << it->second << endl;
            }
        }
        return 0;
    }

    void printOut(ostream &s)
    {
        typename SparseRow::const_iterator it;
        for (size_t i = 0; i < nrows; i++) {
            for (it = data[i].begin(); it != data[i].end(); ++it)
                s << i << " " << it->first << " " << it->second << endl;
        }
    }
    T  default_value;
    size_t nrows, ncols;
    vector<SparseRow> data;
};


template<class T>
class JGeneralSparseMatrix
{
    typedef JSparseMatrixRep<T> MatRep;
public:
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef T& reference;
    typedef const T& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef typename MatRep::SparseRow SparseRow;
    typedef typename MatRep::const_row_iterator const_row_iterator;

    JGeneralSparseMatrix(size_t m = 0, size_t n = 0) {
        rep = boost::shared_ptr<MatRep> (new MatRep);
        assert(rep);
        rep->setSize(m, n);
    }

    JGeneralSparseMatrix<T> getClone() const {
        JGeneralSparseMatrix<T> sm;
        sm.rep->nrows = rep->nrows;
        sm.rep->ncols = rep->ncols;
        sm.rep->data  = rep->data;
        return sm;
    }

    int setSize(size_t m, size_t n) {
        return rep->setSize(m, n);
    }

    bool isEmpty() const {
        bool val = rep->getNumNonZeros() == 0 ? 1 : 0;
        return val;
    }

    size_t numRows() const {
        return rep->nrows;
    }

    size_t numCols() const {
        return rep->ncols;
    }

    size_t getNumNonZeros() {
        return rep->getNumNonZeros();
    }

    bool isZero(size_t i, size_t j) const {
        return rep->isZero(i, j);
    }

    size_t getDegree(size_t irow) const {
        return rep->getDegree(irow);
    }

    vector<size_t> getDegrees() const {
        return rep->getDegrees();
    }

    void clear() {
        rep->clear();
    }

    SparseRow &getRow(size_t irow) {
        return rep->getRow(irow);
    }

    const SparseRow &getRow(size_t irow) const {
        return rep->getRow(irow);
    }

    void clearRow(size_t irow) {
        rep->clearRow(irow);
    }

    void prune() {
        rep->prune();
    }

    void setValue(size_t i, size_t j, const T &val) {
        rep->setValue(i, j, val);
    }

    reference operator() (size_t i, size_t j) {
        return rep->at(i, j);
    }

    const_reference operator() (size_t i, size_t j) const {
        return rep->at(i, j);
    }

    void erase( size_t i, size_t j) {
        return rep->erase(i,j);
    }


    void
    get_compressed_row_data(vector<T> &buf, vector<size_t> &colindex,
                            vector<size_t>& offset) const {
        rep->get_compressed_row_data(buf, colindex, offset);
    }

    template<class T1>
    friend ostream & operator <<(ostream &s, const JGeneralSparseMatrix<T1> &sm);

    int saveAs(const string &fname) const {
        ofstream ofile(fname.c_str(), ios::out);
        if (ofile.fail()) {
            cout << "Warning: Cann't open file " << fname << endl;
            return 1;
        }
        ofile << "G " << rep->nrows << " " << rep->ncols << endl;
        rep->saveAs(ofile);
        return 0;
    }
private:
    boost::shared_ptr<MatRep> rep;
};

template<class T>
ostream & operator <<(ostream &s, const JGeneralSparseMatrix<T> &sm)
{
    if( sm.rep == nullptr) return s;
    sm.rep->printOut(s);
    return s;
}

template<class T>
class JSymmetricSparseMatrix
{
    typedef JSparseMatrixRep<T> MatRep;
public:
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef T& reference;
    typedef const T& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef typename MatRep::SparseRow SparseRow;
    typedef typename MatRep::const_row_iterator const_row_iterator;

    JSymmetricSparseMatrix(size_t m = 0, size_t n = 0) {
        rep = boost::shared_ptr<MatRep > (new MatRep);
        assert(rep);
        setSize(m, n);
    }

    bool isEmpty() const {
        bool val = rep->getNumNonZeros() == 0 ? 1 : 0;
        return val;
    }

    size_t getDegree(size_t irow) const {
        return rep->getDegree(irow);
    }

    vector<size_t> getDegrees() const {
        return rep->getDegrees();
    }

    int setSize(size_t m, size_t n) {
        return rep->setSize(m, n);
    }

    size_t numRows() const {
        return rep->nrows;
    }

    size_t getCols() const {
        return rep->nrows;
    }

    size_t getNumNonZeros() const {
        return rep->getNumNonZeros();
    }

    void prune() {
        rep->prune();
    }

    void clear() {
        rep->clear();
    }

    SparseRow & getRow(size_t irow) {
        return rep->getRow(irow);
    }

    const SparseRow & getRow(size_t irow) const {
        return rep->getRow(irow);
    }

    void clearRow(size_t irow) {
        rep->clearRow(irow);
    }

    void setValue(size_t i, size_t j, const T &val) {
        if (j >= i)
            rep->setValue(i, j, val);
        else
            rep->setValue(j, i, val);
    }

    bool isZero(size_t i, size_t j) const {
        if (j >= i) return rep->isZero(i, j);
        rep->isZero(j, i);
    }


    reference operator() (size_t i, size_t j) {
        if (i >= j) return rep->at(i, j);
        return rep->at(j, i);
    }

    const_reference operator() (size_t i, size_t j) const {
        if (i >= j) return rep->at(i, j);
        return rep->at(j, i);
    }

    void erase( size_t i, size_t j) {
        rep->erase(i,j);
        rep->erase(j,i);
    }

    JSymmetricSparseMatrix<T> getClone() {
        JSymmetricSparseMatrix<T> sm;
        sm.rep->nrows = rep->nrows;
        sm.rep->ncols = rep->ncols;
        sm.rep->data = rep->data;
        return sm;
    }

    template<class T1>
    friend ostream & operator <<(ostream &s, const JSymmetricSparseMatrix<T1> &sm);

    int saveAs(const string &fname) const {
        ofstream ofile(fname.c_str(), ios::out);
        if (ofile.fail()) {
            cout << "Warning: Cann't open file " << fname << endl;
            return 1;
        }
        ofile << "S " << rep->nrows << " " << rep->ncols << endl;
        rep->saveAs(ofile);
        return 0;
    }

private:
    boost::shared_ptr<MatRep> rep;
};

template<class T>
ostream & operator <<(ostream &s, const JSymmetricSparseMatrix<T> &sm)
{

//  s << sm.rep;
    return s;
}

#ifdef CSV

template<class T>
class GeneralDenseMatrix : public Matrix<T>
{
    typedef DenseMatrixRep<T> MatRep;
public:
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef T& reference;
    typedef const T& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    GeneralDenseMatrix() {
        rep = boost::shared_ptr<MatRep > (new MatRep);
    }

    GeneralDenseMatrix(size_t n, size_t m) {
        rep = boost::shared_ptr<MatRep > (new MatRep);
        setSize(n, m);
    }

    MatrixType getMatrixType() const {
        return GENERAL_DENSE_MATRIX;
    }

    bool isEmpty() const {
        if (getSize1() == 0) return 1;
        if (getSize2() == 0) return 1;
        return 0;
    }

    int setSize(size_t n, size_t m) {
        assert(rep);
        rep->data.resize(n);
        for (size_t i = 0; i < n; i++)
            rep->data[i].resize(m);
        rep->nrows = n;
        rep->ncols = m;
    }

    size_t getSize1() const {
        return rep->nrows;
    }

    size_t getSize2() const {
        return rep->ncols;
    }

    void clear() {
        rep->clear();
    }

    //! Overloading operator

    reference operator() (size_t i, size_t j) {
        return rep->at(i, j);
    }

    //! Overloading operator

    const_reference operator()(size_t i, size_t j) const {
        return rep->at(i, j);
    }

    vector<T> getRowData() const {
        size_t nrows = rep->nrows;
        size_t ncols = rep->ncols;
        vector<T> tmpdata(nrows * ncols);
        size_t indx = 0;
        for (size_t i = 0; i < nrows; i++)
            for (size_t j = 0; j < ncols; j++)
                tmpdata[indx++] = rep->data[i][j];
        return tmpdata;
    }

    vector<T> getColumnData() const {
        size_t nrows = rep->nrows;
        size_t ncols = rep->ncols;
        vector<T> tmpdata(nrows * ncols);
        size_t indx = 0;
        for (size_t j = 0; j < ncols; j++)
            for (size_t i = 0; i < nrows; i++)
                tmpdata[indx++] = rep->data[i][j];
        return tmpdata;
    }

    template<class T1>
    friend ostream & operator <<(ostream &s, const GeneralDenseMatrix<T1> &rp);

    void
    saveAs(const string &s) const {
        rep->saveAs(s);
    }

private:
    boost::shared_ptr<MatRep> rep;
};

template<class T>
ostream & operator <<(ostream &s, const GeneralDenseMatrix<T> &mat)
{
    size_t nrows = mat.getSize1();
    size_t ncols = mat.getSize2();
    s << nrows << " " << ncols << endl;
    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++)
            s << mat(i, j) << " ";
        s << endl;
    }
    return s;
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
class SymmetricDenseMatrix : public Matrix<T>
{
    typedef DenseMatrixRep<T> MatRep;
public:
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef T& reference;
    typedef const T& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    SymmetricDenseMatrix() {
        rep = boost::shared_ptr<MatRep > (new MatRep);
    }

    SymmetricDenseMatrix(size_t n, size_t m) {
        rep = boost::shared_ptr<MatRep > (new MatRep);
        setSize(n, m);
    }

    MatrixType getMatrixType() const {
        return SYMMETRIC_DENSE_MATRIX;
    }

    bool isEmpty() const {
        if (getSize1() == 0) return 1;
        if (getSize2() == 0) return 1;
        return 0;
    }

    //! Set the matrix size.

    int setSize(size_t n, size_t m) {
        assert(n == m);
        rep->data.resize(n);
        for (size_t i = 0; i < n; i++)
            rep->data[i].resize(i + 1);

        rep->nrows = n;
        rep->ncols = m;
    }

    //! get  the size of matrix in a given direction: Row: 0, Column 1

    size_t getSize1() const {
        return rep->nrows;
    }

    size_t getSize2() const {
        return rep->ncols;
    }

    void clear() {
        rep->clear();
    }

    reference operator() (size_t i, size_t j) {
        if (j > i) return rep->at(j, i);
        return rep->at(i, j);
    }

    const_reference operator() (size_t i, size_t j) const {
        if (j > i) return rep->at(j, i);
        return rep->at(i, j);
    }

    //! Return the lower triangulation: See the Lapack Document

    vector<T> getRowData() const {
        size_t nrows = rep->nrows;
        vector<T> tmpdata;
        tmpdata.resize(nrows * (nrows + 1) / 2);
        size_t index = 0;
        for (size_t j = 0; j < nrows; j++)
            for (size_t i = j; i < j + 1; i++)
                tmpdata[index++] = rep->at(i, j);
        return tmpdata;
    }

    //! Return the upper triangulation : See the Lapack Document

    vector<T> getColumnData() const {
        size_t nrows = rep->nrows;
        size_t ncols = rep->ncols;
        vector<T> tmpdata;
        tmpdata.resize(nrows * (nrows + 1) / 2);
        size_t index = 0;
        for (size_t j = 0; j < ncols; j++)
            for (size_t i = j; i < nrows; i++)
                tmpdata[index++] = rep->at(i, j);
        return tmpdata;
    }

    template<class T1>
    friend ostream & operator <<(ostream &s, const SymmetricDenseMatrix<T1> &rp);

    void saveAs(const string &s) const {
        ofstream ofile(s.c_str(), ios::out);
        if (ofile.fail()) {
            cout << "Warning: Cann't open file " << endl;
            return;
        }
        size_t nrows = getSize1();
        size_t ncols = getSize2();

        ofile << nrows << " " << ncols << endl;
        for (size_t i = 0; i < nrows; i++) {
            for (size_t j = 0; j < ncols; j++)
                ofile << rep->at(i, j) << " ";
            ofile << endl;
        }
    }
private:
    boost::shared_ptr<MatRep> rep;
};

template<class T>
ostream & operator <<(ostream &s, const SymmetricDenseMatrix<T> &mat)
{
    size_t nrows = mat.getSize1();
    size_t ncols = mat.getSize2();
    s << nrows << " " << ncols << endl;
    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++)
            s << mat(i, j) << " ";
        s << endl;
    }
    return s;
}

//////////////////////////////////////////////////////////////////////////////
//

template< class T>
GeneralSparseMatrix<T> transpose(const GeneralSparseMatrix<T> &a)
{
    GeneralSparseMatrix<T> result;

    size_t numRows = a.getSize1();
    if (numRows == 0) return result;

    result.setSize(numRows);

    typename SparseMatrixRep<T>::SparseRow arow;
    typename SparseMatrixRep<T>::const_row_iterator jt;

    for (size_t i = 0; i < numRows; i++) {
        arow = a.getConstRow(i);
        for (jt = arow.begin(); jt != arow.end(); ++jt) {
            size_t j = jt->first;
            result.setValue(j, i, jt->second);
        }
    }
    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
bool isSymmetric(GeneralSparseMatrix<T> &a)
{
    typename SparseMatrixRep<T>::SparseRow arow;
    typename SparseMatrixRep<T>::const_row_iterator jt;

    size_t numRows = a.getSize1();
    for (size_t i = 0; i < numRows; i++) {
        arow = a.getConstRow(i);
        for (jt = arow.begin(); jt != arow.end(); ++jt) {
            size_t j = jt->first;
            if ((j > i) && (a(j, i) != a(i, j))) return 0;
        }
    }
    return 1;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
GeneralSparseMatrix<T> multiply(const GeneralSparseMatrix<T> &a,
                                const GeneralSparseMatrix<T> &b)
{
    size_t nrows = a.getSize1();
    assert(nrows == b.getSize1());

    GeneralSparseMatrix<T> result(nrows);

    typename SparseMatrixRep<T>::SparseRow arow, brow;
    typename SparseMatrixRep<T>::const_row_iterator jt, kt;

    for (size_t i = 0; i < nrows; i++) {
        arow = a.getConstRow(i);
        for (jt = arow.begin(); jt != arow.end(); ++jt) {
            size_t j = jt->first;
            T aval = jt->second;
            brow = b.getConstRow(j);
            for (kt = brow.begin(); kt != brow.end(); ++kt) {
                size_t k = kt->first;
                T bval = kt->second;
                T cval = result(i, k) + aval*bval;
                result.setValue(i, k, cval);
            }
        }
    }
    return result;
}
//////////////////////////////////////////////////////////////////////////////

template<class T>
GeneralSparseMatrix<T> multiply(const SymmetricSparseMatrix<T> &a,
                                const SymmetricSparseMatrix<T> &b)
{
    size_t nrows = a.getSize1();
    assert(nrows == b.getSize1());

    GeneralSparseMatrix<T> result(nrows);

    typename SparseMatrixRep<T>::SparseRow arow, brow;
    typename SparseMatrixRep<T>::const_row_iterator jt, kt;

    for (size_t i = 0; i < nrows; i++) {
        arow = a.getConstRow(i);
        for (jt = arow.begin(); jt != arow.end(); ++jt) {
            size_t j = jt->first;
            T aval = jt->second;
            brow = b.getConstRow(j);
            for (kt = brow.begin(); kt != brow.end(); ++kt) {
                size_t k = kt->first;
                T bval = kt->second;
                T cval = result(i, k) + aval*bval;
                result.setValue(i, k, cval);
            }
        }
    }
    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class T1, class T2>
vector<T2> multiply(const GeneralSparseMatrix<T1> &a, const vector<T2> &v)
{
    size_t nrows = a.getSize1();
    assert(nrows == v.size());

    vector<T2> result(nrows);

    typename SparseMatrixRep<T1>::SparseRow row;
    typename SparseMatrixRep<T1>::const_row_iterator it;

    for (size_t i = 0; i < nrows; i++) {
        row = a.getConstRow(i);
        T2 sum = 0;
        for (it = row.begin(); it != row.end(); ++it) {
            size_t colindex = it->first;
            T1 val = it->second;
            sum += val * v[colindex];
        }
        result[i] = sum;
    }
    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class T1, class T2, class T3>
void multiply(const SymmetricSparseMatrix<T1> &a, const T2 *x, T3 *y)
{
    size_t nrows = a.getSize1();

    typename SparseMatrixRep<T1>::SparseRow sparse_row;
    typename SparseMatrixRep<T1>::const_row_iterator it;

    for (size_t i = 0; i < nrows; i++) y[i] = 0;

    for (size_t i = 0; i < nrows; i++) {
        sparse_row = a.getConstRow(i);
        for (it = sparse_row.begin(); it != sparse_row.end(); ++it) {
            size_t j = it->first;
            T1 coeff = it->second;
            y[i] += coeff * x[j];
            if (i != j) y[j] += coeff * x[i];
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<class T1, class T2>
vector<T2> multiply(const SymmetricSparseMatrix<T1> &a, const vector<T2> &x)
{
    size_t nrows = a.getSize1();
    assert(nrows == x.size());
    vector<T2> y(nrows);

    Math::multiply(a, &x[0], &y[0]);
    return y;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
GeneralSparseMatrix<T> multiply(const GeneralSparseMatrix<T> &a,
                                const SymmetricSparseMatrix<T> &b)
{
    size_t nrows = a.getSize1();
    assert(nrows == b.getSize1());

    GeneralSparseMatrix<T> result(nrows);

    typename GeneralSparseMatrix<T>::SparseRow arow, brow;
    typename GeneralSparseMatrix<T>::const_row_iterator jt, kt;

    for (size_t i = 0; i < nrows; i++) {
        arow = a.getConstRow(i);
        for (jt = arow.begin(); jt != arow.end(); ++jt) {
            size_t j = jt->first;
            T aval = jt->second;
            brow = b.getConstRow(j);
            for (kt = brow.begin(); kt != brow.end(); ++kt) {
                size_t k = kt->first;
                T bval = kt->second;
                T cval = result(i, k) + aval*bval;
                result.setValue(i, k, cval);
            }
        }
    }
    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
GeneralSparseMatrix<T> general_sparse_matrix(SymmetricSparseMatrix<T> &a,
        bool deleteOriginal = 0)
{
    size_t nrows = a.getSize1();
    size_t ncols = a.getSize2();

    GeneralSparseMatrix<T> result;

    result.setSize(nrows, ncols);

    typename GeneralSparseMatrix<T>::SparseRow row;
    typename GeneralSparseMatrix<T>::const_row_iterator jt;

    for (size_t i = 0; i < nrows; i++) {
        row = a.getConstRow(i);
        for (jt = row.begin(); jt != row.end(); ++jt) {
            size_t j = jt->first;
            T val = jt->second;
            result.setValue(i, j, val);
            result.setValue(j, i, val);
        }
        if (deleteOriginal) a.clearRow(i);
    }

    if (deleteOriginal) a.clear(); // will make nrows = ncols = 0

    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
GeneralDenseMatrix<T> general_dense_matrix(GeneralSparseMatrix<T> &a,
        bool deleteOriginal = 0)
{
    size_t nrows = a.getSize1();
    size_t ncols = a.getSize2();

    GeneralDenseMatrix<T> result;
    result.setSize(nrows, ncols);

    typename GeneralSparseMatrix<T>::SparseRow row;
    typename GeneralSparseMatrix<T>::const_row_iterator jt;

    for (size_t i = 0; i < nrows; i++) {
        row = a.getConstRow(i);
        for (jt = row.begin(); jt != row.end(); ++jt) {
            size_t j = jt->first;
            double val = jt->second;
            result(i, j) = val;
        }
        if (deleteOriginal) a.clearRow(i);
    }

    if (deleteOriginal) a.clear(); // will make nrows = ncols = 0

    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
GeneralDenseMatrix<T> general_dense_matrix(SymmetricSparseMatrix<T> &a,
        bool deleteOriginal = 0)
{
    size_t nrows = a.getSize1();
    size_t ncols = a.getSize2();

    GeneralDenseMatrix<T> result;
    result.setSize(nrows, ncols);

    typename GeneralSparseMatrix<T>::SparseRow row;
    typename GeneralSparseMatrix<T>::const_row_iterator jt;

    for (size_t i = 0; i < nrows; i++) {
        row = a.getConstRow(i);
        for (jt = row.begin(); jt != row.end(); ++jt) {
            size_t j = jt->first;
            double val = jt->second;
            result(i, j) = val;
            result(j, i) = val;
        }
        if (deleteOriginal) a.clearRow(i);
    }

    if (deleteOriginal) a.clear(); // will make nrows = ncols = 0

    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
SymmetricSparseMatrix<T> symmetric_sparse_matrix(GeneralSparseMatrix<T> &a,
        bool deleteOriginal = 0)
{
    size_t nrows = a.getSize1();
    size_t ncols = a.getSize2();

    SymmetricSparseMatrix<T> result;

    result.setSize(nrows, ncols);

    typename GeneralSparseMatrix<T>::SparseRow row;
    typename GeneralSparseMatrix<T>::const_row_iterator jt;

    for (size_t i = 0; i < nrows; i++) {
        row = a.getConstRow(i);
        for (jt = row.begin(); jt != row.end(); ++jt) {
            size_t j = jt->first;
            result.setValue(i, j, jt->second);
        }
        if (deleteOriginal) a.clearRow(i);
    }

    if (deleteOriginal) a.clear(); // will make nrows = ncols = 0

    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class TDST, class TSRC>
int copy(GeneralSparseMatrix<TDST> &dst, GeneralSparseMatrix<TSRC> &src, bool deleteSrc = 0)
{
    size_t nrows = src.getSize1();
    size_t ncols = src.getSize2();

    dst.setSize(nrows, ncols);

    typename GeneralSparseMatrix<TSRC>::SparseRow row;
    typename GeneralSparseMatrix<TSRC>::const_row_iterator jt;

    for (size_t i = 0; i < nrows; i++) {
        row = src.getConstRow(i);
        for (jt = row.begin(); jt != row.end(); ++jt) {
            size_t j = jt->first;
            TDST val = static_cast<TDST> (jt->second);
            dst.setValue(i, j, val);
        }
        if (deleteSrc) src.clearRow(i);
    }

    if (deleteSrc) src.clear(); // will make nrows = ncols = 0

    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<class TDST, class TSRC>
int copy(SymmetricSparseMatrix<TDST> &dst, SymmetricSparseMatrix<TSRC> &src, bool deleteSrc)
{
    size_t nrows = src.getSize1();
    size_t ncols = src.getSize2();

    dst.setSize(nrows, ncols);

    typename SymmetricSparseMatrix<TSRC>::SparseRow row;
    typename SymmetricSparseMatrix<TSRC>::const_row_iterator jt;

    for (size_t i = 0; i < nrows; i++) {
        row = src.getConstRow(i);
        for (jt = row.begin(); jt != row.end(); ++jt) {
            size_t j = jt->first;
            TDST val = static_cast<TDST> (jt->second);
            dst.setValue(i, j, val);
        }
        if (deleteSrc) src.clearRow(i);
    }

    if (deleteSrc) src.clear(); // will make nrows = ncols = 0

    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
GeneralSparseMatrix<T> power(const GeneralSparseMatrix<T> &a, int ntimes)
{
    assert(ntimes > 0);
    if (ntimes == 1) return a;

    GeneralSparseMatrix<T> b, result;
    size_t numRows = a.getSize1();

    result = a.getClone();

    for (int i = 1; i < ntimes; i++) {
        b = multiply(result, a);
        result.clear();
        result = b;
    }

    return result;
}

//////////////////////////////////////////////////////////////////////////////

//! \brief Calculating power of Symmetric Sparse matrix...

template<class T>
SymmetricSparseMatrix<T> power(SymmetricSparseMatrix<T> &a, int ntimes)
{
    assert(ntimes > 0);
    if (ntimes == 1) return a;

    size_t numRows = a.getSize1();

    GeneralSparseMatrix<T> b, gm = general_sparse_matrix(a);

    for (int i = 1; i < ntimes; i++) {
        b = multiply(gm, a);
        gm.clear();
        gm = b;
    }

    typename GeneralSparseMatrix<T>::SparseRow arow;
    typename GeneralSparseMatrix<T>::const_row_iterator jt;

    SymmetricSparseMatrix<T> result;

    result.setSize(numRows, numRows);

    for (size_t i = 0; i < numRows; i++) {
        arow = gm.getConstRow(i);
        for (jt = arow.begin(); jt != arow.end(); ++jt) {
            size_t j = jt->first;
            if (j >= i) result.setValue(i, j, jt->second);
        }
        gm.clearRow(i);
    }
    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
T determinant33(const GeneralDenseMatrix<T> &a)
{
    assert(a.getSize1() == 3 && a.getSize2() == 3);
    T val;
    val = a(0, 0)*(a(1, 1) * a(2, 2) - a(2, 1) * a(1, 2))
          - a(0, 1)*(a(1, 0) * a(2, 2) - a(2, 0) * a(1, 2))
          + a(0, 2)*(a(1, 0) * a(2, 1) - a(2, 0) * a(1, 1));

    return val;
}

//###########################################################################

template <class T>
GeneralDenseMatrix<T> inverse_matrix33(const GeneralDenseMatrix<T> &mat)
{

    //*************************************************************
    // Objective :  Frequently matrix inversion of 3X3 are required.
    // The following is the fast matrix inversion. Generally this
    // approach is not used for large matrices.
    //
    // Matrix a =  a[0]  a[1]  a[2]
    //	       a[3]  a[4]  a[5]
    //	       a[6]  a[7]  a[8]
    //*************************************************************
    assert(mat.getSize1() == 3 && mat.getSize2() == 3);

    T Determinant, Divideby;

    Determinant = determinant33(mat);

    /*
    if( fabs((double)Determinant) < 1.0E-15) {
      cout << "Warning: Matrix is singular: Cann't find inverse " << Determinant << endl;
      abort();
    }
     */

    Divideby = 1.0 / (Determinant);

    T a[9];
    a[0] = mat(0, 0);
    a[1] = mat(0, 1);
    a[2] = mat(0, 2);
    a[3] = mat(1, 0);
    a[4] = mat(1, 1);
    a[5] = mat(1, 2);
    a[6] = mat(2, 0);
    a[7] = mat(2, 1);
    a[8] = mat(2, 2);

    T cofactor[9];
    cofactor[0] = Divideby * (a[4] * a[8] - a[5] * a[7]);
    cofactor[1] = -Divideby * (a[3] * a[8] - a[5] * a[6]);
    cofactor[2] = Divideby * (a[3] * a[7] - a[4] * a[6]);

    cofactor[3] = -Divideby * (a[1] * a[8] - a[2] * a[7]);
    cofactor[4] = Divideby * (a[0] * a[8] - a[2] * a[6]);
    cofactor[5] = -Divideby * (a[0] * a[7] - a[1] * a[6]);

    cofactor[6] = Divideby * (a[1] * a[5] - a[2] * a[4]);
    cofactor[7] = -Divideby * (a[0] * a[5] - a[2] * a[3]);
    cofactor[8] = Divideby * (a[0] * a[4] - a[1] * a[3]);

    GeneralDenseMatrix<T> result(3, 3);
    result(0, 0) = cofactor[0];
    result(0, 1) = cofactor[3];
    result(0, 2) = cofactor[6];

    result(1, 0) = cofactor[1];
    result(1, 1) = cofactor[4];
    result(1, 2) = cofactor[7];

    result(2, 0) = cofactor[2];
    result(2, 1) = cofactor[5];
    result(2, 2) = cofactor[8];

    return result;
}

template<class T>
void inverse_matrix2X2(T *s, T *result)
{
    //////////////////////////////////////////////////////////////
    // Matrix a =  a[0]  a[1]
    //	         a[2]  a[3]
    //////////////////////////////////////////////////////////////
    //
    T a = s[0];
    T b = s[1];
    T c = s[2];
    T d = s[3];

    T coeff = 1.0 / (a * d - b * c);

    result[0] = coeff*d;
    result[1] = -coeff*b;
    result[2] = -coeff*c;
    result[3] = coeff*a;
}

template<class T>
vector<T> multiply(const Math::GeneralDenseMatrix<T> &A, const vector<T> &x)
{
    size_t nRows = A.getSize1();
    size_t nCols = A.getSize2();
    assert(x.size() == nRows);

    vector<T> y;
    y.resize(nRows);
    for (size_t i = 0; i < nRows; i++) {
        T sum = 0;
        for (size_t j = 0; j < nCols; j++)
            sum += A(i, j) * x[j];
        y[i] = sum;
    }
    return y;
}

template<class T>
vector<T> multiply(const Math::SymmetricDenseMatrix<T> &A, const vector<T> &x)
{
    size_t nRows = A.getSize1();
    size_t nCols = A.getSize2();
    assert(x.size() == nRows);

    vector<T> y;
    y.resize(nRows);
    for (size_t i = 0; i < nRows; i++) {
        T sum = 0;
        for (size_t j = 0; j < nCols; j++)
            sum += A(i, j) * x[j];
        y[i] = sum;
    }
    return y;
}

template<class T>
void multiply(const Math::GeneralDenseMatrix<T> &A, const T *x, T *y)
{
    size_t nRows = A.getSize1();
    size_t nCols = A.getSize2();

    for (size_t i = 0; i < nRows; i++) {
        T sum = 0;
        for (size_t j = 0; j < nCols; j++)
            sum += A(i, j) * x[j];
        y[i] = sum;
    }
}

template<class T>
void transform(GeneralDenseMatrix<T> &dmatrix, const GeneralSparseMatrix<T> &smatrix)
{
    size_t i, j, nrows, ncols;

    nrows = dmatrix.getSize1();
    ncols = dmatrix.getSize2();

    smatrix.setSize(nrows, ncols);

    for (i = 0; i < nrows; i++)
        for (j = 0; j < ncols; j++)
            dmatrix(i, j) = smatrix(i, j);
}

template<class T>
void transform(GeneralSparseMatrix<T> &smatrix, const GeneralDenseMatrix<T> &dmatrix)
{
    size_t i, j, nrows, ncols;

    nrows = dmatrix.getSize1();
    ncols = dmatrix.getSize2();

    smatrix.setSize(nrows, ncols);

    for (i = 0; i < nrows; i++)
        for (j = 0; j < ncols; j++)
            smatrix(i, j) = dmatrix(i, j);
}

template<class T>
void identity(GeneralDenseMatrix<T> &dmatrix, size_t nrows, size_t ncols)
{
    dmatrix.setSize(nrows, ncols);

    for (size_t i = 0; i < nrows; i++)
        for (size_t j = 0; j < nrows; j++)
            dmatrix(i, j) = 1;
}

template<class T>
void zero(GeneralDenseMatrix<T> &dmatrix, size_t nrows, size_t ncols)
{
    dmatrix.setSize(nrows, ncols);

    for (size_t i = 0; i < nrows; i++)
        for (size_t j = 0; j < nrows; j++)
            dmatrix(i, j) = 0;
}

template<class T>
void random_matrix(GeneralSparseMatrix<T> &smatrix, size_t nrows, size_t ncols,
                   double fill_ratio, T minval, T maxval)
{
}

template<>
inline void random_matrix(GeneralSparseMatrix<double> &smatrix, size_t nrows,
                          size_t ncols, double fill_ratio, double minval, double maxval)
{
    smatrix.setSize(nrows, ncols);

    size_t colid, nfills = ncols*fill_ratio, zero = 0;

    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < nfills; j++) {
            colid = Math::random_value(zero, ncols);
            smatrix(i, colid) = Math::random_value(minval, maxval);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

template<class T1, class T2>
inline void apply_boundary_conditions(GeneralSparseMatrix<T1> &matrix,
                                      const std::map<size_t, T2> &boundmap, vector<T2> &rhs)
{
    assert(!boundmap.empty());

    size_t numRows = matrix.getSize1();

    rhs.resize(numRows);
    for (int i = 0; i < numRows; i++) rhs[i] = (T2) 0;

    typename GeneralSparseMatrix<T1>::SparseRow sparserow;
    typename GeneralSparseMatrix<T1>::const_row_iterator it;

    typename std::map<size_t, T2>::const_iterator biter;

    T1 coeff;

    size_t i, j;
    for (i = 0; i < numRows; i++) {
        sparserow = matrix.getConstRow(i);
        for (it = sparserow.begin(); it != sparserow.end(); ++it) {
            j = it->first;
            coeff = it->second;
            biter = boundmap.find(j);
            if (biter != boundmap.end()) {
                if (i != j) {
                    rhs[i] -= coeff * (biter->second);
                    matrix.setValue(i, j, 0);
                }
            }
        }

        biter = boundmap.find(i);
        if (biter != boundmap.end()) {
            matrix.clearRow(i);
            matrix.setValue(i, i, 1);
            rhs[i] = biter->second;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

template<class T1, class T2>
inline void apply_boundary_conditions(SymmetricSparseMatrix<T1> &matrix,
                                      std::map<size_t, T2> &boundmap, vector<T2> &rhs)
{
    assert(!boundmap.empty());

    GeneralSparseMatrix<T1> genmatrix = general_sparse_matrix(matrix, 1);
    apply_boundary_conditions(genmatrix, boundmap, rhs);
    assert(Math::isSymmetric(genmatrix));
    matrix = symmetric_sparse_matrix(genmatrix, 1);
}

/////////////////////////////////////////////////////////////////////////////

#define FULL_SYMMETRIC_SPARSE_MATRIX 1

template <class T>
class CSR_Matrix
{
public:

    int apply(SymmetricSparseMatrix<T> &matrix, bool fullmatrix = 0, bool deleteOrg = 0) {
        deleteOriginal = deleteOrg;

        if (fullmatrix)
            sym_full_matrix(matrix);
        else
            sym_off_matrix(matrix);

        return 0;
    }

    void apply(GeneralSparseMatrix<T> &matrix, bool deleteOrg = 0) {
        deleteOriginal = deleteOrg;
        gen_matrix(matrix);
    }

    vector<T> data;
    vector<size_t> colindex;
    vector<size_t> offset;
private:
    bool deleteOriginal;
    void sym_full_matrix(SymmetricSparseMatrix<T> &m);
    void sym_off_matrix(SymmetricSparseMatrix<T> &m);
    void gen_matrix(GeneralSparseMatrix<T> &m);
};

/////////////////////////////////////////////////////////////////////////////

template<class T>
void CSR_Matrix<T>::gen_matrix(GeneralSparseMatrix<T> &matrix)
{
    size_t numCols, numRows = matrix.getSize1();

    typename SparseMatrixRep<T>::SparseRow arow;
    typename SparseMatrixRep<T>::const_row_iterator jt;

    data.clear();
    colindex.clear();

    size_t nnz = matrix.getNumNonZeros();

    data.reserve(nnz);
    colindex.reserve(nnz);
    offset.resize(numRows + 1);

    offset[0] = 0;

    for (size_t i = 0; i < numRows; i++) {
        arow = matrix.getConstRow(i);
        numCols = 0;
        for (jt = arow.begin(); jt != arow.end(); ++jt) {
            colindex.push_back(jt->first);
            data.push_back(jt->second);
            numCols++;
        }
        if (deleteOriginal) matrix.clearRow(i);
        offset[i + 1] = offset[i] + numCols;
    }

    assert(data.size() == nnz);
};

template<class T>
void CSR_Matrix<T> ::sym_full_matrix(SymmetricSparseMatrix<T> &smatrix)
{
    GeneralSparseMatrix<T> genmat = general_sparse_matrix(smatrix);
    apply(genmat, 1);
};

template<class T>
void CSR_Matrix<T> ::sym_off_matrix(SymmetricSparseMatrix<T> &matrix)
{
    size_t numCols, numRows = matrix.getSize1();

    typename SparseMatrixRep<T>::SparseRow arow;
    typename SparseMatrixRep<T>::const_row_iterator jt;

    data.clear();
    colindex.clear();

    size_t nnz = matrix.getNumNonZeros();

    data.reserve(nnz);
    colindex.reserve(nnz);
    offset.resize(numRows + 1);

    offset[0] = 0;

    for (size_t i = 0; i < numRows; i++) {
        arow = matrix.getConstRow(i);
        numCols = 0;
        for (jt = arow.begin(); jt != arow.end(); ++jt) {
            colindex.push_back(jt->first);
            data.push_back(jt->second);
            numCols++;
        }
        if (deleteOriginal) matrix.clearRow(i);
        offset[i + 1] = offset[i] + numCols;
    }

    assert(data.size() == nnz);
};

#endif

