#pragma once

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <values.h>
#include <iostream>
#include <vector>
using namespace std;

#define ANGLE_IN_DEGREES  0
#define ANGLE_IN_RADIANS  1

#include <boost/utility.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>
#include <boost/array.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef boost::array<int,2>    Point2I;
typedef boost::array<int,3>    Point3I;
typedef boost::array<float,2>  Point2F;
typedef boost::array<float,3>  Point3F;
typedef boost::array<float,4>  Point4F;

typedef boost::array<double,2> Point2D;
typedef boost::array<double,3> Point3D;
typedef boost::array<double,4> Point4D;

typedef boost::array<float,3> Array3F;
typedef boost::array<float,4> Array4F;

typedef boost::array<double,2> Array2D;
typedef boost::array<double,3> Array3D;
typedef boost::array<double,4> Array4D;

typedef boost::array<int,2> Array2I;
typedef boost::array<int,3> Array3I;
typedef boost::array<int,4> Array4I;

typedef boost::array<float,2>  Vec2F;
typedef boost::array<float,3>  Vec3F;

typedef boost::array<double,2> Vec2D;
typedef boost::array<double,3> Vec3D;
typedef boost::array<double,4> Vec4D;

#include "FastArea.hpp"

struct LinearSystem
{
    static int solve(Eigen::SparseMatrix<double> &K, vector<double> &b, vector<double> &x);
};

namespace JMath
{

template<class T>
inline T max_value( const T &a, const T &b, const T &c)
{
    return max(a,max(b, c));
}

inline double Radian2Degree( double r )
{
    return 180.0*r/M_PI;
}
inline double Degree2Radian( double r )
{
    return M_PI*r/180.0;
}
///////////////////////////////////////////////////////////////////////////////

template<class T>
inline T  high_frequency_item( const vector<T> &v)
{
  std::map<T,size_t> itemCount;
  for( size_t i = 0; i < v.size(); i++) {
       if( itemCount.find(v[i]) == itemCount.end() ) 
           itemCount[v[i]] = 1;
       else
           itemCount[v[i]]++;
  }

  size_t maxCount = 0;
  for( auto keyVal : itemCount) 
        if( keyVal.second > maxCount ) maxCount++;

  for( auto keyVal : itemCount) 
        if( keyVal.second == maxCount ) return keyVal.first;
}
///////////////////////////////////////////////////////////////////////////////

template<class T>
inline T  low_frequency_item( const vector<T> &v)
{

  std::map<T,size_t> itemCount;
  for( size_t i = 0; i < v.size(); i++) {
       if( itemCount.find(v[i]) == itemCount.end() ) 
           itemCount[v[i]] = 1;
       else
           itemCount[v[i]]++;
  }

  size_t minCount = v.size();
  for( auto keyVal : itemCount) 
        if( keyVal.second < minCount ) minCount++;

  for( auto keyVal : itemCount) 
        if( keyVal.second == minCount ) return keyVal.first;
}

///////////////////////////////////////////////////////////////////////////////

template<class T>
inline T min_value( const T &a, const T &b, const T &c)
{
    return min(a,min(b, c));
}

inline int iround( double r)
{
    return (r > 0.0) ? (r+0.5) : (r-0.5);
}

template<class T>
inline T min_value( const T &a, const T &b, const T &c, const T &d)
{
    return min(d, min(a,min(b, c)));
}

inline double length( const Point3D &A, const Point3D &B)
{
    double dx = A[0] - B[0];
    double dy = A[1] - B[1];
    double dz = A[2] - B[2];
    return sqrt( dx*dx + dy*dy + dz*dz );
}

inline double length2( const Point3D &A, const Point3D &B)
{
    double dx = A[0] - B[0];
    double dy = A[1] - B[1];
    double dz = A[2] - B[2];
    return dx*dx + dy*dy + dz*dz;
}

template<class T>
inline double magnitude( const boost::array<T,3> &A )
{
    return sqrt( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] );
}

template<class T>
inline double dot_product( const boost::array<T,3> &A, const boost::array<T,3> &B)
{
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

template<class T>
inline double dot_product( const boost::array<T,2> &A, const boost::array<T,2> &B)
{
    return A[0]*B[0] + A[1]*B[1];
}

template<class T>
inline void cross_product( const boost::array<T,3> &A, const boost::array<T,3> &B, boost::array<T,3> &C)
{
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
}

template<class T>
T random_value(T minVal, T maxVal)
{
    return minVal + drand48()*(maxVal - minVal);
}

inline void getTriAngles(const Point3D &pa, const Point3D &pb, const Point3D &pc, Point3D &angles, int measure =ANGLE_IN_DEGREES)
{
    angles[0] = 0.0;
    angles[1] = 0.0;
    angles[2] = 0.0;

    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );
    double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );
    double cosB =  (a2 + c2 - b2)/(2*sqrt(a2*c2) );
    double cosC =  (a2 + b2 - c2)/(2*sqrt(a2*b2) );

    if( cosA >  1.0) cosA =  1.0;
    if( cosA < -1.0) cosA = -1.0;
    angles[0] = acos(cosA);

    if( cosB >  1.0) cosB =  1.0;
    if( cosB < -1.0) cosB = -1.0;
    angles[1] = acos(cosB);

    if( cosC >  1.0) cosC =  1.0;
    if( cosC < -1.0) cosC = -1.0;
    angles[2] = acos(cosC);

    if( measure ==  ANGLE_IN_DEGREES) {
        angles[0] *= 180/M_PI;
        angles[1] *= 180/M_PI;
        angles[2] *= 180/M_PI;
    }
}

///////////////////////////////////////////////////////////////////////////////

inline double getTriAngle(const Point3D &pa, const Point3D &pb, const Point3D &pc,
                          int measure = ANGLE_IN_DEGREES)
{
    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );
    double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );

    if( cosA >  1.0) cosA =  1.0;
    if( cosA < -1.0) cosA = -1.0;

    double angle;
    if( measure == ANGLE_IN_DEGREES)
        angle = 180*acos(cosA)/M_PI;
    else
        angle = acos(cosA);
    return angle;
}

////////////////////////////////////////////////////////////////////////////////

inline
int getMaxTriAngle(const Point3D &pa, const Point3D &pb, const Point3D &pc, double &angle,
                   int measure = ANGLE_IN_DEGREES)
{
    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );

    double maxlen = max_value(a2,b2,c2);

    if( maxlen == a2) {
        double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );
        if( cosA >  1.0) cosA =  1.0;
        if( cosA < -1.0) cosA = -1.0;
        angle = acos(cosA);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        return 0;
    }

    if( maxlen == b2 ) {
        double cosB =  (a2 + c2 - b2)/(2*sqrt(a2*c2) );
        if( cosB >  1.0) cosB =  1.0;
        if( cosB < -1.0) cosB = -1.0;
        angle  = acos(cosB);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        return 1;
    }

    if( maxlen == c2 ) {
        double cosC =  (a2 + b2 - c2)/(2*sqrt(a2*b2) );
        if( cosC >  1.0) cosC =  1.0;
        if( cosC < -1.0) cosC = -1.0;
        angle = acos(cosC);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        return 2;
    }

    return -1;
}

////////////////////////////////////////////////////////////////////////////////

inline int getMinTriAngle(const Point3D &pa, const Point3D &pb, const Point3D &pc, double &angle,
                          int measure = ANGLE_IN_DEGREES)
{
    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );

    double minlen = min_value(a2,b2,c2);

    if( minlen == a2) {
        double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );
        if( cosA >  1.0) cosA =  1.0;
        if( cosA < -1.0) cosA = -1.0;
        angle = acos(cosA);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        return 0;
    }

    if( minlen == b2 ) {
        double cosB =  (a2 + c2 - b2)/(2*sqrt(a2*c2) );
        if( cosB >  1.0) cosB =  1.0;
        if( cosB < -1.0) cosB = -1.0;
        angle  = acos(cosB);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        return 1;
    }

    if( minlen == c2 ) {
        double cosC =  (a2 + b2 - c2)/(2*sqrt(a2*b2) );
        if( cosC >  1.0) cosC =  1.0;
        if( cosC < -1.0) cosC = -1.0;
        angle = acos(cosC);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        return 2;
    }

    return -1;
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
T mean_value( const vector<T> &v)
{
    assert( !v.empty() );
    vector<T> tmp(v);
    sort( tmp.begin(), tmp.end() );
    return tmp[v.size()/2];
}


template<class T>
T average_value( const vector<T> &v)
{
    size_t nsize = v.size();

    T sum = 0.0;
    for( size_t i = 0; i < nsize; i++)
        sum += v[i];
    T avg = sum/(double)nsize;
    return avg;
}

template<class T>
T standard_deviation( const vector<T> &v)
{
    assert( !v.empty() );

    size_t nsize = v.size();

    T avg = average_value( v );

    T sum = 0.0;
    for( size_t i = 0; i < nsize; i++) {
        double vm = v[i] - avg;
        sum += vm*vm;
    }
    return sqrt(1.0/(double)(nsize-1)*sum );
}

inline void make_vector( const Point3D &head, const Point3D &tail, Vec3D &xyz)
{
    xyz[0] = head[0] - tail[0];
    xyz[1] = head[1] - tail[1];
    xyz[2] = head[2] - tail[2];
}


inline void interpolate( const Point3D &p0, const Point3D &p1, Point3D &pmid, double r = 0.5)
{
    assert( r >= 0.0 && r <= 1.0);
    pmid[0] = (1-r)*p0[0] + r*p1[0];
    pmid[1] = (1-r)*p0[1] + r*p1[1];
    pmid[2] = (1-r)*p0[2] + r*p1[2];
}


////////////////////////////////////////////////////////////////////////////////

/*
inline double polyArea(const vector<Point3D> &p )
{
    // Formula from wikipedia....
    int nSize = p.size();
    assert( nSize >= 3);

    double sum = 0.0;
    for( int i  = 0; i < nSize; i++) {
        double x0 = p[i][0];
        double y0 = p[i][1];
        double x1 = p[(i+1)%nSize][0];
        double y1 = p[(i+1)%nSize][1];
        sum += x0*y1 - x1*y0;
        assert(p[i][2] == 0.0);
    }
    return 0.5*sum;
}

////////////////////////////////////////////////////////////////////////////////

inline double polyArea(const vector<Point2D> &p )
{
    // Formula from wikipedia....
    double sum = 0.0;
    int nSize = p.size();

    assert( nSize >= 3);

    for( int i  = 0; i < nSize; i++) {
        double x0 = p[i][0];
        double y0 = p[i][1];
        double x1 = p[(i+1)%nSize][0];
        double y1 = p[(i+1)%nSize][1];
        sum += x0*y1 - x1*y0;
    }
    return 0.5*sum;
}
////////////////////////////////////////////////////////////////////////////////

inline void polyCentroid( const vector<Point2D> &p, Point2D &c )
{
    int nSize = p.size();
    assert( nSize >= 3);

    // Formula from wikipedia. Note that centroid does not depend on the
    // orientation of the polygon. The division by Area will take care
    // of correct value.

    // For convex bodies, centroid is always inside the region.

    double cx = 0.0;
    double cy = 0.0;
    double cf;

    for( int i  = 0; i < nSize; i++) {
        cf  = p[i][0]*p[(i+1)%nSize][1] - p[(i+1)%nSize][0]*p[i][1];
        cx +=  (p[i][0]+p[(i+1)%nSize][0])*cf;
        cy +=  (p[i][1]+p[(i+1)%nSize][1])*cf;
    }

    double A = polyArea( p );

    c[0] = cx/(6.0*A);
    c[1] = cy/(6.0*A);
}
*/

////////////////////////////////////////////////////////////////////////////////

struct DimensionSorter {
    DimensionSorter( int d) {
        dim = d;
    }
    bool operator() ( const Point3D &pa, const Point3D &pb) {
        return pa[dim] < pb[dim];
    }
private:
    int dim;
};


struct SpatialSort {
    bool operator() ( const Point3D &pa, const Point3D &pb) {
        if( pa[2] < pb[2] ) return 1;
        if( pa[1] < pb[1] ) return 1;
        if( pa[0] < pb[0] ) return 1;
        return 0;
    }
};


////////////////////////////////////////////////////////////////////////////////

inline void normal( const Point3D &p0, const Point3D &p1, const Point3D &p2,
                    Vec3D &normal)
{
    Vec3D p1p0, p2p0;
    JMath::make_vector( p2, p0, p2p0);
    JMath::make_vector( p1, p0, p1p0);
    JMath::cross_product( p2p0, p1p0, normal);

    double mag = JMath::magnitude( normal );
    normal[0] /= mag;
    normal[1] /= mag;
    normal[2] /= mag;
}

///////////////////////////////////////////////////////////////////////////////
inline Vec3D unit_vector( const Vec3D &vec)
{
    double dl  = magnitude(vec);
    Vec3D  uvec;
    uvec[0] = vec[0]/dl;
    uvec[1] = vec[1]/dl;
    uvec[2] = vec[2]/dl;
    return uvec;
}
///////////////////////////////////////////////////////////////////////////////

inline int unit_vector( const Point3D &head, const Point3D &tail, Vec3D &uvec)
{
    make_vector( head, tail, uvec);
    double dl  = magnitude(uvec);

    uvec[0]  /=  dl;
    uvec[1]  /=  dl;
    uvec[2]  /=  dl;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
template<class T>
inline T Angle2D( T x1, T y1, T x2, T y2)
{
    double theta1 = atan2( (double)y1, (double)x1);
    double theta2 = atan2( (double)y2, (double)x2);

    double dtheta = theta2-theta1;

    if( dtheta >  M_PI) dtheta -= 2.0*M_PI;
    if( dtheta < -M_PI) dtheta += 2.0*M_PI;

    return dtheta;
}
///////////////////////////////////////////////////////////////////////////////

template<class T>
inline double getVecAngle( const boost::array<T,3> &A, const boost::array<T,3> &B, int measure)
{
    double AB = dot_product(A,B);
    double Am = magnitude(A);
    double Bm = magnitude(B);

    if( Am < 1.0E-15 || Bm < 1.0E-15) return 0.0;

    double x = AB/(Am*Bm);

    if( x > 1.0) x = 1.0;
    if( x < -1.0) x = -1.0;

    if( measure == ANGLE_IN_DEGREES ) return 180*acos(x)/M_PI;
    return acos(x);
}

} //End of namespace Math

#ifdef CSV
#ifndef HAVE_BOOST
template <class T, class n>
inline double getAngle(const Array<T, n> &VecA, const Array<T, n> &VecB,
                       int unit_measure)
{
    double Abar, Bbar, theta;
    Abar = magnitude(VecA);
    Bbar = magnitude(VecB);

    if (Abar < 1.0E-10 || Bbar < 1.0E-10) {
        cout << " Warning: Error in Angle calculation " << endl;
        cout << " Magnitude of Vector A is " << Abar << endl;
        cout << " Magnitude of Vector B is " << Bbar << endl;
        return 0.0;
    }

    double value = dot_product(VecA, VecB) / (Abar * Bbar);

    if (value > +1.0) value = +1.0;
    if (value < -1.0) value = -1.0;

    theta = acos(value);

    if (unit_measure == ANGLE_IN_DEGREES) theta *= (180.0 / M_PI);

    return theta;
}

//////////////////////////////////////////////////////////////////////////////

template <class T, size_t n>
inline T getAngle(const Array<T, n> &pa, const Array<T, n> &pb,
                  const Array<T, n> &pc, int unit_measure = ANGLE_IN_DEGREES)
{
    Array<T, n> VecA = make_vector(pb, pa);
    Array<T, n> VecB = make_vector(pc, pa);
    return getAngle(VecA, VecB, unit_measure);
}

#endif

///////////////////////////////////////////////////////////////////////////////
#endif

