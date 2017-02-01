#pragma once

#include "basic_math.hpp"
#include <boost/smart_ptr.hpp>

#include <vector>
#include <ANN/ANN.h>
using namespace std;

class UniqueCoordinates
{
   struct ComparePoints {
        bool operator()(const Point3D &a, const Point3D &b) const
        {   
            if( (fabs(a[0] - b[0]) <= eps ) && 
                (fabs(a[1] - b[1]) <= eps ) && 
                (fabs(a[2] - b[2]) <= eps ))  return 0;
            return 1;
        }   
        double eps = 0.0;
    };

    bool areEqual( const Point3D &a, const Point3D &b)  
    {
          if( a[0] != b[0] ) return 0;
          if( a[1] != b[1] ) return 0;
          if( a[2] != b[2] ) return 0;
          return 1;
    }

public:
    static const int  USE_STD_SORTING = 0;
    static const int  USE_BOOST_SORTING = 1;
    static const int  USE_ANN     = 2;
    static const int  USE_NAIVE   = 3;

    void setPoints( const vector<Point3D> &p) {
         pCloud = p;
         uniquePoints.clear();
    }

    void  setAlgorithm(int a) { algorithm = a; 
          uniquePoints.clear();
    }

    void  setTolerance (double e ) {
        epsilon = e;
    }

    vector<Point3D> getPoints();
    int  getID( const Point3D &p, size_t &id);

    double  getMinDistance() const;
private:
    vector<Point3D>  pCloud;
    vector<Point3D>  uniquePoints;
    double  epsilon   = 0.0;
    int     algorithm = USE_STD_SORTING;

    ANNpointArray     dataPts;
    ANNpoint          queryPoint;
    ANNidxArray       nnIdx;
    ANNdistArray      dists;
    boost::scoped_ptr<ANNkd_tree> kdTree;

    void buildTree();
    void stdSort();
    void boostSort();
    void useANN();
    void useNaive();
};
