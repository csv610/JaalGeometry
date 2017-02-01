#pragma once

#include "Mesh.hpp"
#include <vector>
#include <verdict.h>

using namespace std;
using namespace Jaal;
///////////////////////////////////////////////////////////////////////////

class JStatisticalInfo
{
public:
    JStatisticalInfo()
    {
        numSamples = 0;
        sum        = 0.0;
        minVal     = 0.0;
        maxVal     = 0.0;
        avgVal     = 0.0;
        medianVal  = 0.0;
    }
    void   setData(vector<double> &v);

    double getSum() const {
        return sum;
    }
    double getMinimum() const {
        return minVal;
    }
    double getMaximum() const {
        return maxVal;
    }
    double getAverage() const {
        return avgVal;
    }
    double getMedian()  const {
        return medianVal;
    }

private:
    size_t numSamples;
    double sum, minVal, maxVal, avgVal, medianVal;
};

///////////////////////////////////////////////////////////////////////////

class JMeshQuality
{
    static JLogger *logger;
public:
    static const int AREA = 0;
    static const int ASPECT_RATIO = 1;
    static const int ASPECT_BETA  = 2;
    static const int ASPECT_GAMMA = 3;
    static const int CONDITION_NUMBER = 4;
    static const int DIAGONAL_RATIO = 5;
    static const int EDGE_LENGTH  = 6;
    static const int EDGE_WEIGHT = 7;
    static const int MIN_ANGLE = 8;
    static const int MAX_ANGLE = 9;
    static const int DISTORTION = 10;
    static const int JACOBIAN = 11;
    static const int ODDY   = 12;
    static const int SCALED_JACOBIAN = 13;
    static const int RELATIVE_SIZE_SQUARED = 14;
    static const int SHAPE = 15;
    static const int SHAPE_AND_SIZE = 16;
    static const int SHEAR   = 17;
    static const int SHEAR_AND_SIZE  = 18;
    static const int SKEW  = 19;
    static const int STRETCH  = 20;
    static const int TAPER  = 21;
    static const int VOLUME  = 22;
    static const int WARPAGE = 23;

    JMeshQuality() {}

    JMeshQuality( const JMeshPtr &m ) {
        mesh = m;
        attribname = "Quality";
    }

    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    // You can add the quality attribute to the entity, if there is
    // some need to make them persistent...

    void addAttribute(const string &s) {
        attribname = s;
    }

    JStatisticalInfo  getBoundaryEdgesLength() const;
    JStatisticalInfo  getBoundaryFacesArea() const;
    JStatisticalInfo  getEdgesLength() const;

    JStatisticalInfo  getFacesQuality( int q) const;
    JStatisticalInfo  getCellsQuality( int q) const;

    vector<double> getEdgesQuality(int what, int pos, int sorted = 0);
    vector<double> getFacesQuality(int what, int pos, int sorted = 0);
    vector<double> getCellsQuality(int what, int pos, int sorted = 0);

    size_t   getAcceptableEdges( int what, int pos);
    size_t   getAcceptableFaces( int what, int pos);
    size_t   getAcceptableCells( int what, int pos);

    // Edge Quality
    double getEdgeLength( const JEdgePtr &e);
    double getEdgeWeight( const JEdgePtr &e, int dir = 1);

    // Face Quality
    double getArea( const JFacePtr &f);
    double getAspectRatio( const JFacePtr &f);
    double getConditionNumber( const JFacePtr &f);
    double getDistortion( const JFacePtr &f);
    double getJacobian(  const JFacePtr &f);
    double getMinAngle(  const JFacePtr &f);
    double getMaxAngle(  const JFacePtr &f);
    double getOddy( const JFacePtr &f);
    double getRelativeSize2(  const JFacePtr &f);
    double getScaledJacobian(  const JFacePtr &f);

    double getShear( const JFacePtr &f);
    double getShearSize( const JFacePtr &f);
    double getShape( const JFacePtr &f);
    double getShapeSize( const JFacePtr &f);
    double getSkew( const JFacePtr &f);
    double getStretch( const JFacePtr &f);
    double getTaper( const JFacePtr &f);
    double getWarpage( const JFacePtr &f);

    // Face Quality
    double getAspectRatio( const JCellPtr &f);
    double getAspectBeta( const JCellPtr &f);
    double getAspectGamma( const JCellPtr &f);
    double getConditionNumber( const JCellPtr &f);
    double getDimension( const JCellPtr &f);
    double getDiagonalRatio( const JCellPtr &f);
    double getDistortion( const JCellPtr &f);
    double getOddy( const JCellPtr &f);
    double getJacobian( const JCellPtr &f);
    double getRelativeSize2( const JCellPtr &f);
    double getShape( const JCellPtr &f);
    double getShapeSize( const JCellPtr &f);
    double getScaledJacobian( const JCellPtr &f);
    double getShear( const JCellPtr &f);
    double getShearSize( const JCellPtr &f);
    double getSkew( const JCellPtr &f);
    double getStretch( const JCellPtr &f);
    double getTaper( const JCellPtr &f);
    double getVolume(const JCellPtr &c);

private:
    JMeshPtr mesh;
    Point3D xyz;
    double  coords[8][3];
    string  attribname;
    vector<double> quality;
    void registerQuality( const JFacePtr &f, double val, int pos);
    void registerQuality( const JEdgePtr &e, double val, int pos);
    vector<double> minQuadAcceptable, maxQuadAcceptable;
    vector<double> minTriAcceptable,  maxTriAcceptable;
    void setRange();
};

