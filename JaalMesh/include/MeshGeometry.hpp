#pragma once

#include <iomanip>
#include "basic_math.hpp"
#include <iterator>

#include "Mesh.hpp"
#include "MinDiam.hpp"
#include "Miniball.hpp"
#include "MeshQuality.hpp"
#include "BasicShapes.hpp"

namespace Jaal
{
class JMeshGeometry 
{
    static JLogger* logger;
public:
    static const int  NORMAL_NO_WEIGHT   = 0;
    static const int  NORMAL_AREA_WEIGHT = 1;

    static double getVolume( const JCellSequence &fs);
    static JCellPtr getMinBox( const vector<Point3D> &pnts);
    static JCellPtr getMinBox( const JNodeSequence &nodes);

    static JFacePtr getMinRectangle( const vector<Point2D> &nodes);
    static JFacePtr getMinRectangle( const JNodeSequence &nodes);

    static vector<double> getEuclideanDistance( const JMeshPtr &a, const JMeshPtr &b);

    static double getMaxDistance( const JMeshPtr &m1, const JMeshPtr &m2);

    explicit JMeshGeometry( const JMeshPtr &m)
    {
        mesh = m;
    }

    int getDimension() const;

    void update();

    JMeshPtr explode(double alpha);

    JBoundingBox getBoundingBox() const;

    JCellPtr getMinBox();
    JSphere  getMinSphere();

    void  addNoise( double maxVal);

    bool  hasPlanarEmbedding() const;

    int   spatialSort();

    // Get the center of entire mesh ...
    Point3D getCenter() const;

    // Get the surface area of the mesh ...
    double getSurfaceArea() const;

    // Get the volume of the mesh ( 3D only ).
    double getVolume() const;

    // Calculate the normals at the faces. ( Assumming planer faces ).
    void   setFacesNormal();

    // Calculate the normal at each vertex. It requires FaceNormals..
    void   setNodesNormal(int weight = NORMAL_NO_WEIGHT);

    // Calculate the angular defect at each vertex. This can be used to
    // identify the corners in the mesh...
    size_t setAngleDefects(double angle = 5);

    // Detect sharp edges in the mesh...
    size_t setSharpEdges( double angle = 30); // in Degrees

    // How many concave faces in the mesh..
    size_t count_concave_faces() const;

    // How many invert faces.
    size_t getNumOfInvertedElements() const;
    void   getInvertedFaces(JFaceSequence &f);
    void   getInvertedCells(JCellSequence &f);

    // void surfaceSegmentation();

    // Backup the coordinates of nodes. ( Connnectivity is backed up by Topology module ..
    void backup();

    // Retrieve the nodes coordinates from the backup data...
    void fromBackup();

    // Just freeup the coordinates arrary used for backup...
    void clearBackup();

    // Get the nodes coordinates as a stream ( Used for interface with other
    // software. example Mesquite )
    int getCoordsArray( vector<double> &a, vector<size_t> &l2g);

    // Reset the coordinates of nodes, ( Generally it comes from optimization
    // modules. The array size must match in order to reflect changes.
    int setCoordsArray(const vector<double> &v, const vector<size_t> &l2g);

    void   setFeatureLength();

    void   setBoundaryConstraints( bool preserve );
    JNodePtr  getNearest(const Point3D &p);

    // Get all the neighbours within dist from the vtx
    void   getEuclideanNeigbours( const JNodePtr &vtx, double dist,  JNodeSequence &seq);
    void   getGeodesicNeigbours(  const JNodePtr &vtx, double dist,  JNodeSequence &seq);

    JNodeSequence getBoundaryCorners(double angle) const;
    JNodeSequence getConvexCorners(double angle) const;
    JNodeSequence getConcaveCorners(double angle) const;
    double  getMeanEdgeLength() const;

private:
    JMeshPtr mesh;

    void setNormal(const JFacePtr &f);
    void setNormal(const JNodePtr &f, int weight = NORMAL_NO_WEIGHT);
    void getEdgeLengths( vector<double> &l);
};
};

