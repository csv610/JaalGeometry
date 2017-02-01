#include "CurveShorteningFlow.hpp"

/////////////////////////////////////////////////////////////////

void JCurveShorteningFlow :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    mesh->getTopology()->getBoundary(boundLoops);
    int nloops = boundLoops.size();
    boundNodes.resize(nloops);
    for( int i = 0; i < nloops; i++)
        JEdgeTopology::getChainNodes( boundLoops[i], boundNodes[i] );
}

/////////////////////////////////////////////////////////////////
double JCurveShorteningFlow :: getCurveLength() const
{
    int nloops = boundLoops.size();
    double sum = 0.0;
    for( int i = 0; i < nloops; i++)
        sum += JEdgeGeometry::getLength(boundLoops[i] );

    return sum;
}

/////////////////////////////////////////////////////////////////
void JCurveShorteningFlow :: performOneStep()
{
    int nloops = boundLoops.size();
    for( int i = 0; i < nloops; i++)
        JEdgeGeometry::makeUniform( boundLoops[i] );

    vector<Point3D> newPoints;
    Point3D pnew;
    for( int i = 0; i < nloops; i++) {
        int np = boundNodes[i].size();
        newPoints.resize(np);
        for( int j = 0; j < np; j++) {
            const Point3D &p0 = boundNodes[i][(j+0)%np]->getXYZCoords();
            const Point3D &p1 = boundNodes[i][(j+2)%np]->getXYZCoords();
            pnew[0] = 0.5*(p0[0] + p1[0]);
            pnew[1] = 0.5*(p0[1] + p1[1]);
            pnew[2] = 0.5*(p0[2] + p1[2]);
            newPoints[(j+1)%np] = pnew;
        }
        for( int j = 0; j < np; j++)
            boundNodes[i][j]->setXYZCoords( newPoints[j] );
    }
}
/////////////////////////////////////////////////////////////////
