#include "UniqueCoordinates.hpp"

///////////////////////////////////////////////////////////////////////////////////////

vector<Point3D> UniqueCoordinates:: getPoints()
{
    if( uniquePoints.empty() ) {
        useNaive();
    }
    return  uniquePoints;
}
///////////////////////////////////////////////////////////////////////////////////////

void UniqueCoordinates:: useNaive()
{
    size_t nPoints = pCloud.size();
    if (nPoints < 2) return;
    vector<bool> accept(nPoints);
    for( size_t i = 0; i < nPoints; i++)
        accept[i] = 1;

    ComparePoints cmp;
    cmp.eps = epsilon;
    for( int i = 0; i < nPoints; i++) {
        if( accept[i] ) {
            for( int j = i+1; j < nPoints; j++) {
                int dist = cmp( pCloud[i], pCloud[j] );
                if( dist == 0 ) {
                     cout << "Do not accept " << endl;
                     cout << i << " 1st Point " << pCloud[i][0] << " " << pCloud[i][1] << " " << pCloud[i][2] << endl;
                     cout << j << " 2nd Point " << pCloud[j][0] << " " << pCloud[j][1] << " " << pCloud[j][2] << endl;
                     accept[j] = 0;
                }
            }
        }
    }
    uniquePoints.reserve(nPoints);
    for( int i = 0; i < nPoints; i++)
        if( accept[i] ) uniquePoints.push_back( pCloud[i] );
}
///////////////////////////////////////////////////////////////////////////////////////
int UniqueCoordinates:: getID( const Point3D &queryPoint, size_t &id)
{
   if( uniquePoints.empty() )  getPoints();

   size_t nSize = uniquePoints.size();
   for(size_t i = 0; i < nSize; i++) {
       if( areEqual( uniquePoints[i], queryPoint) ) {
           id = i;
           return 0;
       }
   }
   return 1;
}
///////////////////////////////////////////////////////////////////////////////////////


