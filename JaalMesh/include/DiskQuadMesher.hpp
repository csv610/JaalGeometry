#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>

#include <vector>
#include <map>
#include <algorithm>

#include "tfiblend.hpp"

using namespace std;

struct DiskQuadMesher {

    void setBoundary(const vector<double> &p) {
        boundCoords = p;
    }

    int getQuadMesh( int nr, vector<Point2D> &points, vector<int> &quads) {
        nR = nr;
        return execute();
    }

    int saveAs( const string &s);

private:
//   Input:
    int   nR, nT;  // Number of nodes along radial and tangential directions.
    int   nBounds;
    vector<Point2D>  boundCoords;  // Coordinates along the circumfernce in CCW dir.

    struct Edge {
        int connect[2];
        bool isSame(int n0, int n1) const {
            if( connect[0] == n0 && connect[1] == n1 ) return  1;
            if( connect[0] == n1 && connect[1] == n0 ) return -1;
            return 0;
        }
        vector<int>  inserted_nodes;
    };

    int global_id;
    vector<Edge> edges;

    vector<Point2D>  quadCoords; // Coordinate of quadmesh. (Boundary included);
    vector<int>      quadConnect;// Connectivity of quadmesh.

    int  landmarks[17]; // There are 17 landmark nodes in the full circle.
    int  execute();
    void getCentroid();// Centroid of the boundary nodes.

    void mid_point( int n0, int n1, int &nw, double r = 0.0);
    void split_boundary( int n0, int n1);
    int  linear_interpolation( int n0, int n1, int np, int edgeid );
    int  quad_template( int n0, int n1, int n2, int n3 );
};

////////////////////////////////////////////////////////////////////////////////////

inline
void DiskQuadMesher :: mid_point( int n0, int n1, int &nw, double r)
{

    double xm = TFI::linear_interpolation(r, quadCoords[2*n0 + 0], quadCoords[2*n1+0] );
    double ym = TFI::linear_interpolation(r, quadCoords[2*n0 + 1], quadCoords[2*n1+1] );

    nw = quadCoords.size()/2;

    quadCoords.push_back( xm );
    quadCoords.push_back( ym );
}

////////////////////////////////////////////////////////////////////////////////////

inline
int DiskQuadMesher :: linear_interpolation( int n0, int n1, int np, int edgeid )
{
    assert( np >= 1 );
    vector<int>  edgenodes(np);

    double dt =  2.0/(double)(np-1);

    edgenodes[0] = landmarks[n0];
    for( int i = 1; i < np-1; i++) {
        double r = -1.0 + i*dt;
        mid_point( landmarks[n0], landmarks[n1], edgenodes[i], r);
    }
    edgenodes[np-1] = landmarks[n1];

    Edge edge;
    edge.connect[0] = n0;
    edge.connect[1] = n1;
    edge.inserted_nodes =  edgenodes;
    edges[edgeid] = edge;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline
void DiskQuadMesher::split_boundary(int n0, int n1)
{
    vector<int> edgenodes;
    assert( nT > 0);

    int dr = nT - 1;
    edgenodes.resize(nT);
    for(int i = 0; i < nT; i++) {
        edgenodes[i] = (n0*dr  + i)%nBounds;
    }

    Edge edge;
    edge.connect[0] = n0;
    edge.connect[1] = n1;
    edge.inserted_nodes =  edgenodes;
    edges[n0] = edge;
}

///////////////////////////////////////////////////////////////////////////////

inline
int DiskQuadMesher :: saveAs( const string &s)
{
    ofstream ofile( s.c_str(), ios::out);

    ofile << "OFF" << endl;

    int numNodes = quadCoords.size()/2;
    int numFaces = quadConnect.size()/4;

    ofile <<  numNodes << " " << numFaces << "  0 " << endl;

    for( int i = 0; i < numNodes; i++)
        ofile << quadCoords[2*i] << " " << quadCoords[2*i+1] <<  " 0.0 " << endl;

    for( int i = 0; i < numFaces; i++)
        ofile << " 4 " << quadConnect[4*i+0] << " "
              << quadConnect[4*i+1] << " "
              << quadConnect[4*i+2] << " "
              << quadConnect[4*i+3] << endl;
}

///////////////////////////////////////////////////////////////////////////////

inline
int DiskQuadMesher :: quad_template( int n0, int n1, int n2, int n3 )
{
    int err, vid;

    vector<int> ab = edges[n0].inserted_nodes;
    vector<int> bc = edges[n1].inserted_nodes;
    vector<int> dc = edges[n2].inserted_nodes;
    vector<int> ad = edges[n3].inserted_nodes;

    if( ab.size() != dc.size() ) return 2;
    if( ad.size() != bc.size() ) return 2;

    if( bc.back()  == ab.back()  ) std::reverse( bc.begin(), bc.end() );
    if( bc.back()  == dc.front() ) std::reverse( dc.begin(), dc.end() );
    if( dc.front() == ad.front() ) std::reverse( ad.begin(), ad.end() );

    int nx = ab.size();
    int ny = ad.size();

    assert( ab.front() == ad.front() );
    assert( ab.back()  == bc.front() );
    assert( bc.back()  == dc.back() );
    assert( ad.back()  == dc.front() );

    vector<double> x(nx*ny), y(nx*ny);
    vector<int> facenodes(nx*ny);

    for( int j = 0; j < ny; j++) {
        for( int i = 0; i < nx; i++) {
            vid = j*nx + i;
            x[vid] =  0.0;
            y[vid] =  0.0;
        }
    }

    for( int i = 0; i < nx; i++) {
        vid = i;
        facenodes[vid] =  ab[i];

        x[vid] =  quadCoords[2*ab[i] + 0];
        y[vid] =  quadCoords[2*ab[i] + 1];

        vid = (ny-1)*nx + i;
        facenodes[vid] =  dc[i];
        x[vid] =  quadCoords[2*dc[i] + 0];
        y[vid] =  quadCoords[2*dc[i] + 1];
    }

    for( int j = 1; j < ny-1; j++) {
        vid = j*nx;
        facenodes[vid] =  ad[j];
        x[vid] =  quadCoords[2*ad[j] + 0];
        y[vid] =  quadCoords[2*ad[j] + 1];

        vid = j*nx + (nx-1);
        facenodes[vid] =  bc[j];
        x[vid] =  quadCoords[2*bc[j] + 0];
        y[vid] =  quadCoords[2*bc[j] + 1];
    }

    TFI::blend_from_edges( &x[0], nx, ny);
    TFI::blend_from_edges( &y[0], nx, ny);

    global_id = quadCoords.size()/2;
    for( int j = 1; j < ny-1; j++) {
        for( int i = 1; i < nx-1; i++) {
            vid = j*nx + i;
            facenodes[vid] = global_id++;
            quadCoords.push_back( x[vid] );
            quadCoords.push_back( y[vid] );
        }
    }

    for( int j = 0; j < ny-1; j++) {
        for( int i = 0; i < nx-1; i++) {
            vid = j*nx + i;
            quadConnect.push_back( facenodes[vid] );

            vid = j*nx + i+1;
            quadConnect.push_back( facenodes[vid] );

            vid = (j+1)*nx + i+1;
            quadConnect.push_back( facenodes[vid] );

            vid = (j+1)*nx + i;
            quadConnect.push_back( facenodes[vid] );

        }
    }
}

///////////////////////////////////////////////////////////////////////////////
inline
int  DiskQuadMesher:: execute()
{
    int vid;
    quadCoords.clear();
    quadConnect.clear();

    edges.resize(28);

    nBounds = boundCoords.size()/2;

    if( nBounds  < 8 ) return 1;
    if( nBounds%8 )    return 2;

    int dr = nBounds/8;
    for( int i = 0; i < 8; i++) {
        vid = i*dr;
        landmarks[i] = vid;
    }
    nT = dr + 1;

    for( int i = 0; i < nBounds; i++) {
        quadCoords.push_back( boundCoords[2*i + 0] );
        quadCoords.push_back( boundCoords[2*i + 1] );
    }

    double centroid[] = { 0.0, 0.0};

    for( size_t i = 0; i < nBounds; i++)  {
        centroid[0] += boundCoords[2*i + 0];
        centroid[1] += boundCoords[2*i + 1];
    }

    centroid[0] /= ( double) nBounds;
    centroid[1] /= ( double) nBounds;

    landmarks[8] = nBounds;
    quadCoords.push_back( centroid[0] );
    quadCoords.push_back( centroid[1] );

    global_id = nBounds + 1;
    mid_point( landmarks[8], landmarks[0], landmarks[9]  );
    mid_point( landmarks[8], landmarks[1], landmarks[10] );
    mid_point( landmarks[8], landmarks[2], landmarks[11] );
    mid_point( landmarks[8], landmarks[3], landmarks[12] );
    mid_point( landmarks[8], landmarks[4], landmarks[13] );
    mid_point( landmarks[8], landmarks[5], landmarks[14] );
    mid_point( landmarks[8], landmarks[6], landmarks[15] );
    mid_point( landmarks[8], landmarks[7], landmarks[16] );

    for( int i = 0; i < 8; i++)
        split_boundary( i, (i+1)%8 );

    assert( nT >= 1 );

    linear_interpolation( 9,  10, nT, 16 );
    linear_interpolation( 10, 11, nT, 17 );
    linear_interpolation( 11, 12, nT, 18 );
    linear_interpolation( 12, 13, nT, 19 );
    linear_interpolation( 13, 14, nT, 20 );
    linear_interpolation( 14, 15, nT, 21 );
    linear_interpolation( 15, 16, nT, 22 );
    linear_interpolation( 16,  9, nT, 23 );

    assert( nR > 1 );
    linear_interpolation( 9,   0, nR,  8 );
    linear_interpolation( 10,  1, nR,  9 );
    linear_interpolation( 11,  2, nR, 10 );
    linear_interpolation( 12,  3, nR, 11 );
    linear_interpolation( 13,  4, nR, 12 );
    linear_interpolation( 14,  5, nR, 13 );
    linear_interpolation( 15,  6, nR, 14 );
    linear_interpolation( 16,  7, nR, 15 );

    linear_interpolation( 8,  9,  nT, 24 );
    linear_interpolation( 8,  11, nT, 25 );
    linear_interpolation( 8,  13, nT, 26 );
    linear_interpolation( 8,  15, nT, 27 );

    // First quarter circle
    quad_template(0, 9, 16, 8  );
    quad_template(1, 10, 17, 9 );
    quad_template(16, 17, 25, 24);

    // Second quarter circle
    quad_template(2, 11, 18, 10 );
    quad_template(3, 12, 19, 11 );
    quad_template(18, 19, 26, 25 );

    // Third quarter circle
    quad_template(4, 13, 20, 12 );
    quad_template(5, 14, 21, 13 );
    quad_template(20, 21, 27, 26 );

    // Forth quarter circle
    quad_template(6, 15, 22, 14 );
    quad_template(7, 8,  23, 15 );
    quad_template(22, 23, 24, 27);
}
///////////////////////////////////////////////////////////////////////////////

#ifdef DEF_MAIN
int main()
{
    vector<double> points;

    int nr = 8*3;
    double dt = 2.0*M_PI/(double)nr;

    points.resize( 2*nr );
    for( int i = 0; i < nr; i++) {
        points[2*i]   = cos( i*dt );
        points[2*i+1] = sin( i*dt );
    }

    DiskQuadMesher dqm;
    dqm.setBoundary( points );

    vector<double>  qPoints;
    vector<int>     quads;

    dqm.getQuadMesh(5, qPoints, quads);

    dqm.saveAs("tmp.off");
}
#endif
