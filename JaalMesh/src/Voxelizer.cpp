#include "MeshVoxelizer.hpp"

/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(double boxcenter[3],      */
/*          double boxhalfsize[3],double triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/
#include <math.h>
#include <stdio.h>

#define X 0
#define Y 1
#define Z 2

#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2];

#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;

static int planeBoxOverlap(double normal[3],double d, double maxbox[3])
{
    int q;
    double vmin[3],vmax[3];
    for(q=X; q<=Z; q++)
    {
        if(normal[q]>0.0f)
        {
            vmin[q]=-maxbox[q];
            vmax[q]=maxbox[q];
        }
        else
        {
            vmin[q]=maxbox[q];
            vmax[q]=-maxbox[q];
        }
    }
    if(DOT(normal,vmin)+d>0.0f) return 0;
    if(DOT(normal,vmax)+d>=0.0f) return 1;

    return 0;
}


/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)             \
    p0 = a*v0[Y] - b*v0[Z];                    \
    p2 = a*v2[Y] - b*v2[Z];                    \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fa, fb)              \
    p0 = a*v0[Y] - b*v0[Z];                    \
    p1 = a*v1[Y] - b*v1[Z];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
    if(min>rad || max<-rad) return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)             \
    p0 = -a*v0[X] + b*v0[Z];                   \
    p2 = -a*v2[X] + b*v2[Z];                       \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fa, fb)              \
    p0 = -a*v0[X] + b*v0[Z];                   \
    p1 = -a*v1[X] + b*v1[Z];                       \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
    if(min>rad || max<-rad) return 0;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)             \
    p1 = a*v1[X] - b*v1[Y];                    \
    p2 = a*v2[X] - b*v2[Y];                    \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fa, fb)              \
    p0 = a*v0[X] - b*v0[Y];                \
    p1 = a*v1[X] - b*v1[Y];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
    if(min>rad || max<-rad) return 0;

static int triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3])
{

    /*    use separating axis theorem to test overlap between triangle and box */
    /*    need to test for overlap in these directions: */
    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
    /*       we do not even need to test these) */
    /*    2) normal of the triangle */
    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
    /*       this gives 3x3=9 more tests */
    double v0[3],v1[3],v2[3];
    double min,max,d,p0,p1,p2,rad,fex,fey,fez;
    double normal[3],e0[3],e1[3],e2[3];

    /* This is the fastest branch on Sun */
    /* move everything so that the boxcenter is in (0,0,0) */
    SUB(v0,triverts[0],boxcenter);
    SUB(v1,triverts[1],boxcenter);
    SUB(v2,triverts[2],boxcenter);

    /* compute triangle edges */
    SUB(e0,v1,v0);      /* tri edge 0 */
    SUB(e1,v2,v1);      /* tri edge 1 */
    SUB(e2,v0,v2);      /* tri edge 2 */

    /* Bullet 3:  */
    /*  test the 9 tests first (this was faster) */
    fex = fabs(e0[X]);
    fey = fabs(e0[Y]);
    fez = fabs(e0[Z]);
    AXISTEST_X01(e0[Z], e0[Y], fez, fey);
    AXISTEST_Y02(e0[Z], e0[X], fez, fex);
    AXISTEST_Z12(e0[Y], e0[X], fey, fex);

    fex = fabs(e1[X]);
    fey = fabs(e1[Y]);
    fez = fabs(e1[Z]);
    AXISTEST_X01(e1[Z], e1[Y], fez, fey);
    AXISTEST_Y02(e1[Z], e1[X], fez, fex);
    AXISTEST_Z0(e1[Y], e1[X], fey, fex);

    fex = fabs(e2[X]);
    fey = fabs(e2[Y]);
    fez = fabs(e2[Z]);
    AXISTEST_X2(e2[Z], e2[Y], fez, fey);
    AXISTEST_Y1(e2[Z], e2[X], fez, fex);
    AXISTEST_Z12(e2[Y], e2[X], fey, fex);

    /* Bullet 1: */
    /*  first test overlap in the {x,y,z}-directions */
    /*  find min, max of the triangle each direction, and test for overlap in */
    /*  that direction -- this is equivalent to testing a minimal AABB around */
    /*  the triangle against the AABB */

    /* test in X-direction */
    FINDMINMAX(v0[X],v1[X],v2[X],min,max);
    if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;

    /* test in Y-direction */
    FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
    if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;

    /* test in Z-direction */
    FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
    if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;

    /* Bullet 2: */
    /*  test if the box intersects the plane of the triangle */
    /*  compute plane equation of triangle: normal*x+d=0 */
    CROSS(normal,e0,e1);
    d=-DOT(normal,v0);  /* plane eq: normal.x+d=0 */
    if(!planeBoxOverlap(normal,d,boxhalfsize)) return 0;

    return 1;   /* box and triangle overlaps */
}

/////////////////////////////////////////////////////////////////////////////////////////////

JVoxelMeshPtr JMeshVoxelizer :: genVoxels( const JMeshPtr &mesh, int nsize)
{
    if( mesh == nullptr) return nullptr;

    JBoundingBox box;
    box = mesh->getGeometry()->getBoundingBox();

    length[0] = 1.1*box.getLength(0);
    length[1] = 1.1*box.getLength(1);
    length[2] = 1.1*box.getLength(2);

    double minlen = *min_element( length, length + 3);
    double maxlen = *max_element( length, length + 3);
    double avglen = 0.5*(minlen + maxlen);

    double dx = avglen/(double)nsize;

    cellDim[0]  = max(1, (int) (length[0]/dx));
    cellDim[1]  = max(1, (int) (length[1]/dx));
    cellDim[2]  = max(1, (int) (length[2]/dx));

    nodeDim[0]  = cellDim[0] + 1;
    nodeDim[1]  = cellDim[1] + 1;
    nodeDim[2]  = cellDim[2] + 1;

    Point3D  center = box.getCenter();

    origin[0] = center[0] - 0.5*length[0];
    origin[1] = center[1] - 0.5*length[1];
    origin[2] = center[2] - 0.5*length[2];

    gridSpacing[0] = length[0]/(double)cellDim[0];
    gridSpacing[1] = length[1]/(double)cellDim[1];
    gridSpacing[2] = length[2]/(double)cellDim[2];

    voxmesh.reset( new JVoxelMesh );
    bgmesh = voxmesh->getBackgroundMesh(cellDim, length, origin);

    if( bgmesh == nullptr) {
        cout << "Warning: Voxel background mesh not created " << endl;
        return nullptr;
    }

    size_t numCells = bgmesh->getSize(3);

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = bgmesh->getCellAt(i);
        cell->setStatus(JMeshEntity::INACTIVE);
    }

    size_t numfaces = mesh->getSize(2);

//#pragma omp for
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        checkOverlap(face);
    }

    return voxmesh;
}
/////////////////////////////////////////////////////////////////////////////////////////////
void JMeshVoxelizer :: checkOverlap( const JFacePtr &face)
{
    if( !face->isActive() )  return;

    double triVerts[3][3];

    for( int i = 0; i < 3; i++) {
        const JNodePtr &vtx = face->getNodeAt(i);
        const Point3D &xyz  = vtx->getXYZCoords();
        triVerts[i][0]  = xyz[0];
        triVerts[i][1]  = xyz[1];
        triVerts[i][2]  = xyz[2];
    }

    double xmin = triVerts[0][0];
    double ymin = triVerts[0][1];
    double zmin = triVerts[0][2];

    double xmax = xmin;
    double ymax = ymin;
    double zmax = zmin;

    for( int i = 0; i < 3; i++) {
        xmin = min( xmin, triVerts[i][0] );
        ymin = min( ymin, triVerts[i][1] );
        zmin = min( zmin, triVerts[i][2] );

        xmax = max( xmax, triVerts[i][0] );
        ymax = max( ymax, triVerts[i][1] );
        zmax = max( zmax, triVerts[i][2] );
    }

    int i0 = (xmin-origin[0])/gridSpacing[0] -1.0;
    int j0 = (ymin-origin[1])/gridSpacing[1] -1.0;
    int k0 = (zmin-origin[2])/gridSpacing[2] -1.0;

    int i1 = (xmax-origin[0])/gridSpacing[0]+1.0;
    int j1 = (ymax-origin[1])/gridSpacing[1]+1.0;
    int k1 = (zmax-origin[2])/gridSpacing[2]+1.0;

    for( int k = k0; k <= k1; k++) {
        for( int j = j0; j <= j1; j++) {
            for( int i = i0; i <= i1; i++) {
                checkOverlap(triVerts, i,j,k);
            }
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////

void JMeshVoxelizer :: checkOverlap( double triVerts[3][3], int i, int j, int k)
{
    double boxCenter[3];

    boxCenter[0] = origin[0] + (i+0.5)*gridSpacing[0];
    boxCenter[1] = origin[1] + (j+0.5)*gridSpacing[1];
    boxCenter[2] = origin[2] + (k+0.5)*gridSpacing[2];

    if( i < 0 || i > cellDim[0]-1) return;
    if( j < 0 || j > cellDim[1]-1) return;
    if( k < 0 || k > cellDim[2]-1) return;

    if( triBoxOverlap(boxCenter, gridSpacing, triVerts) ) {
        size_t vid = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
        const JCellPtr &cell = bgmesh->getCellAt(vid);
        cell->setStatus(JMeshEntity::ACTIVE);
    }
}

//////////////////////////////////////////////////////////////////////////////
JVoxelMeshPtr JMeshVoxelizer :: genVoxels( const JNodeSequence &nodes, int nsize)
{
    if( nodes.empty() ) return nullptr;

    JBoundingBox box = JNodeGeometry::getBoundingBox(nodes);

    length[0] = 1.5*box.getLength(0);
    length[1] = 1.5*box.getLength(1);
    length[2] = 1.5*box.getLength(2);

    double minlen = *min_element( length, length + 3);
    double maxlen = *max_element( length, length + 3);
    double avglen = 0.5*(minlen + maxlen);

    double dx = avglen/(double)nsize;

    cellDim[0]  = length[0]/dx;
    cellDim[1]  = length[1]/dx;
    cellDim[2]  = length[2]/dx;

    nodeDim[0]  = cellDim[0] + 1;
    nodeDim[1]  = cellDim[1] + 1;
    nodeDim[2]  = cellDim[2] + 1;

    Point3D  center = box.getCenter();

    origin[0] = center[0] - 0.5*length[0];
    origin[1] = center[1] - 0.5*length[1];
    origin[2] = center[2] - 0.5*length[2];

    gridSpacing[0] = length[0]/(double)cellDim[0];
    gridSpacing[1] = length[1]/(double)cellDim[1];
    gridSpacing[2] = length[2]/(double)cellDim[2];

    voxmesh.reset( new JVoxelMesh );
    bgmesh = voxmesh->getBackgroundMesh(cellDim, length, origin);

    size_t numCells = bgmesh->getSize(3);
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = bgmesh->getCellAt(i);
        cell->setStatus(JMeshEntity::INACTIVE);
    }

    size_t nSize = nodes.size();
    for( size_t inode = 0; inode < nSize; inode++) {
        const Point3D &xyz = nodes[inode]->getXYZCoords();
        int i = (xyz[0]-origin[0])/gridSpacing[0];
        int j = (xyz[1]-origin[1])/gridSpacing[1];
        int k = (xyz[2]-origin[2])/gridSpacing[2];
        size_t vid = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
        const JCellPtr &cell = bgmesh->getCellAt(vid);
        cell->setStatus(JMeshEntity::ACTIVE);
    }
    return voxmesh;
}
///////////////////////////////////////////////////////////////////////////////
