#include "GeomPredicates.hpp"

void JGeomPredicates::init()
{
    exactinit();
}

int JGeomPredicates :: getPointOrientation( const double *pa, const double *pb, const double *pc, int mode)
{
    //
    // Given a line segment (in 2D: xy plane) defined by two end points i.e.(pa,pb), this
    // module checks if the point "pc" is in clockwise or clockwise direction.
    //
    // This method uses "Jonathan Shewchuks geometric predicates..
    //
    double val = 0.0;

    if( mode == JGeomPredicates::FAST)
        val = orient2dfast(pa, pb, pc);
    else {
        val = orient2d(pa, pb, pc);
    }

    if( val > 0) return JGeomPredicates::ANTI_CLOCKWISE;
    if( val < 0) return JGeomPredicates::CLOCKWISE;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGeomPredicates :: getPointOrientation( const double *pa, const double *pb, const double *pc, const double *pd, int mode)
{
    //
    // Given three points on a plane ( pa,pb,pc) which are anti-clockwise oriented, this function
    // tests if the point "pd" is below the plane. Assume that (pa,pb,pc) on the screen, this checks
    // if the point is inside the monitor, away from your eyes.
    //
    double val = 0.0;

    if( mode == JGeomPredicates::FAST)
        val = orient3dfast( pa, pb, pc, pd);
    else
        val = orient3d( pa, pb, pc, pd);

    if( val > 0) return  JGeomPredicates::OUTSIDE;
    if( val < 0) return  JGeomPredicates::INSIDE;
    return 0;

}

///////////////////////////////////////////////////////////////////////////////

int JGeomPredicates :: inCircle(const double *pa, const double *pb, const double *pc, const double *pd, int mode)
{
    double val = 0.0;

    if( mode == JGeomPredicates::FAST)
        val = incirclefast( pa, pb, pc, pd);
    else
        val = incircle( pa, pb, pc, pd);

    if( val > 0) return  JGeomPredicates::OUTSIDE;
    if( val < 0) return  JGeomPredicates::INSIDE;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGeomPredicates :: inSphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe, int mode)
{
    double val;

    if( mode == JGeomPredicates::FAST)
        val = inspherefast(pa, pb, pc, pd, pe);
    else
        val = insphere( pa, pb, pc, pd, pe);

    if( val > 0) return  JGeomPredicates::OUTSIDE;
    if( val < 0) return  JGeomPredicates::INSIDE;
    return 0;

}
///////////////////////////////////////////////////////////////////////////////

