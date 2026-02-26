#pragma once

// static void exactinit();
extern void exactinit();
extern double orient2d(const double *pa, const double *pb, const double *pc);
extern double orient2dfast(const double *pa, const double *pb, const double *pc);
extern double orient2dexact(const double *pa, const double *pb, const double *pc);

extern double orient3d(const double *pa, const double *pb, const double *pc, const double *pd);
extern double orient3dfast(const double *pa, const double *pb, const double *pc, const double *pd);
extern double orient3dexact(const double *pa, const double *pb, const double *pc, const double *pd);

extern double incircle(const double *pa, const double *pb, const double *pc, const double *pd);
extern double incirclefast(const double *pa, const double *pb, const double *pc, const double *pd);
extern double insphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe);
extern double inspherefast(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe);


class JGeomPredicates
{
    static void init();
public:

    static const int ADAPTIVE = 0;
    static const int INFINITE = 1;
    static const int FAST     = 2;

    static const int ANTI_CLOCKWISE  =  1;
    static const int CLOCKWISE       = -1;

    static const int INSIDE   = -1;
    static const int BOUNDARY =  0;
    static const int OUTSIDE  =  1;

    static int  getPointOrientation( const double *pa, const double *pb, const double *pc, int m = 0);
    static int  getPointOrientation( const double *pa, const double *pb, const double *pc, const double *pd, int m = 0);
    static int  inCircle(const double *pa, const double *pb, const double *pc, const double *pd, int m = 0);
    static int  inSphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe, int m = 0);
};

