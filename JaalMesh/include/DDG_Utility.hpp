#pragma once
#include <cstdlib>
//#include "Utility.h"
#include "DDG_Complex.hpp"

namespace DDG
{
inline double sqr( double x )
{
    return x*x;
}

inline double unitRand( void )
{
    const double rRandMax = 1. / (double) RAND_MAX;

    return rRandMax * (double) rand();
}

inline double seconds( int t0, int t1 )
{
    return (double)(t1-t0) / (double) CLOCKS_PER_SEC;
}
}

namespace DDGConstants
{
static DDG::Complex ii( 0., 1. );
}
