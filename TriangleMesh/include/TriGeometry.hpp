#pragma once

struct TriGeometry
{

    static void getAngles(const Point3D &pa, const Point3D &pb, const Point3D &pc,
                          Point3D &angles, int measure =ANGLE_IN_DEGREES);
    static double getAngle(const Point3D &pa, const Point3D &pb, const Point3D &pc,
                           int measure = ANGLE_IN_DEGREES);

};


inline
void TriGeometry:: getAngles(const Point3D &pa, const Point3D &pb, const Point3D &pc, Point3D &angles, int measure =ANGLE_IN_DEGREES)
{
    angles[0] = 0.0;
    angles[1] = 0.0;
    angles[2] = 0.0;

    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );
    double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );
    double cosB =  (a2 + c2 - b2)/(2*sqrt(a2*c2) );
    double cosC =  (a2 + b2 - c2)/(2*sqrt(a2*b2) );

    if( cosA >  1.0) cosA =  1.0;
    if( cosA < -1.0) cosA = -1.0;
    angles[0] = acos(cosA);

    if( cosB >  1.0) cosB =  1.0;
    if( cosB < -1.0) cosB = -1.0;
    angles[1] = acos(cosB);

    if( cosC >  1.0) cosC =  1.0;
    if( cosC < -1.0) cosC = -1.0;
    angles[2] = acos(cosC);

    if( measure ==  ANGLE_IN_DEGREE) {
        angles[0] *= 180/M_PI;
        angles[1] *= 180/M_PI;
        angles[2] *= 180/M_PI;
    }
}

inline double getTriAngle(const Point3D &pa, const Point3D &pb, const Point3D &pc, int measure = ANGLE_IN_DEGREES)
{
    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );
    double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );

    if( cosA >  1.0) cosA =  1.0;
    if( cosA < -1.0) cosA = -1.0;

    double angle;
    if( measure == ANGLE_IN_DEGREES)
        angle = 180*acos(cosA)/M_PI;
    else
        angle = acos(cosA);
    return angle;
}

////////////////////////////////////////////////////////////////////////////////

inline int getMaxTriAngle(const Point3D &pa, const Point3D &pb, const Point3D &pc, double &angle)
{
    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );

    double maxlen = max_value(a2,b2,c2);

    if( maxlen == a2) {
        double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );
        if( cosA >  1.0) cosA =  1.0;
        if( cosA < -1.0) cosA = -1.0;
        angle = 180*acos(cosA)/M_PI;
        return 0;
    }

    if( maxlen == b2 ) {
        double cosB =  (a2 + c2 - b2)/(2*sqrt(a2*c2) );
        if( cosB >  1.0) cosB =  1.0;
        if( cosB < -1.0) cosB = -1.0;
        angle  = 180*acos(cosB)/M_PI;
        return 1;
    }

    if( maxlen == c2 ) {
        double cosC =  (a2 + b2 - c2)/(2*sqrt(a2*b2) );
        if( cosC >  1.0) cosC =  1.0;
        if( cosC < -1.0) cosC = -1.0;
        angle = 180*acos(cosC)/M_PI;
        return 2;
    }

    return -1;
}

////////////////////////////////////////////////////////////////////////////////

inline int getMinTriAngle(const Point3D &pa, const Point3D &pb, const Point3D &pc, double &angle)
{
    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );

    double minlen = min_value(a2,b2,c2);

    if( minlen == a2) {
        double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );
        if( cosA >  1.0) cosA =  1.0;
        if( cosA < -1.0) cosA = -1.0;
        angle = 180*acos(cosA)/M_PI;
        return 0;
    }

    if( minlen == b2 ) {
        double cosB =  (a2 + c2 - b2)/(2*sqrt(a2*c2) );
        if( cosB >  1.0) cosB =  1.0;
        if( cosB < -1.0) cosB = -1.0;
        angle  = 180*acos(cosB)/M_PI;
        return 1;
    }

    if( minlen == c2 ) {
        double cosC =  (a2 + b2 - c2)/(2*sqrt(a2*b2) );
        if( cosC >  1.0) cosC =  1.0;
        if( cosC < -1.0) cosC = -1.0;
        angle = 180*acos(cosC)/M_PI;
        return 2;
    }

    return -1;
}

////////////////////////////////////////////////////////////////////////////////
