#pragma once

#include "basic_math.hpp"

/*
Plane.h
Written by Matthew Fisher

A standard 3D plane (space plane;) i. e. the surface defined by a*x + b*y + c*z + d = 0
*/

struct JLine3D
{
    Point3D  p0;
    Point3D  p1;
};

struct JPlane
{
    static double dot(const JPlane &P, const Vec4D &V);         //dot product of a plane and a 4D vector
    static double dotCoord(const JPlane &P, const Vec3D &V);   //dot product of a plane and a 3D coordinate
    static double dotNormal(const JPlane &P, const Vec3D &V);  //dot product of a plane and a 3D normal
    static bool   JPlaneJPlaneIntersection(const JPlane &P1, const JPlane &P2, JLine3D &L);

    // loads the plane from a point on the surface and a normal vector
    static JPlane fromPointNormal(const Vec3D &Pt, const Vec3D &Normal);

    //loads the plane from a point on the surface and two vectors in the plane
    static JPlane fromPointVectors(const Vec3D &Pt, const Vec3D &V1, const Vec3D &V2);

    //loads the plane from 3 points on the surface
    static JPlane fromPoints(const Vec3D &V1, const Vec3D &V2, const Vec3D &V3);

    JPlane();
    JPlane(const JPlane &P);
    JPlane(double _a, double _b, double _c, double _d);
    JPlane(const Vec3D &NormalizedNormal, double _d);

    double getUnsignedDistance(const Vec3D &Pt) const;
    double getSignedDistance(const Vec3D &Pt) const;
    bool   fitToPoints(const vector<Vec3D> &Points, double &ResidualError);
    bool   fitToPoints(const vector<Vec4D> &Points, Vec3D &Basis1, Vec3D &Basis2,
                       double &NormalEigenvalue, double &ResidualError);
    Vec3D  getClosestPoint(const Vec3D &Point);

    // determines the intersect of the line defined by the points V1 and V2 with the plane.
    // Returns the point of intersection.  Origin is returned if no intersection exists.
    Vec3D intersectLine(const Vec3D &V1, const Vec3D &V2) const;

    // determines the intersect of the line defined by the points V1 and V2 with the plane.
    // If there is no intersection, Hit will be false.
    Vec3D intersectLine(const Vec3D &V1, const Vec3D &V2, bool &Hit) const;

    // Paramaterize the line with the variable t such that t = 0 is V1 and t = 1 is V2.
    // returns the t for this line that lies on this plane.
    double intersectLineRatio(const Vec3D &V1, const Vec3D &V2);

    Vec3D intersectLine(const JLine3D &Line) const;

    // Normalization
    //
    JPlane normalize();

    /*
        Vec3D normal() const
        {
            return Vec3D(a, b, c);
        }
    */

    JPlane flip()
    {
        JPlane Result;
        Result.a = -a;
        Result.b = -b;
        Result.c = -c;
        Result.d = -d;
        return Result;
    }

    double a, b, c, d;        //the (a, b, c, d) in a*x + b*y + c*z + d = 0.
};

