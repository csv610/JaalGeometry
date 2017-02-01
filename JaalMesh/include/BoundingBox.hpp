#pragma once

#include "basic_math.hpp"

namespace Jaal {
class JBoundingBox {
public:
    JBoundingBox()
    {
        lowerCorner[0] = 0.0;
        lowerCorner[1] = 0.0;
        lowerCorner[2] = 0.0;
        upperCorner[0] = 0.0;
        upperCorner[1] = 0.0;
        upperCorner[2] = 0.0;
    }

    JBoundingBox( const Point3D &p0, const Point3D &p1)
    {
        setPoints(p0,p1);
    }

    int getDimension() const
    {
        int d = 0;
        for( int i = 0; i < 3; i++)
            if( getLength(i) ) d++;
        return d;
    }

    void setCenter( const Point3D &pc, double len)
    {
        lowerCorner[0] = pc[0] - 0.5*len;
        lowerCorner[1] = pc[0] - 0.5*len;
        lowerCorner[2] = pc[0] - 0.5*len;

        upperCorner[0] = pc[0] + 0.5*len;
        upperCorner[1] = pc[0] + 0.5*len;
        upperCorner[2] = pc[0] + 0.5*len;
    }

    void setPoints( const Point3D &p0, const Point3D &p1)
    {
        for(int i = 0; i < 3; i++) {
            lowerCorner[i] = std::min(p0[i], p1[i] );
            upperCorner[i] = std::max(p0[i], p1[i] );
        }
    }

    double getLength(int d) const
    {
        assert(d >= 0 && d < 3);
        return fabs(upperCorner[d] - lowerCorner[d]);
    }

    double getDiameter() const
    {
        return JMath::length( upperCorner, lowerCorner);
    }

    double getMaxLength() const
    {
        double dx = getLength(0);
        double dy = getLength(1);
        double dz = getLength(2);
        return std::max( dx, std::max(dy,dz));
    }

    double getMinLength() const
    {
        double dx = getLength(0);
        double dy = getLength(1);
        double dz = getLength(2);
        return std::min( dx, std::min(dy,dz));
    }

    Point3D getCenter() const
    {
        Point3D pc;
        pc[0] = 0.5*(lowerCorner[0] + upperCorner[0] );
        pc[1] = 0.5*(lowerCorner[1] + upperCorner[1] );
        pc[2] = 0.5*(lowerCorner[2] + upperCorner[2] );
        return pc;
    }

    Point3D getCorner( int ic ) const;

    void expandBy(double l);

    void setLower(const Point3D & p)
    {
        lowerCorner = p;
    }

    const Point3D & getLower() const
    {
        return lowerCorner;
    }

    void setUpper(const Point3D & p)
    {
        upperCorner = p;
    }

    const Point3D & getUpper() const
    {
        return upperCorner;
    }

    void getEdges(vector<Point3D> &edges) const;

    static int  getOrientation( const JBoundingBox &b, const Point3D &p);
    static bool intersect( const JBoundingBox &b1, const JBoundingBox &b2);
    void   setUnion( const JBoundingBox &b);

    friend ostream& operator << (ostream &os, const JBoundingBox &b);

private:
    Point3D lowerCorner;
    Point3D upperCorner;
};

inline
ostream& operator << (ostream &os, const JBoundingBox &box)
{
    Point3D p3d;

    p3d = box.getCorner(0);
    os << "{" << p3d[0] << " " << p3d[1] << " " << p3d[2] << "},";
    p3d = box.getCorner(6);
    os << "{" << p3d[0] << " " << p3d[1] << " " << p3d[2] << "}";
    os << endl;

    return os;
}
}
