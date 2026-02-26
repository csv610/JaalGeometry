#pragma once

#include <string>
#include "Attributes.hpp"

struct JShape
{
    template<class T>
    int setAttribute(const string &s, const T &val)
    {
        return attribManager->setAttribute(s,val);
    }

    template<class T>
    int getAttribute(const string &s, T &val) const
    {
        return attribManager->getAttribute(s,val);
    }

    int getNumAttributes() const
    {
        return attribManager->getNumAttributes();
    }

    int getAttributeNames( vector<string> &names) const
    {
        names.clear();
        return attribManager->getAttributeNames( names );
    }

    int hasAttribute(const string &s) const
    {
        return attribManager->hasAttribute(s);
    }

    void  deleteAttribute(const string &s)
    {
        attribManager->deleteAttribute(s);
    }

    JAttributeManagerPtr attribManager;
};

class JCircle : public JShape
{
public:
    JCircle()
    {
        center[0] = 0.0;
        center[1] = 0.0;
        center[2] = 0.0;
        radius = 0.0;
    }

    JCircle(double r)
    {
        center[0] = 0.0;
        center[1] = 0.0;
        center[2] = 0.0;
        radius = r;
    }
    ~JCircle() {}

    void setCenter( const Point3D &p)
    {
        center = p;
    }
    const Point3D &getCenter() const
    {
        return center;
    }

    void setRadius(double r)
    {
        radius = r;
    }

    double getRadius() const
    {
        return radius;
    }

    double getArea() const
    {
        return M_PI*radius*radius;
    }
    double getDiameter() const
    {
        return 2*radius;
    }
    double getCircumference() const
    {
        return 2*M_PI*radius;
    }
private:
    double radius;
    Point3D center;
};

typedef boost::shared_ptr<JCircle> JCirclePtr;

///////////////////////////////////////////////////////////////////////////////

class JSphere : public JShape
{
public:
    JSphere()
    {
        center[0] = 0.0;
        center[1] = 0.0;
        center[2] = 0.0;
        radius = 0.0;
    }

    JSphere(double r)
    {
        center[0] = 0.0;
        center[1] = 0.0;
        center[2] = 0.0;
        radius = r;
    }

    void setCenter( const Point3D &p)
    {
        center = p;
    }

    const Point3D &getCenter() const
    {
        return center;
    }

    void setRadius(double r)
    {
        radius = r;
    }

    double getRadius()
    {
        return radius;
    }

    double getArea()
    {
        return 4*M_PI*radius*radius;
    }

    double getVolume()
    {
        return 4*M_PI*radius*radius*radius/3.0;
    }

private:
    double radius;
    Point3D center;
};
typedef boost::shared_ptr<JSphere> JSpherePtr;

///////////////////////////////////////////////////////////////////////////////

class JCylinder : public JShape
{
public:
    void setEndPoints( const Point3D &p0, const Point3D &p1)
    {
        pStart = p0, pEnd = p1;
    }

    void setRadius( double r) {
        radius = r;
    }

    void setOrigin( const Point3D &p) {
        pStart = p;
    }

    void setAxis( const Vec3D &v, double h) {
        pEnd[0] = pStart[0] + h*v[0];
        pEnd[1] = pStart[1] + h*v[1];
        pEnd[2] = pStart[2] + h*v[2];
    }
private:
    Point3D pStart, pEnd;
    double  radius, height;
};
typedef boost::shared_ptr<JCylinder> JCylinderPtr;

///////////////////////////////////////////////////////////////////////////////

class JCone : public JShape
{
public:
    void setOrigin( const Point3D &p) {
        origin = p;
    }
    void setHeight( double h)  {
        height = h;
    }
    void setRadius( double r)  {
        radius  = r;
    }
    void setAxis( const Vec3D &v) {
        axis = v;
    }
private:
    bool cap = 1;
    Point3D origin;
    Vec3D   axis;
    double  radius, height;
};
typedef boost::shared_ptr<JCone> JConePtr;

///////////////////////////////////////////////////////////////////////////////

class JBox : public JShape
{
public:
    void setHeight( double h)  {
        height = h;
    }
    void setWidth( double w)   {
        width  = w;
    }
    void setCenter( const Point3D &p) {
        center = p;
    }
    void setAxis( const Vec3D &v)     {
        axis   = v;
    }
private:
    Point3D center;
    Vec3D   axis;
    double  height, width;
};
typedef boost::shared_ptr<JBox> JBoxPtr;

///////////////////////////////////////////////////////////////////////////////

class JTube : public JShape
{
public:
    void setObject( const vector<Point3D> &seq) {
        points = seq;
    }
    void setRadius( double r) {
        radius = r;
    }
private:
    vector<Point3D> points;
    double  radius;
};
typedef boost::shared_ptr<JTube> JTubePtr;

///////////////////////////////////////////////////////////////////////////////
