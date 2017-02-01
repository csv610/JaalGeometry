#pragma once

#ifndef BARYCOORDS_H
#define BARYCOORDS_H

#include "Mesh.hpp"

////////////////////////////////////////////////////////////////////////////////

struct JBarycentricCoordinates
{

    /*
    void setGreenBary2GlobalCoords( const vector<Vertex*> &cageNodes,
                                    const vector<double>  &orgseglen,
                                    const vector<double>  &newseglen,
                                    const vector<Point2D> &cageNormals,
                                    Vertex *anyNode );
    void setBary2GlobalCoords( const vector<Vertex*> &cageNodes, Vertex *anyNode );
    void GreenCoordinates( const vector<Vertex*>    & modelNodes,
                           const vector<Vertex*>    & cageNodes , bool genfield = 0);
    void MeanValueCoordinates(const vector<Vertex*> & modelNodes,
                              const vector<Vertex*> & cageNodes, bool genfield = 0 );
    void HarmonicCoordinates(const vector<Vertex*>   & modelNodes,
                             const vector<Vertex*>   & cageNodes, bool genfield = 1);
    */

    int setXYCoords(const JMeshPtr &m, const vector<Point2D> &cageCoords);
    int getXYCoords(const vector<Point2D> &cageCoords, const vector<double> &baryCoords, Point2D &xy);

    int setFieldValues(const vector<Point2D> &cageCoords, const vector<double> &cageValues, const JMeshPtr &mesh, const string &s);

    virtual int     getCoords( const vector<Point2D> &cageCoords, const Point2D &qPoint, vector<double> &c) = 0;
    virtual double  getFieldValueAt( const vector<Point2D> &cageCoords, const vector<double> &scalarfield, const Point2D &qPoint) = 0;

    int  quadShapeFunc( const vector<Point2D> &quadPoints, const Point2D &uv, vector<double> &shapefunc);
};

////////////////////////////////////////////////////////////////////////////////

struct JMeanValueCoordinates : public JBarycentricCoordinates
{
    typedef boost::shared_ptr<JMeanValueCoordinates>  SharedPtr;
    static  SharedPtr  newObject();

    int     getCoords( const vector<Point2D> &cageCoords, const Point2D &qPoint, vector<double> &c);
    double  getFieldValueAt( const vector<Point2D> &cageCoords, const vector<double> &scalarfield, const Point2D &qPoint);
    int     getShapeGradients( const vector<Point2D> & cageCoords, const Point2D &qPoint, vector<Vec2D> &grad);

    int     getW( const vector<Point2D> &cageCoords, const Point2D &qPoint, vector<double> &w);
    int     getGradW( const vector<Point2D> &cageCoords, const Point2D &qPoint, vector<Vec2D> &grad);
    void    setBaryCoords(const JMeshPtr &m, const vector<Point2D> &cageCoords);
};

struct JGreenCoordinates : public JBarycentricCoordinates
{
    typedef boost::shared_ptr<JGreenCoordinates>  SharedPtr;
    static  SharedPtr  newObject();

    int     getCoords( const vector<Point2D> &cageCoords, const Point2D &qPoint, vector<double> &c);
    double  getFieldValueAt( const vector<Point2D> &cageCoords, const vector<double> &scalarfield, const Point2D &qPoint);

};
////////////////////////////////////////////////////////////////////////////////

/*
struct JHarmonicCoordinates : public JBarycentricCoords
{
    static int  getCoords( const vector<Point2D> &cageCoords, const Point2D &qPoint, vector<double> &c);
    static double  getFieldValueAt( const vector<Point2D> &cageCoords, const vector<double> &scalarfield, const Point2D &qPoint);
};
*/

#endif
