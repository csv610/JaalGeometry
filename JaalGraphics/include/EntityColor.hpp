#pragma once

#include <boost/assign/std/vector.hpp>

#include "Mesh.hpp"
using namespace std;
using namespace Jaal;
using namespace boost::assign;

#include <GL/gl.h>

typedef Array4F JColor;

// The following colors are obstained from
// http://en.wikipedia.org/wiki/Web_colors

const static GLfloat White[]  = { 1.00, 1.00, 1.00, 1.0};
const static GLfloat SnowWhite[]  = { 1.00, 0.97, 0.97, 1.0};
const static GLfloat Black[]  = { 0.00, 0.00, 0.00, 1.0};
const static GLfloat Red[]    = { 1.00, 0.00, 0.00, 1.0};
const static GLfloat Blue[]   = { 0.00, 0.00, 1.00, 1.0};
const static GLfloat Green[]  = { 0.00, 1.00, 0.00, 1.0};
const static GLfloat Lime[]   = { 0.00, 0.50, 0.00, 1.0};
const static GLfloat Silver[] = { 0.75, 0.75, 0.75, 1.0};
const static GLfloat Gray[]   = { 0.50, 0.50, 0.50, 1.0};
const static GLfloat Moroon[] = { 0.50, 0.00, 0.00, 1.0};
const static GLfloat Yellow[] = { 1.00, 1.00, 0.00, 1.0};
const static GLfloat Olive[]  = { 0.50, 0.50, 0.00, 1.0};
const static GLfloat Aqua[]   = { 0.00, 1.00, 1.00, 1.0};
const static GLfloat Teal[]   = { 0.00, 0.50, 0.50, 1.0};
const static GLfloat Navy[]   = { 0.00, 0.00, 0.50, 1.0};
const static GLfloat Fuchsia[]= { 1.00, 0.00, 1.00, 1.0};
const static GLfloat Purple[] = { 0.50, 0.00, 0.50, 1.0};
const static GLfloat DarkGreen[] = { 1.00, 0.19, 0.125, 1.0};

const static GLfloat LightGray[] = { 0.46,  0.53, 0.60, 1.0};
const static GLfloat LightBlue[] = { 0.675, 0.843, 0.898, 1.0};

///////////////////////////////////////////////////////////////////////////////
struct JColorMap
{
    static void jet( const vector<double> &a, vector<JColor> &c);
};

class JEntityColor {
public:
    static JColor getRandomColor();
    static JColor getColor( const string &s);
    static JColor getMinColor( int id );

    JEntityColor()
    {
        minColors = 0;
        alpha     = 1.0;
    }

    virtual ~JEntityColor() {}

    void   setAlpha( float  a ) {
        alpha = a;
    }

    float  getAlpha() const {
        return alpha;
    }

    void setHighlightColor( const JColor &clr) {
        highlightColor[0] = clr[0];
        highlightColor[1] = clr[1];
        highlightColor[2] = clr[2];
        highlightColor[3] = clr[3];
    }

    void setDefaultColor( const JColor &clr) {
        defaultColor[0] = clr[0];
        defaultColor[1] = clr[1];
        defaultColor[2] = clr[2];
        defaultColor[3] = clr[3];
    }

    void useMinColor( bool v) {
        minColors = v;
    }

protected:
    bool  minColors = 0;
    float alpha     = 0.5;
    JColor color, highlightColor, defaultColor;
};

///////////////////////////////////////////////////////////////////////////////

class JNodeColor;
typedef boost::shared_ptr<JNodeColor>  JNodeColorPtr;

class JNodeColor : public JEntityColor {
public:
    virtual ~JNodeColor() {}

    virtual string getName() const = 0;
    virtual bool isPerNode() const = 0;

    virtual int  assign(const JNodePtr &v) = 0;

    int  setMesh(const JMeshPtr &m);
    int  setScalarField(const JMeshPtr &m, const string &s);

    virtual int  buildRelations(const JMeshPtr &) {
        return 0;
    }
};
///////////////////////////////////////////////////////////////////////////////

class JEdgeColor;
typedef boost::shared_ptr<JEdgeColor>  JEdgeColorPtr;

struct JEdgeColor : public JEntityColor {
    static int   assign_scalar_field(const JMeshPtr &m, const string &s);
    virtual ~JEdgeColor() {}

    virtual string getName() const = 0;
    virtual bool isPerEdge() const = 0;
    virtual int  assign( const JEdgePtr &e)  = 0;

    int   setMesh(const JMeshPtr &m);

    virtual int buildRelations( const JMeshPtr &) {
        return 0;
    }
};
///////////////////////////////////////////////////////////////////////////////

class JFaceColor;
typedef boost::shared_ptr<JFaceColor>  JFaceColorPtr;

class JFaceColor : public JEntityColor {
public:
    static int  assign_scalar_field(const JMeshPtr &m, const string &s);

    virtual ~JFaceColor() {}

    virtual string getName() const = 0;
    virtual bool isPerFace() const = 0;
    virtual int  assign( const JFacePtr &e)  = 0;

    int  setMesh(const JMeshPtr &m);

    virtual int buildRelations( const JMeshPtr &) {
        return 0;
    }
};
///////////////////////////////////////////////////////////////////////////////

class JCellColor;
typedef boost::shared_ptr<JCellColor>  JCellColorPtr;

class JCellColor : public JEntityColor {
public:
    static int  assign_scalar_field(JMeshPtr m, const string &s);

    virtual ~JCellColor() {}

    virtual string getName() const = 0;
    virtual bool isPerCell() const = 0;
    virtual int  assign( const JCellPtr &c)  = 0;

    int  setMesh(const JMeshPtr &m);
};

///////////////////////////////////////////////////////////////////////////////

class JMeshNodeColor : public JNodeColor {
public:
    string getName() const {
        return "MeshNodeColor";
    }

    bool isPerNode() const {
        return 0;
    }

    JMeshNodeColor();

    void setInternalColor( float *rgb );

    const JColor &getInternalColor() const {
        return color;
    }

    void setBoundaryColor( float *rgb);
    const JColor &getBoundaryColor() const {
        return boundColor;
    }

    void setFeatureColor( float *rgb);
    const JColor &getFeatureColor() const {
        return featureColor;
    }

    void setInterfaceColor( float *rgb);
    const JColor &getInterfaceColor() const {
        return interfaceColor;
    }

    void displayFeatures( bool v )  {
        display_features  = v;
    }
    void displayInterface( bool v ) {
        display_interface = v;
    }
    void displayBoundary( bool v )  {
        display_boundary  = v;
    }
    void displayInternal( bool v )  {
        display_internal  = v;
    }

    int assign( const JNodePtr &vertex);

    int  operator() (const JNodePtr &v) {
        return assign(v);
    }

protected:
    bool  display_features, display_interface, display_boundary, display_internal;
    JColor boundColor, featureColor, interfaceColor;
};

///////////////////////////////////////////////////////////////////////////////

class IrregularNodeColor : public JNodeColor {
public:
    string getName() const {
        return "IrregularNodeColor";
    }

    bool isPerNode() const {
        return 1;
    }
    IrregularNodeColor();

    int buildRelations( const JMeshPtr &m);

    int assign( const JNodePtr &v);

private:
    JFaceSequence   faceneighs;
    void getFaceType( const  JNodePtr v, int &etype, int &degree ) {
        JNode::getRelations(v, faceneighs);
        if( faceneighs.empty() ) {
            cout << "Warning: Relation 02 absent " << endl;
            return;
        }
        degree = faceneighs.size();

        etype = faceneighs[0]->getTypeID();
        for( size_t  i = 1; i < faceneighs.size(); i++) {
            if( faceneighs[i]->getTypeID() != etype) {
                etype = 0;
                return;
            }
        }
    }
};

///////////////////////////////////////////////////////////////////////////////

class AngleDefectColor : public JNodeColor {
public:
    string getName() const {
        return "AngleDefectColor";
    }

    bool isPerNode() const {
        return 1;
    }

    AngleDefectColor();

    int assign(const JNodePtr &v);

    int  operator() (const JNodePtr &v) {
        return assign(v);
    }
};

///////////////////////////////////////////////////////////////////////////////

class EigenColor : public JNodeColor {
public:
    string getName() const {
        return "EigenColor";
    }

    bool isPerNode() const {
        return 1;
    }

    int assign(const JNodePtr &v);

    int  operator() (const JNodePtr &v) {
        return assign(v);
    }
};

///////////////////////////////////////////////////////////////////////////////

class JNodeDegreeColor : public JNodeColor {
public:
    JNodeDegreeColor() {
        lowDegree = 0;
        equalDegree = 0;
        highDegree = 0;
    }
    ~JNodeDegreeColor() {}

    string getName() const {
        return "NodeDegreeColor";
    }

    bool isPerNode() const {
        return 1;
    }

    void setLowDegree( int v )   {
        lowDegree = v;
    }
    void setEqualDegree( int v ) {
        equalDegree = v;
    }
    void setHighDegree( int v )  {
        highDegree  = v;
    }

    int assign(const JNodePtr &node);

    int  operator() (const JNodePtr &v) {
        return assign(v);
    }
private:
    int lowDegree, equalDegree, highDegree;
};

///////////////////////////////////////////////////////////////////////////////

class NodeHeightColor : public JNodeColor {
public:
    string getName() const {
        return "NodeHeightColor";
    }

    bool isPerNode() const {
        return 1;
    }

    void setMesh(const JMeshPtr &m, int dir );

    int assign(const JNodePtr &v);

    int  operator() (const JNodePtr &v) {
        return assign(v);
    }
private:
    double minH, maxH;
    int    dir;
};

///////////////////////////////////////////////////////////////////////////////
class ManifoldEdgeColor : public JEdgeColor {
public:
    string getName() const {
        return "ManifoldEdgeColor";
    }
    bool isPerEdge() const {
        return 1;
    }

    int  assign( const JEdgePtr &edge);

    int  operator() (const JEdgePtr &e) {
        return assign(e);
    }
};



///////////////////////////////////////////////////////////////////////////////
class JMeshEdgeColor : public JEdgeColor {

public:
    JMeshEdgeColor();

    string getName() const {
        return "MeshEdgeColor";
    }
    bool isPerEdge() const {
        return 0;
    }

    void setInternalColor( float *rgb );

    const JColor &getInternalColor() const {
        return color;
    }

    void setBoundaryColor( float *rgb);
    const JColor &getBoundaryColor() const {
        return boundColor;
    }

    void setFeatureColor( float *rgb);
    const JColor &getFeatureColor() const {
        return featureColor;
    }

    void setInterfaceColor( float *rgb);
    const JColor &getInterfaceColor() const {
        return interfaceColor;
    }

    void displayFeatures( bool v ) {
        display_features  = v;
    }
    void displayInterface( bool v ) {
        display_interface = v;
    }
    void displayBoundary( bool v ) {
        display_boundary  = v;
    }
    void displayInternal( bool v ) {
        display_internal  = v;
    }

    int assign( const JEdgePtr &e);

    int  operator() (const JEdgePtr &e) {
        return assign(e);
    }

protected:
    JColor boundColor, featureColor, interfaceColor;
    bool  display_internal, display_boundary, display_interface, display_features;

    bool isFeature( const JEdgePtr edge ) {
        if( edge->hasAttribute( "Feature")  ) return 1;
        if( edge->hasAttribute( "SharpEdge")) return 1;
        return 0;
    }
};

///////////////////////////////////////////////////////////////////////////////

class QuadParamLinesColor : public JEdgeColor {

public:
    string getName() const {
        return "QuadParamLinesColor";
    }

    bool isPerEdge() const {
        return 1;
    }


    QuadParamLinesColor();

    int assign( const JEdgePtr &edge);

    int  operator() (const JEdgePtr &e) {
        return assign(e);
    }

protected:
    JColor xcolor, ycolor, defaultcolor;
};

////////////////////////////////////////////////////////////////////////////////

class EdgeInterfaceColor : public JEdgeColor {
public:
    string getName() const {
        return "EdgeInterfaceColor";
    }

    bool isPerEdge() const {
        return 0;
    }

    void assignColors( const JMeshPtr &m);

    int assign(const JEdgePtr &edge);

    int  operator() (const JEdgePtr &e) {
        return assign(e);
    }

private:
    typedef std::pair<int,int> ipair_t;
    ipair_t ipair;
    map<ipair_t, JColor>  colormap;
};

///////////////////////////////////////////////////////////////////////////////
class SharpEdgeColor : public JEdgeColor {
public:
    string getName() const {
        return "EdgeFeatureColor";
    }

    bool isPerEdge() const {
        return 0;
    }

    int assign(const JEdgePtr &edge);
};

///////////////////////////////////////////////////////////////////////////////

class JMeshFaceColor : public JFaceColor {
public:
    string getName() const {
        return "MeshFaceColor";
    }

    bool isPerFace() const {
        return 0;
    }

    JMeshFaceColor();

    void setFrontColor( float *rgb);
    const JColor &getFrontColor() const {
        return frontColor;
    }

    void setBackColor( float *rgb);
    const JColor &getBackColor() const {
        return backColor;
    }

    void setInternalColor( float *rgb );
    void setInternalColor( const JColor &c) {
        color = c;
    }

    JColor getInternalColor() const {
        return color;
    }

    void setBoundaryColor( float *rgb);
    void setBoundaryColor( const JColor &c) {
        boundColor = c;
    }

    const JColor &getBoundaryColor() const {
        return boundColor;
    }

    void setInterfaceColor( float *rgb);
    const JColor &getInterfaceColor() const {
        return interfaceColor;
    }

    void displayInterface( bool v ) {
        display_interface = v;
    }

    void displayBoundary( bool v )  {
        display_boundary  = v;
    }
    void displayInternal( bool v )  {
        display_internal  = v;
    }

    void setCell( JCellPtr c ) {
        cell = c;
    }
    JCellPtr getCell() const {
        return cell;
    }

    int assign( const JFacePtr &face);

    int  operator() (const JFacePtr &f) {
        return assign(f);
    }
private:
    JColor frontColor, backColor, boundColor, featureColor, interfaceColor;
    bool  display_interface, display_boundary, display_internal;

protected:
    JCellPtr cell;
};
///////////////////////////////////////////////////////////////////////////////
class JMeshCellColor : public JCellColor {
public:
    JMeshCellColor();

    string getName() const {
        return "MeshCellColor";
    }
    bool isPerCell() const {
        return 0;
    }

    void setDefaultColor( float *rgb );

    int assign(const JCellPtr &c);

    int  operator() (const JCellPtr &c) {
        return assign(c);
    }
};

///////////////////////////////////////////////////////////////////////////////

class PenroseTileColor : public JFaceColor {
public:
    string getName() const {
        return "PenroseColor";
    }

    bool isPerFace() const {
        return 1;
    }

    int assign( const JFacePtr &face);

    int  operator() (const JFacePtr &f) {
        return assign(f);
    }
};

///////////////////////////////////////////////////////////////////////////////

class QuadDiamondColor : public JFaceColor {
public:

    string getName() const {
        return "QuadDiamondColor";
    }

    bool isPerFace() const {
        return 1;
    }

    QuadDiamondColor();

    int assign(const JFacePtr &face);

    int  operator() (const JFacePtr &f) {
        return assign(f);
    }
};

