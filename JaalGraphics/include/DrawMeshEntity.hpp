#pragma once

#include <GL/glu.h>
#include <boost/shared_ptr.hpp>

#include "Material.hpp"
#include "FontsManager.hpp"
#include "DrawMeshEntity.hpp"

class JaalViewer;

///////////////////////////////////////////////////////////////////////////////
struct JRender
{
    static const int   POINTCLOUD    = 0;
    static const int   WIREFRAME     = 1;
    static const int   HIDDENLINES   = 2;
    static const int   FLAT_SHADE    = 3;
    static const int   SMOOTH_SHADE  = 4;
    static const int   SURFACE_AND_WIREFRAME  = 5;

    static const int HIDDENLINE_WITH_DEPTH_TEST      = 0;
    static const int HIDDENLINE_WITH_BACKFACES_CULL  = 1;
    static const int HIDDENLINE_WITH_FRONTLINES      = 2;
    static const int HIDDENLINE_WITH_STENCIL_BUFFER  = 3;

    static const int TRANSPARENT_SCREENDOOR  = 0;
    static const int TRANSPARENT_BLEND       = 1;
};

//////////////////////////////////////////////////////////////////////////////

struct JMeshEntityRender {

    JMeshEntityRender()
    {
        fontAngles[0] = 0.0;
        fontAngles[1] = 0.0;
        fontAngles[2] = 0.0;
    }

    virtual ~JMeshEntityRender() {}

    void setViewManager( JaalViewer *v) {
        viewManager = v;
    }

    void setDisplay(bool v)
    {
        display = v;
    }

    bool getDisplay() const
    {
        return display;
    }

    void setMesh( const JMeshPtr &m)
    {
        mesh = m;
    }

    void setFont( FTFont *f, const JColor &c, float s)
    {
        font = f;
        fontColor = c;
        fontScale = s;
    }

    void setRotateFonts( const Point3F &angles)
    {
        fontAngles = angles;
    }

    void drawTextAt( char *text, const float *p) const;

    void setMaterial( const string &s)
    {
        materialName = s;
    }

    void setMaterial(JMaterial *jm)
    {
        material = jm;
    }
    void   setAlpha(double a)
    {
        alpha = a;
    }
    double getAlpha() const
    {
        return alpha;
    }

    FTFont *font = nullptr;
    JColor   fontColor;
    Point3F fontAngles;
    string  materialName;

    bool    display = 1;
    float   fontScale = 1.0;
    double  alpha     = 1.0;
    JMeshPtr  mesh;

    JMaterial *material = nullptr;
    JaalViewer *viewManager = nullptr;

    mutable JColor  color;
    bool   pickable = 1;
};

///////////////////////////////////////////////////////////////////////////////

struct JNodeRender 
{
    static boost::shared_ptr<JNodeRender> newObject();

    static const int NODE_AS_POINT   = 0;
    static const int NODE_AS_SPHERE  = 1;
    static const int NODE_AS_SPLAT   = 2;

    JNodeRender()
    {
        color[0]   = 0.8;
        color[1]   = 0.0;
        color[2]   = 0.0;
        color[3]   = 1.0;
    }

    void setGlyph( int g )
    {
        if( g >= 0 && g  <= 3) glyph = g;
    }

    bool   display  = 1;
    short  glyph    = 0;
    float  scale   = 1.0;
    float  pointSize  = 1.0;
    float  ballRadius = 1.0;
    size_t index;
    JColor color;
};

typedef boost::shared_ptr<JNodeRender> JNodeRenderPtr;

///////////////////////////////////////////////////////////////////////////////

struct JNodeBuffers {
    void clearAll()
    {
        pointGroup.clear();
        sphereGroup.clear();
        splatGroup.clear();
        rootGroup.clear();

        coordsArray.clear();
        colorArray.clear();
        normalArray.clear();
    }

    std::map<double,vector<int> > pointGroup;
    std::map<double,vector<int> > sphereGroup;
    std::map<double,vector<int> > splatGroup;
    vector<int>    rootGroup;
    vector<float>  coordsArray;
    vector<float>  colorArray;
    vector<float>  normalArray;
};

///////////////////////////////////////////////////////////////////////////////

class JNodeDraw : public JMeshEntityRender {
public:
    static GLUquadricObj *sphereObj;
    static GLUquadricObj *diskObj;

    JNodeDraw();

    virtual ~JNodeDraw() {}

    JNodeRenderPtr initRenderAttribute(const JNodePtr &v);

    void  initRenderAttributes( const JMeshPtr &m);

    void  setOffset( bool v )
    {
        offset = v;
    }

    void setTransparency( double a)
    {
        alpha = a;
    }
    void setBorderColor( const JColor &c)
    {
        borderColor = c;
    }

    void setBorderThickness( double t)
    {
        borderThickness = t;
    }

    void  setOffset( float factor, float unit )
    {
        offset_factor = factor;
        offset_unit   = unit;
    }

    void setAntiAliasing( bool v )
    {
        antiAlias = v;
    }

    void setGlyph( int g )
    {
        glyph = g;
    }

    int  getGlyph() const
    {
        return glyph;
    }

    void setPointSize( double r )
    {
        pointSize = r;
    }

    double getPointSize() const
    {
        return pointSize;
    }

    void setBallRadius( double r )
    {
        ballRadius  = r;
    }

    double getBallRadius() const
    {
        return ballRadius;
    }

    void setSphereResolution( int nslices, int mstacks )
    {
        numSlices = nslices;
        numStacks = mstacks;
    }

    // What is the length of the normals, if they are displayed ...
    void setNormalsLength( float s )
    {
        normalLength = s;
    }

    double getNormalsLength() const
    {
        return normalLength;
    }

    // Set the color of the nodes...
    void setNormalsColor( const JColor &c)
    {
        normalsColor = c;
    }

    const JColor &getNormalsColor() const
    {
        return normalsColor;
    }

    void  draw( const JNodePtr &vertex );
    void  drawID(const JNodePtr &v);
    void  drawNormal( const JNodePtr &v);
    void  withName(const JNodePtr &vertex, size_t pos);

    void draw(const JMeshPtr &mesh);
    void drawIDs(const JMeshPtr &mesh);
    void drawNormals(const JMeshPtr &mesh);
    void withName( const JMeshPtr &mesh );

    void preRender();
    void postRender();

    void updateBuffers(const JMeshPtr &mesh);

private:
    int    glyph;
    bool   offset;
    bool   lights;
    bool   antiAlias;
    int    numSlices, numStacks;  // For the circle and sphere
    float  offset_factor, offset_unit;
    float  pointSize, ballRadius;
    float  normalLength;
    float  alpha;
    float  borderThickness;
    JColor  normalsColor;
    JColor  borderColor;

    bool isDrawable( const JNodePtr &vertex) const;
    void setAttributes(const JNodePtr &vertex) const;
    void calculate_normal( const JNodePtr &v);
};

////////////////////////////////////////////////////////////////////////////////

struct JEdgeRender
{
    static boost::shared_ptr<JEdgeRender> newObject();

    JEdgeRender()
    {
        glyph    = 0;
        lineWidth = 1.0;
        cylinderRadius = 1.0;
        scale    = 1.0;
        color[0] = 0.1;
        color[1] = 0.1;
        color[2] = 0.1;
        color[3] = 1.0;
    }

    bool   display = 1;
    bool   pickable = 0;
    short  glyph = 0;
    float  scale = 1;
    float  cylinderRadius;
    float  lineWidth = 1.0;
    JColor color;
};

////////////////////////////////////////////////////////////////////////////////

struct JEdgeBuffers
{
    void clearAll()
    {
        lineGroup.clear();
        cylGroup.clear();
        nodesArray.clear();
        coordsArray.clear();
        colorArray.clear();
        lineGroup.clear();
        cylGroup.clear();
    }

    map<int,int>  edgePos;
    vector<int>   nodesArray;
    vector<float> coordsArray;
    vector<float> colorArray;
    map<double,vector<int> > lineGroup, cylGroup;
};

typedef boost::shared_ptr<JEdgeRender> JEdgeRenderPtr;

////////////////////////////////////////////////////////////////////////////////

class JEdgeDraw : public JMeshEntityRender {
public:
    static const int EDGE_AS_LINE     = 0;
    static const int EDGE_AS_CYLINDER = 1;

    JEdgeDraw();

    JEdgeRenderPtr initRenderAttribute( const JEdgePtr &v);
    void  initRenderAttributes( const JMeshPtr &m);

    virtual ~JEdgeDraw() {};

    void setGlyph( int g )
    {
        glyph = g;
    }
    int  getGlyph() const
    {
        return glyph;
    }

    void  setOffset( bool v )
    {
        offset = v;
    }

    void setOffset( float factor, float unit )
    {
        offset_factor = factor;
        offset_unit   = unit;
    }

    void setAntiAliasing( bool v )
    {
        antiAlias = v;
    }

    void setLineWidth( float r )
    {
        lineWidth = r;
    }

    void setCylinderRadius( float r )
    {
        cylRadius = r;
    }

    void setNumCylinderSides( int n )
    {
        numCylSides = n;
    }
    int getNumCylinderSides() const
    {
        return numCylSides;
    }

    double getLineWidth() const
    {
        return lineWidth;
    }

    double getCylinderRadius() const
    {
        return cylRadius;
    }

    double getRadius() const
    {
        if( glyph == EDGE_AS_LINE) return lineWidth;
        return cylRadius;
    }

    void  drawCylinder( const float *p0, const float *p1, float radius) const;
    void  drawCylinder( const JEdgePtr &edge, double radius) const;

    void  draw( const JEdgePtr &edge);
    void  drawID( const JEdgePtr &e);
    void  withName(const JEdgePtr &edge, size_t pos);

    void  draw( const JMeshPtr &mesh);
    void  drawIDs(const JMeshPtr &mesh);
    void  withName(const JMeshPtr &mesh);

    void  updateGeometryBuffers( const JMeshPtr &mesh);
    void  updateTopologyBuffers( const JMeshPtr &mesh);
    void  updateBuffers( const JMeshPtr &mesh);

private:
    bool  offset;
    bool  antiAlias;
    bool  withSteiner;
    int   numCylSides;
    int   glyph;
    float offset_factor, offset_unit;
    float lineWidth, cylRadius;

    boost::shared_ptr<JEdgeColor> colorMethod;
    boost::shared_ptr<JEdgeColor> defaultColorMethod;

    void  preRender();
    void  postRender();
    void  activateLowerEntities(const JEdgePtr &c);
};

////////////////////////////////////////////////////////////////////////////////

struct JFaceRender
{
    static boost::shared_ptr<JFaceRender> newObject();

    JFaceRender()
    {
        color[0] = 0.2;
        color[1] = 0.2;
        color[2] = 0.8;
        color[3] = 1.0;
    }

    bool       useMaterial = 0;
    bool       pickable    = 1;
    bool       display     = 1;
    bool       displayID   = 0;
    short int  faceSide    = 0;

    JColor    color, backColor;
    JMaterial material, backMaterial;
};

typedef boost::shared_ptr<JFaceRender> JFaceRenderPtr;

struct JFaceBuffers {
    void clearAll()
    {
        idArray.clear();
        triArray.clear();
        quadArray.clear();
        polyArray.clear();
        nodeCoords.clear();
        faceColor.clear();
        faceNormal.clear();
        faceNormalTail.clear();
        faceNormalHead.clear();
    }

    size_t getSize() const {
        return idArray.size();
    }

    vector<bool> textID;
    vector<int>  idArray;
    vector<int>  triArray, quadArray;
    vector<vector<int> >  polyArray;
    vector<float> nodeCoords, uvCoords;
    vector<float> faceColor;
    vector<float> faceNormal, faceNormalTail, faceNormalHead;

    vector<float> nodeColor;
    vector<float> nodeNormal;
};

////////////////////////////////////////////////////////////////////////////////

class JFaceDraw : public JMeshEntityRender {

public:
    JFaceDraw();

    virtual ~JFaceDraw() {};

    JFaceRenderPtr initRenderAttribute( const JFacePtr &v);
    void  initRenderAttributes( const JMeshPtr &m);
    void setLights( bool l ) {
        lights = l;
    }
    void  setOffset( bool v )
    {
        offset = v;
    }
    void setFrontFace( int f) {
        frontface = f;
    }

    void  setOffset( float factor, float unit )
    {
        offset_factor = factor;
        offset_unit   = unit;
    }

    void setScreendoorPattern( int p )
    {
        screendoorPattern = p;
    }

    void setCulling( bool v )
    {
        culling = v;
    }

    void preRender();
    void postRender();

    void useNormal( bool v )
    {
        use_normal = v;
    }

    void setNormalsLength(  float s )
    {
        normalLength = s;
    }
    double getNormalsLength() const
    {
        return normalLength;
    }

    void setNormalsColor( const JColor &c)
    {
        normalsColor = c;
    }

    const JColor &getNormalsColor() const
    {
        return normalsColor;
    }

    void  setRenderSide( int s = 1) {
        renderSide = s;
    }

    void  draw(const JFacePtr &face);
    void  drawID( const JFacePtr &v);
    void  withName(const JFacePtr &face, size_t pos);
    void  drawNormal( const JFacePtr &v);

    void draw( const JMeshPtr &mesh);
    void drawIDs( const JMeshPtr &mesh);
    void drawNormals( const JMeshPtr &mesh);
    void drawHiddenLines(const JMeshPtr &m);
    void withName( const JMeshPtr &mesh);

    void  updateGeometryBuffers( const JMeshPtr &mesh);
    void  updateTopologyBuffers( const JMeshPtr &mesh);
    void  updateBuffers( const JMeshPtr &m);

private:
    bool    offset;
    float   offset_factor, offset_unit;
//  int     shade;
    bool    use_normal;
    float   normalLength;
    bool    culling;
    int     transparencyMethod;
    int     screendoorPattern;
    JColor  normalsColor;
    bool    lights;
    int     renderSide;
    int     frontface;

    vector<int>    idArray;
    vector<int>    triArray, quadArray, polyArray, polyOffsetArray;
    vector<float>  posArray;
    vector<float>  colorArray;
    vector<float>  normalArrowArray;

    boost::shared_ptr<JFaceColor> colorMethod;

    void setAttributes( const JFacePtr &face )const;
    void setNormal( const JFacePtr &face )const;

    void drawPolySurface(const JMeshPtr &m);
    void drawFlatSurface(const JMeshPtr &m);
    void drawSmoothSurface(const JMeshPtr &m);
    void drawHiddenlines(const JMeshPtr &m);
    void drawTexturedSurface(const JMeshPtr &m);

private:
    void draw_tri( const JFacePtr &face, int ori);
    void draw_quad( const JFacePtr &face,int ori);
    void draw_poly( const JFacePtr &face,int ori);
    void activateLowerEntities(const JFacePtr &c);
};

typedef boost::shared_ptr<JFaceRender> JFaceRenderPtr;

/////////////////////////////////////////////////////////////////////////////

class JCellDraw : public JMeshEntityRender {
public:
    JCellDraw();

    virtual ~JCellDraw() {};

    void  draw(const JCellPtr &cell);
    void  withName(const JCellPtr &cell, size_t pos);
    void  drawID(const JCellPtr &cell);

    void draw( const JMeshPtr &mesh);
    void withName( const JMeshPtr &mesh);
    void drawIDs( const JMeshPtr &mesh);

    void  initRenderAttributes( const JMeshPtr &m);
    void  updateBuffers( const JMeshPtr &m);

private:
    JFaceSequence faces;
    JEdgeSequence edges;
    boost::shared_ptr<JCellColor> colorMethod;

    void setAttributes( const JCellPtr &c )const;
    bool isDrawable( const JCellPtr &c) const;
    void activateLowerEntities(const JCellPtr &c);
};

struct JCellRender {
    static boost::shared_ptr<JCellRender> newObject();
    JCellRender()
    {
        color[0] = 0.2;
        color[1] = 0.2;
        color[2] = 0.8;
        color[3] = 1.0;
        display  = 1;
        pickable = 1;
    }
    bool       display;
    bool       pickable;
    JColor color,  backColor;
};

typedef boost::shared_ptr<JCellRender> JCellRenderPtr;

/////////////////////////////////////////////////////////////////////////////
struct JMeshRender
{
    JMeshRender()
    {
        displayEntity[0] = 1;
        displayEntity[1] = 1;
        displayEntity[2] = 1;
        displayEntity[3] = 1;
        displayIDs[0] = 0;
        displayIDs[1] = 0;
        displayIDs[2] = 0;
        displayIDs[3] = 0;
        displayNormals[0] = 0;
        displayNormals[1] = 0;
        displayNormals[2] = 0;
        displayNormals[3] = 0;
        nodeBuffers.reset(new JNodeBuffers);
        edgeBuffers.reset(new JEdgeBuffers);
        faceBuffers.reset(new JFaceBuffers);
    }

    static int setNodeScalarFieldColor( const JMeshPtr &m, const vector<double> &c);
    static int setNodeScalarFieldColor( const JMeshPtr &m, const string &name);

    void setAlpha( double a) {
        alpha = a;
    }

    void setTransparency( bool t, int m = 0)
    {
        transparent = t;
        transparencyMethod = m;
    }

    void setRenderStyle(int s)
    {
        renderStyle = s;
    }

    int  getRenderStyle() const
    {
        return renderStyle;
    }

    void setSurfaceShade(int s)
    {
        surfShade = s;
    }

    int getSurfaceShade() const
    {

        return surfShade;
    }
    void setUniformColor( const JMeshPtr &mesh, int entity, JColor &c);

    int    renderStyle;
    int    surfShade = JRender::FLAT_SHADE;
    bool   faceCulling;
    bool   useNormals;
    bool   useMaterial;
    bool   transparent = 0;
    int    transparencyMethod = 0;
    JColor surfDefaultColor;

    int   pickableEntity = -1;
    bool  display = 1;
    bool  displayEntity[4];
    bool  displayIDs[4];
    bool  displayNormals[4];
    double alpha = 1.0;
    boost::shared_ptr<JNodeBuffers>   nodeBuffers;
    boost::shared_ptr<JEdgeBuffers>   edgeBuffers;
    boost::shared_ptr<JFaceBuffers>   faceBuffers;

    int   displayAll( int entity, bool val);
};

typedef boost::shared_ptr<JMeshRender> JMeshRenderPtr;
#include "EntityColor.hpp"
