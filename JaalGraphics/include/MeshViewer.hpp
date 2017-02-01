#pragma once

#include <qapplication.h>
#include <QGLViewer/qglviewer.h>
#include <QtWidgets>
#include <QMessageBox>
#include <QColorDialog>

#include <string>
#include <vector>

#ifdef HAVE_LOG4
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/helpers/exception.h>
#endif

#include "JaalHeaders.hpp"

#include <GL/gl.h>
#include <GL/gle.h>
#include <GL/glu.h>
#include "Lights.hpp"

using namespace std;
using namespace qglviewer;
using namespace Jaal;

class JMeshViewer;

#include "EntityColor.hpp"
#include "JaalViewer.hpp"
#include "FontsManager.hpp"
#include "MeshEntityPicker.hpp"
#include "MeshEntityAttribListDialog.hpp"
#include "DrawMeshEntity.hpp"

////////////////////////////////////////////////////////////////////////////

/*
struct JTransformNode
{
    JTransformNode ()
    {
       scale[0] = 0.0;
       scale[1] = 0.0;
       scale[2] = 0.0;
       rotate[0] = 0.0;
       rotate[1] = 0.0;
       rotate[2] = 0.0;
       translate[0] = 0.0;
       translate[1] = 0.0;
       translate[2] = 0.0;
    }
    Point3D scale;
    Point3D rotate;
    Point3D translate;
};
typedef std::unique_ptr<JTransformNode> JTransformNodePtr;

struct JMeshEnclosure
{
   int  type = 0;
   JBoundingBox  aabox;
   Hexahedron    minBox;
   JSphere       minSphere;
   JColor        color;
};

typedef std::unique_ptr<JMeshEnclosure> JMeshEnclosurePtr;
*/

class JMeshViewer : public JViewComponent {
public:

    // We can use four ways to get hidden line removal effect.
    static const int AXIS_ALIGNED_BOUNDING_BOX = 0;
    static const int MINIMUM_BOUNDING_BOX      = 1;
    static const int MINIMUM_SPHERE            = 2;
    static const int CONVEX_HULL               = 3;

    static JViewComponentPtr registerComponent(JaalViewer *root);

    explicit JMeshViewer( JaalViewer *p);

    ~JMeshViewer();

    boost::shared_ptr<JMeshEntityPicker>  getEntityPicker() const
    {
        return entityPicker;
    }

    size_t getNumVisible( int entity );
    size_t getNumVisible( const JMeshPtr &m, int entity );

    JNodeDraw *getNodeDraw() const
    {
        return nodeDraw.get();
    }

    JEdgeDraw *getEdgeDraw() const
    {
        return edgeDraw.get();
    }

    JFaceDraw *getFaceDraw() const
    {
        return faceDraw.get();
    }

    JCellDraw *getCellDraw() const
    {
        return cellDraw.get();
    }

    void setRenderMode(int m);

    void fastDraw();
    void draw();
    void refreshDisplay();

/*
    void displayAll(bool val)
    {
        displayAll(0,val);
        displayAll(1,val);
        displayAll(2,val);
        displayAll(3,val);
    }
    void displayAll( int entity, bool val);
*/

    void displayAll(const JMeshPtr &mesh, bool val)
    {
        displayAll(mesh, 0, val);
        displayAll(mesh, 1, val);
        displayAll(mesh, 2, val);
        displayAll(mesh, 3, val);
    }
    void displayAll( const JMeshPtr &mesh, int entity, bool val);

    // Update all the buffer of every mesh object enlisted.
    void initBuffers( const JMeshPtr &mesh);

    void updateBuffers();
    void updateBuffers( const JMeshPtr &mesh);
    void updateGeometryBuffers( const JMeshPtr &mesh);
    void updateTopologyBuffers( const JMeshPtr &mesh);

    void displayEnclosure(bool val, int t = AXIS_ALIGNED_BOUNDING_BOX);

    void displayNormals( int entity, bool val );

    int   getSize() const
    {
        return meshdb.size();
    }

    JMeshPtr  loadMesh( const string &s);

    int  addObject( const JMeshPtr &m);

    // Return the mesh by name ...
    JMeshPtr getMesh(const string &s)
    {
        for( const JMeshPtr &mesh: meshdb)
            if(mesh->getName() == s) return mesh;
        JMeshPtr nullPtr;
        return nullPtr;
    }

    // Return the mesh by id ...
    JMeshPtr getMesh(int id) const
    {
        if( !meshdb.empty() ) return meshdb[id];
        JMeshPtr nullPtr;
        return nullPtr;

    }

    void setCurrentMesh( const JMeshPtr &m) {
        if( find( meshdb.begin(), meshdb.end(), m) != meshdb.end() )
            currMesh = m;
    }

    JMeshPtr getCurrentMesh() const {
        return currMesh;
    }

    int  removeObject(const JMeshPtr &m);

    void lookAt( const JMeshPtr &m);
    void lookAt( const JNodePtr &v);

    void alignAlong(const JMeshPtr &mesh, const Vec3F &srcVec, const Vec3F &dstVec);
    void alignAlong(const JMeshPtr &mesh, const JEdgePtr &edge, int axis = 0);
    void alignAlong(const JMeshPtr &mesh, const JFacePtr &face, int axis = 2);

    void actionMouseEvent(int id);
    JBoundingBox  getAxisAlignedBox();

private:
    size_t currCounter;
    vector<JMeshPtr>  meshdb;
    JMeshPtr  currMesh;

    boost::shared_ptr<JMeshEntityPicker> entityPicker;
    boost::scoped_ptr<JNodeDraw>  nodeDraw;
    boost::scoped_ptr<JEdgeDraw>  edgeDraw;
    boost::scoped_ptr<JFaceDraw>  faceDraw;
    boost::scoped_ptr<JCellDraw>  cellDraw;

    bool draw_axis;
    int  edgemeshStyle, facemeshStyle;
    int  hiddenlinesMethod;
    int  enclosure_type;
    int  display_boundary;
    int  display_enclosure;
    int  display_dual_graph;
    int  display_primal_mesh;

    JMeshPtr readData( const string &f);

    void init();
    void animate();
    void drawWithNames();

    void initMesh( const JMeshPtr &m);

    void drawEnclosure( const JMeshPtr &m);
    void draw_ids();
    void setNodesGlyph( int glyph, double r = 0.0 );
    void setEdgesGlyph( int glyph, double r = 0.0 );
};

typedef boost::shared_ptr<JMeshViewer>  JMeshViewerPtr;

///////////////////////////////////////////////////////////////////////////////

class JSimpleMeshViewer : public JViewComponent {
public:
    explicit JSimpleMeshViewer( JaalViewer *p);

    ~JSimpleMeshViewer();

    void fastDraw();
    void draw();
    void refreshDisplay();

    void loadNewMesh( const string &s);
private:
    JMeshPtr  mesh;
    vector<JMeshPtr>  meshdb;
};
///////////////////////////////////////////////////////////////////////////////
