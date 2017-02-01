#pragma once

#ifndef JAAL_MESHVIEWER
#define JAAL_MESHVIEWER

#include <QtWidgets>
#include <QGLViewer/qglviewer.h>
#include <qapplication.h>

#include <string>
#include <vector>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/helpers/exception.h>

#include "JaalHeaders.hpp"

#include <GL/gl.h>
#include <GL/gle.h>
#include <GL/glu.h>
#include "Lights.hpp"

using namespace std;
using namespace qglviewer;
using namespace Jaal;

#include "EntityColor.hpp"
#include "JaalViewer.hpp"
#include "FontsManager.hpp"

#include "MeshEntityAttribListDialog.hpp"
#include "DrawMeshEntity.hpp"

class JMeshViewer;

class JMeshEntityPicker {
public:
    static const int SINGLE_SELECTION   = 1;
    static const int MULTIPLE_SELECTION = 2;

    JMeshEntityPicker();

    void setMeshViewer( JMeshViewer *m );

    void setStatus( bool v ) {
        picking_on = v;
    }

    void setOp( int v ) {
        pickOp = v;
    }

    void setMode( int mode ) {
        selection_mode = mode;
    }

    int  getMode() const {
        return selection_mode;
    }

    void setPickableEntity( int e = -1) {
        pick_entity = e;

        picking_on  = 0;
        if( e >= 0) picking_on  = 1;
    }

    int  getPickableEntity() const {
        return pick_entity;
    }

    void setPickedEntityColor( float *rgb, int e);

    void select_entity();

    void draw();

    void clearAll();

    JNodeSequence getPickedNodes() const;
    JEdgeSequence getPickedEdges() const;
    JFaceSequence getPickedFaces() const;
    JCellSequence getPickedCells() const;

    size_t getSize( int e ) const {
        switch(e) {
        case 0:
            return  picked_nodes.size();
        case 1:
            return  picked_edges.size();
        case 2:
            return  picked_faces.size();
        case 3:
            return  picked_cells.size();
        }
        return 0;
    }

    PickNodeColor *getPickNodeColor() const {
        return pickNodeColor;
    }
    PickEdgeColor *getPickEdgeColor() const {
        return pickEdgeColor;
    }
    PickFaceColor *getPickFaceColor() const {
        return pickFaceColor;
    }
    PickCellColor *getPickCellColor() const {
        return pickCellColor;
    }

private:
    JMeshViewer *meshViewer;
    Mesh *mesh;

    bool   picking_on;
    int    pick_entity;
    int    pickOp;
    int    selection_mode;

    Vertex *currPickedNode;
    Edge   *currPickedEdge;
    Face   *currPickedFace;
    Cell   *currPickedCell;

    set<int> picked_nodes, picked_edges, picked_faces, picked_cells;

    PickNodeColor *pickNodeColor;
    PickEdgeColor *pickEdgeColor;
    PickFaceColor *pickFaceColor;
    PickCellColor *pickCellColor;

    NodeColor *preNodeColor;
    EdgeColor *preEdgeColor;
    FaceColor *preFaceColor;
    CellColor *preCellColor;

    void init();
};

///////////////////////////////////////////////////////////////////////////////

class JMeshViewer : public ViewComponent
{
public:
    static log4cxx::LoggerPtr logger;
    static void init_logger();

    // How an edge mesh be rendered: (1) Simple wireframe (2) Hidden lines removed.
    static const int MESH_EDGES_LINE    = 0;
    static const int MESH_EDGES_HIDDENLINES_REMOVED  = 1;

    // We can use four ways to get hidden line removal effect.
    static const int HIDDENLINE_WITH_DEPTH_TEST      = 0;
    static const int HIDDENLINE_WITH_BACKFACES_CULL  = 1;
    static const int HIDDENLINE_WITH_FRONTLINES      = 2;
    static const int HIDDENLINE_WITH_STENCIL_BUFFER  = 3;

    static const int AXIS_ALIGNED_BOUNDING_BOX = 0;
    static const int MINIMUM_BOUNDING_BOX      = 1;
    static const int MINIMUM_SPHERE            = 2;
    static const int CONVEX_HULL               = 3;

    void alignAlong(Mesh *mesh, const Vec3D &srcVec, const Vec3D &dstVec);

    explicit JMeshViewer( JaalViewer *p);

    ~JMeshViewer();

    JaalViewer*  getViewManager() const {
        return viewManager;
    }

//   void addObserver( boost::any o) { observers.push_back(o); }

    JMeshEntityPicker *getEntityPicker() const {
        return entityPicker;
    }

    void addMeshFile( const string &f) {
        meshFileName = f;
    }

    Mesh *getMesh() const {
        return mesh;
    }

    void setNewMesh( Mesh *m ) {
        mesh = m;
        init_mesh();
        viewManager->updateGL();
    }

    void attach( JCurve *c );
    void detach( JCurve *c );

    size_t getNumVisible( int entity );

    JLights *getLights() {
        if(viewManager)
            return viewManager->getLights();
        return NULL;
    }

    Mesh *getDualGraph() const {
        return dualGraph;
    }

    void setDualGraph( Mesh *m ) {
        dualGraph = m;
        viewManager->updateGL();
    }

    DrawNode *getDrawNode() const {
        return drawNode;
    }

    DrawEdge *getDrawEdge() const {
        return drawEdge;
    }

    DrawFace *getDrawFace() const {
        return drawFace;
    }

    DrawCell *getDrawCell() const {
        return drawCell;
    }

    void draw();

    void refreshDisplay();

    void displayEntity( int entity, bool val) {
        display_entity[entity] = val;
    }

    bool isEnabled(int e) const {
        return display_entity[e];
    }

    void setCenter( const Point3D &p) {
        if( mesh == NULL ) return;
        size_t numnodes = mesh->getSize(0);
        for( size_t i = 0; i < numnodes; i++) {
            Vertex *v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            xyz[0] -= p[0];
            xyz[1] -= p[1];
            xyz[2] -= p[2];
            v->setXYZCoords(xyz);
        }
    }

    void setFaceMeshStyle(int type ) {
        facemeshStyle = type;
    }

    void setEdgeMeshStyle( int v ) {
        edgemeshStyle = v;
        viewManager->updateGL();
    }

    void setHiddenLinesMethod( int m ) {
        hiddenlinesMethod = m;
    }

    void  displayDualGraph( bool v ) {
        display_dual_graph = v;
    }

    void displayEnclosure(bool val, int t = AXIS_ALIGNED_BOUNDING_BOX) {
        display_enclosure = val;
        enclosure_type = t;

        if( val && mesh != NULL) {
            switch( enclosure_type )
            {
            case AXIS_ALIGNED_BOUNDING_BOX:
                box = mesh->getGeometry()->getBoundingBox();
                break;
            case MINIMUM_BOUNDING_BOX:
                if( minBox ) delete minBox;
                minBox = mesh->getGeometry()->getMinimumBox();
                break;
            }
        }
    }

    void displayNormals( int entity, bool val );

    void loadNewMesh( const string &s);
    void saveMesh( const string &s);

    void updateMouseReleasedEvent();
    void updateGeometry() ;

    void resetAll( bool v );
    void resetAll( int e, bool v );

    void alignAlong( Edge *edge, int axis, bool refresh);
    void alignAlong( Face *face, int axis, bool refresh);

private:
    string  meshFileName;
    int    enclosure_type;

    Mesh  *mesh,  *dualGraph;
    vector<Mesh*>  deletedmesh;
    vector<JCurve*> curves;
    // Observer design pattern: Which dialog boxes will like to
    // receive a notification, if certain changes take place
    //    vector<boost::any> observers;
    void notify();

    JLights    *lights;    // Controlled and held by viewmanager.
    JaalViewer *viewManager;
    JMeshEntityPicker *entityPicker;

    DrawNode  *drawNode;
    DrawEdge  *drawEdge;
    DrawFace  *drawFace;
    DrawCell  *drawCell;

    MeshNodeColor *meshNodeColor;
    MeshEdgeColor *meshEdgeColor;
    MeshFaceColor *meshFaceColor;
    MeshCellColor *meshCellColor;

    bool display_entity[4];

    size_t currCounter;

    bool   draw_axis;

    int  edgemeshStyle, facemeshStyle;
    int  hiddenlinesMethod;

    bool  display_normals[4];

    int  display_boundary, display_enclosure;
    int  display_dual_graph;

    void readData( const string &f);

    Jaal::BoundingBox  box;
    Hexahedron *minBox;

    void init();
    void animate();
    void drawWithNames();

    void init_mesh();

    void draw_ids();
    void draw_curves();

    void draw_nodes();
    void draw_edges();
    void draw_faces();
    void draw_cells();

    void draw_graph_nodes();
    void draw_graph_edges();

    void draw_hidden_lines();
    void draw_filled_faces();
    void draw_wired_faces();
    void draw_all_faces();

    void draw_primal_mesh();
    void draw_dual_graph();

    void draw_faces_normal();
    void draw_nodes_normal();

    void draw_boundary();
    void draw_enclosure();
    void draw_bounding_box();
    void draw_minimum_box();
};

///////////////////////////////////////////////////////////////////////////////

inline
void JMeshViewer :: displayNormals( int entity, bool val )
{
    if( mesh == NULL ) return;

    if( entity < 0 || entity > 2 ) return;

    display_normals[entity] = val;

    if( entity == 0 ) {
        mesh->getGeometry()->setNodesNormal();
        drawNode->displayNormals( val );
    }

    if( entity == 2 ) {
        mesh->getGeometry()->setFacesNormal();
        drawFace->displayNormals( val );
    }

    viewManager->updateGL();
}
///////////////////////////////////////////////////////////////////////////////

#endif
