#pragma once

#include "JaalViewer.hpp"
#include "DrawMeshEntity.hpp"
#include "EntityColor.hpp"

////////////////////////////////////////////////////////////////////////////////
class JNodePickColor : public JNodeColor {
public:
    string getName() const {
        return "PickNode";
    }

    bool isPerNode() const {
        return 0;
    }

    JNodePickColor();

    int  assign(const JNodePtr &v);

    int  operator() ( const JNodePtr &v) {
        return assign(v);
    }
};
////////////////////////////////////////////////////////////////////////////////

class JEdgePickColor : public JEdgeColor {
public:
    string getName() const {
        return "PickEdgeColor";
    }
    bool isPerEdge() const {
        return 0;
    }

    JEdgePickColor();

    int  assign(const JEdgePtr &);

    int  operator() (const JEdgePtr &e) {
        return assign(e);
    }
};
////////////////////////////////////////////////////////////////////////////////

class JFacePickColor : public JFaceColor {
public:
    string getName() const {
        return "PickFaceColor";
    }
    bool isPerFace() const {
        return 0;
    }

    JFacePickColor();

    int  assign( const JFacePtr &);

    int  operator() (const JFacePtr &f) {
        return assign(f);
    }
};
////////////////////////////////////////////////////////////////////////////////

class JCellPickColor : public JCellColor {
public:
    string getName() const {
        return "PickCellColor";
    }
    bool isPerCell() const {
        return 0;
    }

    JCellPickColor();

    int  assign( const JCellPtr &);
    int  operator() (const JCellPtr &c) {
        return assign(c);
    }
private:
    JFacePickColor *faceColor;
};

////////////////////////////////////////////////////////////////////////////////

class JMeshEntityPicker {
public:
    static const int SINGLE_SELECTION   = 1;
    static const int MULTIPLE_SELECTION = 2;

    JMeshEntityPicker();

    void setViewManager( JaalViewer *v );

    void setMesh( const JMeshPtr &m)
    {
        mesh = m;
    }

    void setActive( bool v )
    {
        active = v;
    }

    void setOp( int v )
    {
        pickOp = v;
    }

    void setMode( int mode )
    {
        selection_mode = mode;
    }

    int  getMode() const
    {
        return selection_mode;
    }

    void setHighlightColor( const JColor &c, int e);

    int  select_entity();

    void draw();

    void clearAll();

    JNodeSequence getPickedNodes() const;
    JEdgeSequence getPickedEdges() const;
    JFaceSequence getPickedFaces() const;
    JCellSequence getPickedCells() const;

    size_t getSize( int e ) const
    {
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

private:
    JaalViewer  *viewManager;
    boost::shared_ptr<JMeshViewer> meshViewer;
    JMeshPtr mesh;

    int    pickOp;
    int    selection_mode;
    int    active;
    int    pickWidth;

    JNodePtr currPickedNode;
    JEdgePtr currPickedEdge;
    JFacePtr currPickedFace;
    JCellPtr currPickedCell;

    set<int> picked_nodes, picked_edges, picked_faces, picked_cells;

    boost::scoped_ptr<JNodePickColor> pickNodeColor;
    boost::scoped_ptr<JEdgePickColor> pickEdgeColor;
    boost::scoped_ptr<JFacePickColor> pickFaceColor;
    boost::scoped_ptr<JCellPickColor> pickCellColor;

    QRect selectRect;

    void append( int id )
    {

        switch(id) {
        case 0:
            picked_nodes.insert(id);
            break;
        case 1:
            picked_edges.insert(id);
            break;
        case 2:
            picked_faces.insert(id);
            break;
        case 3:
            picked_cells.insert(id);
            break;
        }
    }

    void revert_colors(int );
    void init();
    void drawPickRegion();
};

typedef boost::shared_ptr<JMeshEntityPicker> JMeshEntityPickerPtr;

