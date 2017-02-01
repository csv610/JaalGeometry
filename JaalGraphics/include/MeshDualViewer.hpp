#pragma once

#include <GL/gl.h>

#include "MeshDual.hpp"
#include "MeshViewer.hpp"
#include "EntityColor.hpp"

#include "JaalViewer.hpp"

using namespace Jaal;

class JaalViewer;

class JMeshDualViewer  : public JViewComponent {
public:
    typedef boost::shared_ptr<JMeshDualViewer> shared_ptr;

    explicit JMeshDualViewer(JaalViewer *p) {
        assert( p != nullptr);
        viewManager = p;
        init();
    }

    ~JMeshDualViewer() { }

    JaalViewer* getViewManager() const {
        return viewManager;
    }

    void setPrimalMesh(const JMeshPtr &m) {
        primalMesh = m;
    }

    JMeshPtr getPrimalMesh() {
        return primalMesh;
    }

    void setDualMesh(const JMeshPtr &) {
        /*
                  dualMesh = m;
                  dualMeshViewer->setNewMesh(dualMesh);
        */
    }

    JMeshViewerPtr getMeshViewer() const {
        return dualMeshViewer;
    }

    void displayGraph(bool v) {
        display_graph = v;
    }

    void displayDuals(bool v) {
        display_duals = v;
    }

    void displayParallelEdges(bool v) {
        parallel_edges = v;
    }

    // Display the chord either as simple line or tube..
    void setChordStyle( bool v ) {
        chordStyle = v;
//        refreshDisplay();
    }

    void displayBeads( bool v ) {
        chordBeads = v;
        refreshDisplay();
    }

    // Display one particular chord.
    void setChord(const JDualChordPtr &ch);
    void setChords( const vector<JDualChordPtr> &vch);

    void setSheet(JDualSheetPtr &ch) {
        sheets.resize(1);
        sheets[0] = ch;
    }

    void setSheets( const vector<JDualSheetPtr> &vch) {
        sheets = vch;
    }

    // Display the elements making chord or sheet with different color and
    // transparency so that user can visualize the dual effectively...
    void displayElements( bool v) {
        display_elements = v;
        refreshDisplay();
    }

    // Display self intersecting elements with different color.
    void displayIntersectingElements (int e, bool v ) {
        if( e == 2) display_intersecting_faces = v;
        if( e == 3) display_intersecting_cells = v;
//        refreshDisplay();
    }

    // Display self touching nodes. edges, and faces of a chord/sheet.
    void displayTouchingElements (int e, bool v ) {
        if( e == 0) display_touching_nodes = v;
        if( e == 1) display_touching_edges = v;
        if( e == 2) display_touching_faces = v;
//        refreshDisplay();
    }

    // Completely clear chords.
    void clearChords() {
        chords.clear();
        refreshDisplay();
    }

    // Completely clear sheets...
    void clearSheets() {
        for( size_t i = 0; i < sheets.size(); i++)
            sheets[i]->clear();
        sheets.clear();
        refreshDisplay();
    }

    void deleteAll();

    void draw();

    void refreshDisplay() {
        if( viewManager )
            viewManager->refreshDisplay();
    }

private:
    JMeshPtr  primalMesh, dualMesh;
    JaalViewer *viewManager;
    JMeshViewerPtr dualMeshViewer, primalMeshViewer;

    bool   chordBeads;
    bool   chordStyle;
    size_t currCounter;

    vector<Tube>  chordTubes;
    vector<JDualChordPtr>  chords;
    vector<JDualSheetPtr>  sheets;

    JFaceSequence sheetFaces; // get the Dual Sheet faces ...
    JEdgeSequence parEdges;   // get the parallel edges making the sheet ...

    JNodeSequence chordNodes;
    JEdgeSequence chordEdges;
    JFaceSequence chordFaces, intersectingFaces;
    JCellSequence chordCells, intersectingCells;

    boost::shared_ptr<JFaceColor> chordFaceColor;
    boost::shared_ptr<JFaceColor> sheetColor;

    int  curr_sheet_dir;
    bool display_graph;
    bool display_duals;
    bool display_chord_faces, display_chord_cells;
    bool parallel_edges;

    void resetAll( bool val );
    void activate( const JFacePtr &f);
    void activate( const JCellPtr &c);

    bool display_elements;
    bool display_intersecting_faces, display_intersecting_cells;
    bool display_touching_nodes, display_touching_edges, display_touching_faces;

    void init();
    void draw_sheet( const JDualSheetPtr &c, int id);
    void draw_chord( const JDualChordPtr &c, int id);
};

