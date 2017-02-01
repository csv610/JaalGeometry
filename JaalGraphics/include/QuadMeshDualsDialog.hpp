#pragma once

#include <QDialog>

#include "Ui_QuadMeshDualsDialog.hpp"

#include "MeshViewer.hpp"
#include "QuadDual.hpp"
#include "MeshOptBoundaryLayerDialog.hpp"

class JChordNodeColor : public JNodeColor {
public:
    string getName() const {
        return "ChordNode";
    }

    bool isPerNode() const {
        return 0;
    }
    void generateColors(int n);

    int assign(const JNodePtr &vertex);

    int  operator() (const JNodePtr &f) {
        return assign(f);
    }

private:
    map<int, JColor> mapColor;
};
///////////////////////////////////////////////////////////////////////////////

class JChordEdgeColor : public JEdgeColor {
public:
    string getName() const {
        return "ChordEdge";
    }

    bool isPerEdge() const {
        return 0;
    }

    void generateColors(int n);

    int assign( const JEdgePtr &edge);

    int  operator() (const JEdgePtr &e) {
        return assign(e);
    }

private:
    map<int, JColor> mapColor;
};
///////////////////////////////////////////////////////////////////////////////


class JChordFaceColor : public JFaceColor {
public:
    string getName() const {
        return "ChordFace";
    }

    bool isPerFace() const {
        return 0;
    }
    void generateColors( int n );

    int assign( const JFacePtr &face );

    int  operator() (const JFacePtr &f) {
        return assign(f);
    }
private:
    map<int, JColor> mapColor;
};
///////////////////////////////////////////////////////////////////////////////

class JQuadMeshDualsDialog : public QDialog, public Ui::QuadMeshDualsDialog {
    Q_OBJECT

public:
    JQuadMeshDualsDialog( QWidget *parent = 0);

    ~JQuadMeshDualsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh(const JMeshPtr &m);

protected:
    virtual void showEvent( QShowEvent *e);
    virtual void  mouseReleaseEvent( QMouseEvent *e);

private slots:

    void enablePicking();
    void nextSeeds();
    void selectDual();
    void getNewDual();
    void modifyDual();
    void diceAllChords();
    void displayCyclicChords();

    void getStyle();
    void checkDisplay();

    void getAllChords();
    void clearChords();
    void closeDialog();

    void keyPressEvent( QKeyEvent *e);
    void getMaxEdgeLength();
    void openBoundaryLayerDialog();

private:
    JMeshPtr mesh;
    JMeshViewerPtr meshViewer;
    JaalViewer *viewManager;
    boost::shared_ptr<JMeshEntityPicker> entityPicker;
    boost::shared_ptr<JChordFaceColor>   chordFaceColor;
    boost::shared_ptr<JChordEdgeColor>   chordEdgeColor;

    boost::scoped_ptr<JMeshOptBoundaryLayerDialog>  boundaryLayerDialog;

    string dualname;

    size_t nextEdgeID, nextFaceID;

    JQuadDualPtr qdual;
    vector<JQuadChordPtr> qChords;

    JMeshNonlinearOptimization nonlinearOpt;

    void setChordID(const JQuadChordPtr &chord, int id);

    void init();
    void initMesh();
    void makeConnections();
    void getEdgeChords(JEdgeSequence &eseq);
    void displayChords();
};

