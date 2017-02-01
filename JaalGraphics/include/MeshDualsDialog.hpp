#pragma once

#include <QDialog>

#include "Ui_MeshDualsDialog.hpp"

#include "QuadMeshDualsDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshDualViewer.hpp"
#include "QuadDual.hpp"

class JChordCellColor : public JCellColor {
public:
    string getName() const {
        return "ChordCell";
    }

    bool isPerCell() const {
        return 0;
    }
    void setChordID( int id );

    int assign(const JCellPtr &c);
    int  operator() (const JCellPtr &c) {
        return assign(c);
    }

private:
    int chordID;
    map<int, JColor> mapColor;
};

///////////////////////////////////////////////////////////////////////////////

class JSheetColor : public JFaceColor {
public:
    string getName() const {
        return "SheetColor";
    }

    bool isPerFace() const {
        return 1;
    }

    void setSheetID( int id );

    JColor getColor(int id ) {
        return mapColor[id];
    }

    int assign( const JFacePtr &face );

    int  operator() (const JFacePtr &f) {
        return assign(f);
    }
private:
    int chordID;
    map<int, JColor> mapColor;
};

///////////////////////////////////////////////////////////////////////////////


class JMeshDualsDialog : public QDialog, public Ui::MeshDualsDialog {
    Q_OBJECT

public:
    JMeshDualsDialog( QWidget *parent = 0);

    ~JMeshDualsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:

    void displayDuals();
    void pickSeeds();
    void nextSeeds();
    void selectDual();
    void getSomeDual();
    void modifyDual();

    void getStyle();
    void checkDisplay();

    void getAllChords();
    void getAllHexSheets();

    void clearChords();
    void resetPrimalMesh();
    void reject();

    void keyPressEvent( QKeyEvent *e);

private:
    JMeshPtr mesh;
    JMeshViewerPtr meshViewer;

    JaalViewer *viewManager;
    boost::shared_ptr<JMeshEntityPicker> entityPicker;

    string dualname;

    size_t nextEdgeID, nextFaceID, nextCellID;

    int  chordDim;    // 2 for quadmesh: 3 for hexmesh
    JQuadDualPtr qdual;
    JHexDualPtr  hdual;
    vector<JQuadChordPtr> qChords;
    vector<JDualSheetPtr> dSheets;
    vector<Tube> chordTubes;

    vector<JQuadChordPtr> quadChords;

    void init();
    void makeConnections();
    void setMesh( const JMeshPtr &m);
    void getQuadChords( JEdgeSequence &s);
    void getHexChords( JFaceSequence &s);
    void getHexSheets( JEdgeSequence &s);
    void getHexChords( JCellSequence &s);
    int  getSheaves(JEdgePtr &e);
    void getAllQuadChords();
    void getAllHexChords();
};

