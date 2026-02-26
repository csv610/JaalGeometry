#pragma once

#include "Mesh.hpp"
#include "MeshTopology.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

class JHexSheet;
class JHexChord;
class JHexDual;

typedef boost::shared_ptr<JHexSheet>  JHexSheetPtr;
typedef boost::shared_ptr<JHexChord>  JHexChordPtr;
typedef boost::shared_ptr<JHexDual>   JHexDualPtr;

struct JHexSheet {

    static JHexSheetPtr joinSheets( const JDualSheetPtr &a, const JDualSheetPtr &b);

    typedef std::pair<JCellPtr, JEdgePtr> Key;

    JHexSheet()
    {
        clear();
    }

    void clear()
    {
        ring              = 0;
        self_touching     = 0;
        start_at_boundary = 0;
        end_at_boundary   = 0;
        self_intersecting = 0;
        faces.clear();
        parEdges.clear();
        dElements.clear();
    }

    bool isSimple() const
    {
        if( self_intersecting || self_touching ) return 0;
        return 1;
    }

    bool hasKey( const Key &k) const
    {
        if( keys.find(k) == keys.end() ) return 0;
        return 1;
    }

    void addKey( Key &k)
    {
        keys.insert(k);
    }


    bool hasFace( const JFacePtr &f) const
    {
        if( find( faces.begin(), faces.end(), f ) == faces.end() ) return 0;
        return 1;
    }

    bool hasCell( const JCellPtr &c) const
    {
        size_t nSize = dElements.size();
        for( size_t i = 0; i < nSize; i++)
            if( dElements[i].cell == c ) return 1;
        return 0;
    }


    void shrink();    // make the parallel edges smaller.

    void  addFace( const JFacePtr &f)
    {
        faces.insert(f);
    }

    void  addElements(const JCellPtr &cell, const JEdgePtr &parEdge0, const JEdgePtr &parEdge1, const JEdgePtr &parEdge2, const JEdgePtr &parEdge3 )
    {
        parEdges.insert( parEdge0 );
        parEdges.insert( parEdge1 );
        parEdges.insert( parEdge2 );
        parEdges.insert( parEdge3 );

        JNodeSequence snodes(4);
        if( !parEdge0->hasAttribute("Steiner") ) {
            snodes[0] = JNodeGeometry::getMidNode( parEdge0->getNodeAt(0), parEdge0->getNodeAt(1));
            parEdge0->setAttribute("Steiner", snodes[0] );
        }
        parEdge0->getAttribute("Steiner", snodes[0] );

        if( !parEdge1->hasAttribute("Steiner") ) {
            snodes[1] = JNodeGeometry::getMidNode( parEdge1->getNodeAt(0), parEdge1->getNodeAt(1));
            parEdge1->setAttribute("Steiner", snodes[1] );
        }
        parEdge1->getAttribute("Steiner", snodes[1] );

        if( !parEdge2->hasAttribute("Steiner") ) {
            snodes[2] = JNodeGeometry::getMidNode( parEdge2->getNodeAt(0), parEdge2->getNodeAt(1));
            parEdge2->setAttribute("Steiner", snodes[2] );
        }
        parEdge2->getAttribute("Steiner", snodes[2] );

        if( !parEdge3->hasAttribute("Steiner") ) {
            snodes[3] = JNodeGeometry::getMidNode( parEdge3->getNodeAt(0), parEdge3->getNodeAt(1));
            parEdge3->setAttribute("Steiner", snodes[3] );
        }
        parEdge3->getAttribute("Steiner", snodes[3] );

        JFacePtr sface = Quadrilateral::newObject( snodes);
        dElements.push_back( Element(cell, sface, parEdge0));
    }

    void  getAffectedNeighbors( JCellSequence &c);

    void  getSheetFaces( JFaceSequence &sheetfaces) const
    {
        sheetfaces.clear();
        for( size_t i = 0; i < dElements.size(); i++)
            sheetfaces.push_back( dElements[i].genFace );
    }

    void  getEdges( JEdgeSequence &e) const
    {
        e.clear();
        JEdgeSet::iterator it;
        for( it = parEdges.begin(); it != parEdges.end(); ++it)
            e.push_back(*it);
    }

    void  getFaces( JFaceSequence &f) const
    {
        f.clear();
        JFaceSet::iterator it;
        for( it = faces.begin(); it != faces.end(); ++it)
            f.push_back(*it);

    }

    void  getCells( JCellSequence &cells) const
    {
        cells.clear();
        for( size_t i = 0; i < dElements.size(); i++) {
            JCellPtr c = dElements[i].cell;
            if( c ) cells.push_back( c );
        }
    }

    void Delete() { }

    bool  self_intersecting, self_touching;
    bool  ring, start_at_boundary, end_at_boundary;
private:
    struct Element {
        Element( JCellPtr c, JFacePtr f, JEdgePtr e)
        {
            cell = c;
            genFace  = f;
            parEdge  = e;
        }

        void clear()
        {
            cout << "EXIT" << endl;
            exit(0);
            /*
                        cell = nullptr;
                        parEdge = nullptr;
                        if( genFace ) {
                            delete genFace->getNodeAt(0);
                            delete genFace->getNodeAt(1);
                            delete genFace->getNodeAt(2);
                            delete genFace->getNodeAt(3);
                            delete genFace;
                        }
                        genFace = nullptr;
            */
        }
        JCellPtr cell;
        JFacePtr genFace;   // Sheet created within the cell.
        JEdgePtr parEdge;   // Parallel Edge that created this Dual sheet
    };
    std::set<Key> keys;
    void  expand();
    JEdgeSet parEdges;  // which mesh edges are in the dual. 4 per cell.
    JFaceSet faces;     // which mesh faces are in the dual: 2 per cell.
    vector<Element> dElements;
};

///////////////////////////////////////////////////////////////////////////////

class JHexChord;
typedef boost::shared_ptr<JHexChord> JHexChordPtr;

class JHexChord : public JDualChord {

public:
    struct Element {
        Element()
        {
            parEdge1 = nullptr;
            parEdge2 = nullptr;
            face  = nullptr;
            cell  = nullptr;
            segments[0] = nullptr;
            segments[1] = nullptr;
        }

        // Chord through Quadmesh.
        Element( JEdgePtr e, JFacePtr f)
        {
            parEdge1 = e;
            parEdge2 = nullptr;
            face  = f;
            cell  = nullptr;
            segments[0] = nullptr;
            segments[1] = nullptr;
        }

        // Chord through hexmesh.
        Element(JFacePtr f, JEdgePtr e1, JEdgePtr e2, JCellPtr c)
        {
            parEdge1 = e1;
            parEdge2 = e2;
            face  = f;
            cell  = c;
            segments[0] = nullptr;
            segments[1] = nullptr;
        }

        JEdgePtr parEdge1, parEdge2;     // Parallel Edge that created this Dual Chord
        JFacePtr face;
        JCellPtr cell;
        JEdgePtr segments[2];  // Two segments generated by each face/cell.
    };
public:
    int initialize();

    void reverse();

    JFacePtr getSeedFace() const
    {
        return seedFace;
    }

    void shrink_parallel_edges();
    int  delete_parallel_edges();
    void delete_dual_segments();

    void shrinkColumn();

    void append(const JDualChordPtr &other);

    void addIntersection( JFacePtr face)
    {
        selfIntersectingFaces.push_back(face);
        self_intersecting = 1;
    }

    void addIntersection( JCellPtr cell)
    {
        selfIntersectingCells.push_back(cell);
        self_intersecting = 1;
    }

    void addElements( JEdgePtr edge, JFacePtr face = nullptr);
    void addElements( JFacePtr face, JEdgePtr edge, JCellPtr cell = nullptr);

    void get_parallel_nodes(JNodeSequence &seq) const;

    void getCells( JCellSequence &seq) const;


    void getIntersectingCells( JCellSequence &seq) const
    {
        seq = selfIntersectingCells;
    }

    void getTouchingEdges( JEdgeSequence &seq) const
    {
        seq = selfTouchingEdges;
    }

    void finalize()
    {
        chordSegments.clear();
        size_t nSize = dElements.size();
        for( size_t i = 0; i < nSize; i++) {
            JEdgePtr e1 = dElements[i].segments[0];
            JEdgePtr e2 = dElements[i].segments[1];
            if( e1 ) chordSegments.push_back( e1 );
            if( e2 ) chordSegments.push_back( e2 );
        }
    }

    bool isRemovable();

    int diceQuads(int npieces = 3);
    int remove_self_intersecting_cycles();
    int remove_self_intersecting_cycle(JFacePtr f);

    bool hasEdge( const JEdgePtr e) const;
    bool hasFace( const JFacePtr f) const;
    bool hasCell( const JCellPtr c) const;

    void print();

    JNodeSequence collapseNodesPair;
    JCellSequence selfIntersectingCells;
    JFaceSequence selfIntersectingFaces;

    bool isTouching() const;

    bool isIntersecting() const
    {
        if( selfIntersectingFaces.empty() ) return 0;
        return 0;
    }

protected:
    int topDim;
    bool start_at_boundary, end_at_boundary;
    bool ring, self_intersecting, self_touching, simple;

    JMeshPtr mesh;
    JEdgePtr seedEdge;
    JFacePtr seedFace;

    vector<Element> dElements;

    JNodeSequence selfTouchingNodes;
    JEdgeSequence selfTouchingEdges;
    JFaceSequence selfTouchingFaces;

    JEdgeSequence chordSegments;

    void expand();
    void isRemovable2D();
    void isRemovable3D();
    int  isMergeable( JNodePtr v1, JNodePtr v2) const;
};

////////////////////////////////////////////////////////////////////////////////////

class JHexSheaf {
public:
    JHexSheaf( JMeshPtr m = nullptr )
    {
        mesh = m;
    }
    void setChords( JDualChordPtr &c1, JDualChordPtr &c2, JMeshPtr m = nullptr)
    {
        chord1 = c1;
        chord2 = c2;
        mesh   = m;
    }

    int  swapEdge(JNodePtr newnode);

private:
    JMeshPtr mesh;
    JDualChordPtr  chord1, chord2;
    bool isTight();
};

///////////////////////////////////////////////////////////////////////////////

struct JMeshDual { };

///////////////////////////////////////////////////////////////////////////////
class  JHexDual : public JMeshDual {
public:
    static void joinSheets( const JDualSheetPtr &a, const JDualSheetPtr &b, JDualSheetPtr &c);
    static void joinChords( const JDualChordPtr &a, const JDualChordPtr &b, JDualChordPtr &c);

    JHexDual()
    {
        mesh = nullptr;
    }

    void setMesh( const JMeshPtr &m )
    {
        mesh = m;
        initialize();
    }

    // Dice the given sheet ntimes. All the cells in it will be refined.
    int diceSheet(JDualSheetPtr s, int n = 1);

    // Delete a sheet from the mesh. Presently only simple sheets are removed..
    int removeColumn( JDualChordPtr c);

    // Delete a sheet from the mesh. Presently only simple sheets are removed..
    int removeSheet( JDualSheetPtr &s);

    JDualSheetPtr getSheet(JEdgePtr e);

    // Get some specific sheet starting from a given face ...
    JDualChordPtr getChord(JFacePtr f, JEdgePtr e = nullptr);
    int getChords(JEdgePtr s, vector<JDualChordPtr> &chords);
    int getChords(JFacePtr s, vector<JDualChordPtr> &chords);
    int getChords(JCellPtr s, vector<JDualChordPtr> &chords);
    int getChords(JFaceSequence &s, vector<JDualChordPtr> &chords);
    int getAllChords( vector<JDualChordPtr> &c);

    // Start a sheet from  a given face and give its two sheets. If  a sheet starts from
    // the boundary, then it must end with the boundary.
    int getSheets(JFacePtr f, vector<JDualSheetPtr> &sheets);

    // Start from a given cell, and give its three sheets ....
    int getSheets(JCellPtr h, vector<JDualSheetPtr> &sheets);

    void getSheets(JEdgeSequence &eseq, vector<JDualSheetPtr> &sheets);

    // Get all the sheets in the mesh. Generally, we call sheets one-by-one and process
    // them.
    int getAllSheets( vector<JDualSheetPtr> &sheets);

    // Get the column from the fiven face. Two sheets will pass through a column.
    int getColumn(JFacePtr f, vector<JDualSheetPtr> &sheet1, vector<JDualSheetPtr> &sheet2);

    int getColumn(JFacePtr f, JCellSequence &cells);

private:
    JMeshPtr mesh;
    void initialize();
    JEdgeSequence par_edges;   // Parallel edges
    void expandChord(const JCellPtr &c, const JFacePtr &f, const JEdgePtr &e, JDualChordPtr &d);
    void expandColumn(const JCellPtr &c, const JFacePtr &f, JCellSequence &cells);
    void expandSheet(const JCellPtr &c, const JFacePtr &f, const JEdgePtr &e1, JDualSheetPtr &d);
    void expandSheet(const JCellPtr &c, const JEdgePtr &e, JDualSheetPtr &d);
};
