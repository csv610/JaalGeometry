#pragma once

#include "Ui_GenSimpleShapeDialog.hpp"
#include "LSystemDialog.hpp"
#include "CurveGenDialog.hpp"
#include "PolygonDialog.hpp"
#include "PolyhedraDialog.hpp"
#include "KnotsDialog.hpp"
#include "ImageContoursDialog.hpp"

#include <QDialog>
#include "MeshViewer.hpp"

class JGenSimpleShapeDialog : public QDialog, public Ui::GenSimpleShapeDialog {
    Q_OBJECT

public:
    static const int  PLANE = 1;
    static const int  DISC  = 2;
    static const int  CUBE  = 3;
    static const int  CYLINDER  = 4;
    static const int  CAPPED_CYLINDER  = 5;
    static const int  CONE  = 6;
    static const int  CAPPED_CONE  = 7;
    static const int  PLATONIC_SOLID  = 8;
    static const int  SPHERE  = 9;
    static const int  SUBDIVIDED_SPHERE  = 10;
    static const int  HELIX  = 11;
    static const int  TORUS  = 12;
    static const int  KLEIN_BOTTLE  = 13;
    static const int  KNOT  = 14;
    static const int  RHOMBIC_TRICONTAHEDRON  = 15;
    static const int  RHOMBIC_DODECAHEDRON  = 16;

    JGenSimpleShapeDialog( QWidget *parent = 0);

    ~JGenSimpleShapeDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    JMeshPtr getShape( int i ) ;

private slots:
    void init();
    void genShape() ;
    void getPolygon();

    void openKnotsDialog();
    void openCurveDialog();
    void openLSystemDialog();
    void openPolygonDialog();
    void openPolyhedraDialog();
    void openImageEdgesDialog();

    void acceptShape() ;
    void cancelShape() ;

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr orgmesh, newmesh;

    boost::scoped_ptr<JKnotsDialog>      knotsDialog;
    boost::scoped_ptr<JCurveGenDialog>   curveDialog;
    boost::scoped_ptr<LSystemDialog>     lsystemDialog;
    boost::scoped_ptr<JPolygonDialog>    polygonDialog;
    boost::scoped_ptr<JPolyhedraDialog>  polyhedraDialog;
    boost::scoped_ptr<JImageContoursDialog>  imageEdgesDialog;

    void setDefault();
    void makeConnections();

    int genHelix();
    int genTorus();
    int genKleinBottle();
    int genKnot();
    int genSphere();
};
