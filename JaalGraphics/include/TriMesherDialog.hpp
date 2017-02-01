#pragma once

#include "Ui_TriMesherDialog.hpp"

#include <QDialog>
#include <QFileDialog>
#include <QColorDialog>

#include "MeshViewer.hpp"
#include "MeshPaver.hpp"
#include "MeshImporter.hpp"
#include "MeshExporter.hpp"
#include "DelaunayMesh.hpp"
#include "IsotropicMesh.hpp"

#include "TriMeshCleanupDialog.hpp"
#include "DelaunayMesherDialog.hpp"
#include "IsotropicTriMeshDialog.hpp"
#include "InstantMeshDialog.hpp"
#include "SurfaceReconstructionDialog.hpp"
#include "ContourEditingDialog.hpp"
#include "MeshHolesFillDialog.hpp"

class JTriDelaunayViewer : public JViewComponent
{
    static GLUquadricObj *diskObj;

public:
    JTriDelaunayViewer() {
        mesh = nullptr;
        kdt  = nullptr;
        displayCircumCircles = 0;
        displayConvexHull    = 0;
        displayMedialAxis    = 0;
    }

    void setMesh( const JMeshPtr &m ) {
        mesh = m;
    }

    void setConstraintMesh( const JMeshPtr &m) {
        kdt = m;
    }

    void drawCircumCircles(bool v)   {
        displayCircumCircles = v;
    }
    void drawConvexHull(bool v)      {
        displayConvexHull = v;
    }
    void drawConstraintMesh(bool v)  {
        displayConvexHull = v;
    }
    void drawMedialAxis(int v )      {
        displayMedialAxis = v;
    }

    void setMedial( const vector<Point3D> &p) {
        medialPoints = p;
    }
    void setConvexHull( const JEdgeSequence &e) {
        convexHull = e;
    }

    void setColor( const JColor &c);

    void draw();
private:
    JMeshPtr mesh;
    JMeshPtr kdt;
    bool  displayCircumCircles, displayConvexHull, displayMedialAxis;
    JEdgeSequence convexHull;

    vector<JCircle>  circumcircles;
    vector<Point3D> medialPoints;
    void setCircles();
    void genAdfrontMesh();
};

/////////////////////////////////////////////////////////////////////////////////////

class JTriMesherDialog : public QDialog, public Ui::TriMesherDialog {
    Q_OBJECT

public:
    JTriMesherDialog( QWidget *parent = 0);
    ~JTriMesherDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

protected:
    bool event( QEvent *e);
    void keyPressEvent( QKeyEvent *e);
    void showEvent( QShowEvent *event);

private slots:
    void openCleanupDialog();
    void openIsotropicDialog();
    void openInstantMeshDialog();
    void openQualityDelaunayDialog();
    void openSurfReconstructionDialog();
    void openContourEditingDialog();
    void openHolesFillDialog();

    void rejectNewMesh();
    void genNewMesh();
    void advancingFront();
    void getIntrinsicDelaunayMesh();

    void rejectMesh();
    void closeDialog();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    JMeshPtr cdt;      // Constrained Delaunay mesh;
    JMeshPtr chull;    // Convex hull;
    JMeshPtr medialAxis; // Medial Axis
    JMeshPtr adfront; // Medial Axis

    double    creaseAngle;
    JDelaunayMesh2D     trimesher;

    boost::scoped_ptr<JIsotropicTriMeshDialog> isotropicDialog;
    boost::scoped_ptr<JTrimeshCleanupDialog> cleanupDialog;
    boost::scoped_ptr<JDelaunayMesherDialog> delaunayDialog;
    boost::scoped_ptr<JInstantMeshDialog> instantMeshDialog;
    boost::scoped_ptr<JSurfaceReconstructionDialog> surfReconDialog;
    boost::scoped_ptr<JContourEditingDialog> contourEditingDialog;
    boost::scoped_ptr<JMeshHolesFillDialog> holesFillDialog;

    boost::scoped_ptr<JEdgeColor>  delEdgeColor;
    boost::scoped_ptr<JTriDelaunayViewer> delaunayViewer;

    void setMesh( const JMeshPtr &m);
    void init();
    void getConvexHull();
    void getMedialAxis();
    void makeConnections();
};

