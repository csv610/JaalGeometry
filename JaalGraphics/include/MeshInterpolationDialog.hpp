#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshInterpolationDialog.hpp"

#include "MeshViewer.hpp"
#include "ImageViewer.hpp"
#include "LocallyInjectiveMap.hpp"
#include "MeshConstraintsDialog.hpp"
#include "StructuredMeshDialog.hpp"
#include "MeshRefine2DDialog.hpp"
#include "LocallyInjectiveMapParamsDialog.hpp"
#include "MeshGeometricQualityDialog.hpp"
#include "TetMesherDialog.hpp"
#include "MeshQuality.hpp"
#include "MeshRefine2DDialog.hpp"
#include "DDG_MeshElasticDeformation.hpp"
#include "ObjectsListDialog.hpp"
#include "MeshInterpolation.hpp"
#include "HarmonicMap.hpp"

class JShapeInterpolationDialog;
class JMeshConstraintsDialog;

////////////////////////////////////////////////////////////////////////////////

class JMeshInterpolationDialog : public QDialog, public Ui::MeshInterpolationDialog
{
    Q_OBJECT

    friend class JMeshDeformViewer;

    struct MyThread : public QThread {
        void run() {
            deformer->solve();
            viewManager->refreshDisplay();
        }
        JaalViewer *viewManager;
        boost::shared_ptr<JLocallyInjectiveMap> deformer;
    };

public:
    JMeshInterpolationDialog( QWidget *parent = 0);
    ~JMeshInterpolationDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setSource( const JMeshPtr &m);
    void setTarget( const JMeshPtr &m);
    JMeshPtr getInterpolatedMesh(double t);

protected:
    void  mouseMoveEvent( QMouseEvent *e);
    void  mousePressEvent( QMouseEvent *e);
    void  mouseReleaseEvent( QMouseEvent *e);
    void  keyPressEvent( QKeyEvent *e);
    void  showEvent(QShowEvent *e);

private slots:
    void checkDisplay();
    void getImageMesh();
    void openLIMDialog();
    void openMeshGeomQualityDialog();
    void closeDialog();
    void loadSource();
    void loadTarget();
    void deformSolver();
    void getHausdorff();


private:
    JaalViewer   *viewManager;
    JMeshViewerPtr  meshViewer;
    JImageViewerPtr imageViewer;

private:
    bool  selectMesh[2];
    JMeshPtr sourceMesh, targetMesh, interpolatedMesh;
    JMeshEntityPickerPtr picker;

    boost::scoped_ptr<MyThread>  thread;
    boost::shared_ptr<JTetMesherDialog>       tetmesherDialog;
    boost::scoped_ptr<JObjectsListDialog>     meshListDialog;
    boost::shared_ptr<JLocallyInjectiveMap>   limDeformer;
    boost::scoped_ptr<DDG::JMeshElasticDeformation> elasticDeform;
    boost::scoped_ptr<JMeshGeometricQualityDialog> meshGeomQualityDialog;
    boost::scoped_ptr<JLocallyInjectiveMapParamsDialog> limParamsDialog;
    boost::scoped_ptr<JObjectsListDialog> meshlistDialog;
    boost::scoped_ptr<JMeshInterpolation> meshInterpolation;
    boost::scoped_ptr<JHarmonicMap>       harmonicMap;

    // Data required for moving  one node ..
    int   ncounter;
    bool  left_button_pressed;
    bool  constrains_moved;

    double  maxDist;
    JNodePtr moveVertex;
    JNodeSequence nodeNeighs;
    Point3D prevCoords, currCoords;
    JMeshNonlinearOptimization meshOpt;

    void   makeConnections();
    double getMaxDistance( const Point3D &p) const;
    void   setMaxTargetDistance();

    void initMesh();
    void genUVCoords();
    void initSolver();
    void genImageMesh();

    void assignColors();

    void init();
    void countNodes();
    void mesquiteSolver();
    void injectiveSolver();
    void getHarmonicMap();
    int  getInvertedElements( const JMeshPtr &m);
    void getWorstQuality();
};

////////////////////////////////////////////////////////////////////////////////

