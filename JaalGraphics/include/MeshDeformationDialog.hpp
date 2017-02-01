#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshDeformationDialog.hpp"

#include "MeshViewer.hpp"
#include "ImageViewer.hpp"
#include "LocallyInjectiveMap.hpp"
#include "MeshConstraintsDialog.hpp"
#include "StructuredMeshDialog.hpp"
#include "LocallyInjectiveMapParamsDialog.hpp"
#include "MeshGeometricQualityDialog.hpp"
#include "TetMesherDialog.hpp"
#include "MeshQuality.hpp"
#include "CurveShorteningFlowDialog.hpp"

class JMeshDeformationDialog;
class JMeshConstraintsDialog;
class JTetMesherDialog;

////////////////////////////////////////////////////////////////////////////////

class JMeshDeformViewer : public JViewComponent
{
public:
    static const int  AFFINE_TRANSLATE = 2;
    static const int  AFFINE_ROTATE    = 3;

    JMeshDeformViewer();

    void setArrows( bool v ) {
        arrows = v;
    }

    void setLasso( const vector<Point3D> &lp) {
        lassoPoints = lp;
    }

    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    /*
        void setBackgroundColor(const Color &bg) {
            backgroundColor[0] = bg[0];
            backgroundColor[1] = bg[1];
            backgroundColor[2] = bg[2];
            backgroundColor[3] = bg[3];
        }

        void setHandle( const JNodePtr &v ) {
            vertexHandle = v;
        }


        void setNodeRadius( double r ) {
            radius = r;
        }

        void setModifier(int b =0) {
            modifier = b;
        }

        void setTranslateVector( const Point3D &vtail, const Point3D &vhead) {
            modifier = AFFINE_TRANSLATE;
            translatePoints[0] = vtail;
            translatePoints[1] = vhead;
        }

        void setRotationArc( const Point3D &vfrom, const Point3D &vto) {
            modifier = AFFINE_ROTATE;
            rotationPoints[0] = vfrom;
            rotationPoints[1] = vto;
        }

        int setConstraint(const JNodePtr &v, bool val = 1);
        int setConstraints(const JMeshPtr &mesh);

        void clearLasso() { lassoPoints.clear(); }


        void deleteGroup(int id ) {
            nodeGroups.erase(id);
        }
    */

    void draw();
private:
    JMeshPtr  mesh;
    vector<Point3D> lassoPoints;
    bool   arrows;

    /*
        Color  srcColor, greenColor, dstColor, handleColor;
        Color  backgroundColor;
        double radius;
        int    modifier;
        JNodePtr vertexHandle;
        map<int,JNodeSequence>  nodeGroups;
        Point3D rotationPoints[2];
        Point3D translatePoints[2];
    */
};

////////////////////////////////////////////////////////////////////////////////

class JMeshDeformationDialog : public QDialog, public Ui::MeshDeformationDialog
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
    JMeshDeformationDialog( QWidget *parent = 0);
    ~JMeshDeformationDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

protected:
    void  mouseMoveEvent( QMouseEvent *e);
    void  mousePressEvent( QMouseEvent *e);
    void  mouseReleaseEvent( QMouseEvent *e);
    void  keyPressEvent( QKeyEvent *e);
    void  showEvent(QShowEvent *e);

private slots:
    void resetData();
    void checkDisplay();
    void getImageMesh();
    void openLIMDialog();
    void openCurveShorteningDialog();
    void openMeshGeomQualityDialog();

    void getNewConstraints();
    void loadSource();
    void loadTarget();
    void runSolver();
    void setHandleAsConstraint();
    void setFixedBoundary();
    void closeDialog();

private:
    JaalViewer   *viewManager;

    JMeshPtr mesh;
    int  entityDim;
    int  entityType[4];
    JMeshPtr simplicialMesh;
    JMeshViewerPtr  meshViewer;
    JImageViewerPtr imageViewer;
    JMeshEntityPickerPtr picker;
    bool initmesh;
    JNodeSet   fixedNodes;

    boost::shared_ptr<JLocallyInjectiveMap>   limDeformer;
    boost::shared_ptr<JMeshConstraintsDialog> meshConstraintsDialog;
    boost::shared_ptr<JMeshDeformViewer>      deformViewer;

    boost::scoped_ptr<JStructuredMeshDialog> structmeshDialog;
//    boost::scoped_ptr<JMeshRefine2DDialog>   refine2dDialog;
    boost::shared_ptr<JTetMesherDialog>      tetmesherDialog;
    boost::scoped_ptr<JLocallyInjectiveMapParamsDialog> limParamsDialog;
    boost::scoped_ptr<JMeshGeometricQualityDialog> meshGeomQualityDialog;
    boost::scoped_ptr<JCurveShorteningFlowDialog> curveShorteningDialog;

    boost::scoped_ptr<MyThread>  thread;

    vector<double> orgCoords;
    vector<size_t> l2g;

    // Node which remain fixed or move toward their specified targets...
    JNodeSet constraintNodes;

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
    void displaceConstraints( Point3D &by);
    void rotateConstraints( double angle);

    void init();
    void countNodes();
    void mesquiteSolver();
    void injectiveSolver();
    int  getInvertedElements();
    void getWorstQuality();
    void dst2src();
    void src2dst();
};

////////////////////////////////////////////////////////////////////////////////

