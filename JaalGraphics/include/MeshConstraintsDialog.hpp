#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshConstraintsDialog.hpp"

#include "Mesh.hpp"
#include "PolygonDialog.hpp"
#include "MeshDeformationDialog.hpp"

class JMeshDeformViewer;
class JConstraintSet;
/////////////////////////////////////////////////////////////////////////////////////////////

class JMeshConstraintsViewer : public JViewComponent
{
public:
    JMeshConstraintsViewer();
    void init();

    void setCanvas( bool v ) {
        canvas = v;
    }

    int setConstraints(const JMeshPtr &mesh);
    int setConstraint(const JNodePtr &v, bool val = 1);

    void setLasso( const vector<Point3D> &lp) {
        lassoPoints = lp;
    }

    void clearAll() {
        lassoPoints.clear();
        nodeGroups.clear();
    }

    void clearLasso() {
        lassoPoints.clear();
    }

    void deleteGroup(int id ) {
        nodeGroups.erase(id);
    }

    void draw();
private:
    JMeshPtr   mesh;
    JColor redColor, greenColor;
    bool   canvas;
    vector<Point3D> lassoPoints;
    map<int,JNodeSequence>  nodeGroups;
};

/////////////////////////////////////////////////////////////////////////////////////////////

class JMeshConstraintsDialog : public QDialog, public Ui::MeshConstraintsDialog
{
    Q_OBJECT

public:
    JMeshConstraintsDialog( QWidget *parent = 0);
    ~JMeshConstraintsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

protected:
    virtual void  mousePressEvent( QMouseEvent *e);
    virtual void  mouseReleaseEvent( QMouseEvent *e);
    virtual void  mouseMoveEvent( QMouseEvent *e);
    virtual void  showEvent( QShowEvent *e);

signals:
    void setNewConstraints();

private slots:
    void init();
    void addBoundaryConstraints();
    void invertConstraints();
    void setSelectionMode();
    void setSelectRegion();
    void openPolygonDialog();
    void getPolygon();
    void saveConstraints();
    void readConstraints();
    void lastDelete();
    void clearAll();
    void closeDialog();

private:
    JaalViewer  *viewManager;

    static int  groupID;
    JMeshPtr  mesh;
    JMeshViewerPtr meshViewer;
    JMeshDeformViewer *deformViewer;
    JMeshEntityPickerPtr picker;


    boost::scoped_ptr<JPolygonDialog> polygonDialog;
    boost::scoped_ptr<JMeshConstraintsViewer>   constraintsViewer;

    vector<boost::shared_ptr<JConstraintSet> > constraintSet;
    vector<Point3D> lassoPoints;
    vector<JNodeSequence>  recordConstraints;

    int   numConstraintNodes;
    bool left_button_pressed;
    JNodePtr moveVertex;
    QString lastSelectedDirectory;

    bool isConstraint( const JNodePtr &v);
    int  newConstraint( const JNodePtr &v);
    void setGlyph();
    void add2Group(const JNodePtr &vtx);
    void eraseRegionPoints(const JNodePtr &vtx);

    void parameterizeLasso();
    void addWithinLasso();
    bool isWithinLasso( const Point3D &xyz);

    void addOnLasso();
    bool isOnLasso( const JEdgePtr edge, const Point3D &p0, const Point3D &p1);
    int  getSign( const Point3D &p0, const Point3D &p1, const Point3D &p2);
    void assignColors();
    size_t  countConstraints();

    void makeConnections();
};

