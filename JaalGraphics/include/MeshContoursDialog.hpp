#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshContoursDialog.hpp"
#include "MeshViewer.hpp"
//#include "PolygonSimplifyDialog.hpp"

class JMeshContoursDialog : public QDialog, public Ui::MeshContoursDialog {
    Q_OBJECT

public:
    JMeshContoursDialog( QWidget *parent = 0);
    ~JMeshContoursDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void displayContour();
    void getOriginal();
    void getSmooth();
    void getReparam();
    void getSimplify();
    void closeDialog();
    void keyPressEvent( QKeyEvent *e);

protected:
    void  mousePressEvent( QMouseEvent *e);
    void  mouseReleaseEvent( QMouseEvent *e);

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshEntityPickerPtr picker;

    JMeshPtr mesh;
    vector<JEdgeSequence>  boundedges;

    vector<double> orgCoords;
    vector<size_t> l2g;

//    boost::scoped_ptr<JPolygonSimplifyDialog> simplifyDialog;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
