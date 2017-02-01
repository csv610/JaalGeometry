#pragma once

#include <QDialog>

#include "Ui_ImageContoursDialog.hpp"
#include "ImageViewer.hpp"
#include "MeshViewer.hpp"
#include "MeshContoursDialog.hpp"

////////////////////////////////////////////////////////////////////////////////

struct JContour {
    JContour() {
        closed = 0;
    }
    bool  closed;
    vector<Point3D> points;
};

////////////////////////////////////////////////////////////////////////////////

class JImageContoursViewer : public JViewComponent
{
public:
    JImageContoursViewer();

    void init();

    void setDrawPoints( bool v) {
        drawPoints = v;
    }

    void setDrawID( bool v) {
        drawIDs = v;
    }

    void append(JContour *c) {
        contours.push_back(c);
    }

    void drawPixelCircle( bool b){ pixelCircle = b; }

    void setNodeRadius( double r) {
        nodeRadius = r;
    }
    void setEdgeWidth( double r)  {
        edgeWidth  = r;
    }
    void setNodeColor( const JColor &c)   {
        nodeColor   = c;
    }
    void setEdgeColor( const JColor &c)   {
        edgeColor   = c;
    }

    void clearAll()  {
        contours.clear();
    }

    void draw();

private:
    bool  drawIDs     = 0;
    bool  drawPoints  = 0;
    bool  pixelCircle = 1;

    double nodeRadius, edgeWidth;
    JColor  nodeColor, edgeColor;
    vector<JContour*> contours;
    void draw( JContour *c);
};

////////////////////////////////////////////////////////////////////////////////

class JImageContoursDialog : public QDialog, public Ui::ImageContoursDialog
{
    Q_OBJECT

public:

    JImageContoursDialog( QWidget *parent = 0);
    ~JImageContoursDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

protected:
    virtual void  mousePressEvent( QMouseEvent *e);
    virtual void  mouseReleaseEvent( QMouseEvent *e);
    virtual void  mouseMoveEvent( QMouseEvent *e);
    virtual void  keyPressEvent( QKeyEvent *e);

private slots:

    void readFile();
    void newContour();
    void setNodeRadius();
    void setEdgeWidth();
    void setNodeColor();
    void setEdgeColor();
    void deleteAll();
    void displayIDs();
    void deleteSegment();
    void closeContour();
    void saveAs();
    void closeDialog();
    void genEdgeMesh();

private:
    JaalViewer *viewManager;

    JImageViewerPtr imageViewer;
    boost::shared_ptr<JImageContoursViewer>   contourViewer;
    bool   leftButton;
    bool   freehand;

    JContour   *currContour;
    vector<JContour*> contours;
    boost::scoped_ptr<JMeshContoursDialog> meshContoursDialog;

    void init();
    void saveFile(string &s);
    void makeConnections();
    void mergeClosePoints();
    void saveOFF( string &s);
    void saveXML( string &s);
};

