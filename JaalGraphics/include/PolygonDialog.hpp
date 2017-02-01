#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_PolygonDialog.hpp"
#include "MeshViewer.hpp"

class JPolygonViewer : public JViewComponent
{
public:
    JPolygonViewer();

    void setNumSides(int n)
    {
        numSides = n;
    }

    void setColor(const JColor &bg)
    {
        color[0] = bg[0];
        color[1] = bg[1];
        color[2] = bg[2];
        color[3] = bg[3];
    }

    void setStartPos(const Point3D &pstart)
    {
        started   = 1;
        startPos  = pstart;
        endPos    = pstart;
    }

    void setEndPos(const Point3D &pend)
    {
        endPos   = pend;
    }

    void setCanvas(bool v ) {
        canvas = v;
    }

    const vector<Point3D> getPoints()
    {
        return polyPoints;
    }

    void draw();

private:
    bool canvas;
    int  numSides;
    bool started, finished;
    double lineWidth;

    JColor  color;
    Point3D startPos;
    Point3D endPos;
    vector<Point3D> polyPoints;

    void  genPoints();
    void  drawRectangle();
    void  drawPolygon();
};

////////////////////////////////////////////////////////////////////////////

class JPolygonDialog : public QDialog, public Ui::PolygonDialog
{
    Q_OBJECT
public:
    JPolygonDialog( QWidget *parent = 0);
    ~JPolygonDialog();

    void setViewManager( JaalViewer *v);
    const vector<Point3D> &getPoints()
    {
        return polyPoints;
    }

protected:
    virtual void  mousePressEvent( QMouseEvent *e);
    virtual void  mouseReleaseEvent( QMouseEvent *e);
    virtual void  mouseMoveEvent( QMouseEvent *e);

signals:
    void setPolygon();

private slots:
    void genShape();
    void setNumSides();
    void closeDialog();
    void keyPressEvent( QKeyEvent *e);
    void setGenMode();

private:
    JaalViewer  *viewManager;
    bool  left_button_pressed;
    vector<Point3D> polyPoints;

    boost::shared_ptr<JPolygonViewer> polyViewer;

    void makeConnections();
};
