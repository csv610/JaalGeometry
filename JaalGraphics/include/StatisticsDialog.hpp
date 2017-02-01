#pragma once

#include <QDialog>
#include <QColorDialog>
#include <QThread>
#include <vector>
#include <numeric>

#include "MeshViewer.hpp"
#include "Ui_StatisticsDialog.hpp"
#include "basic_math.hpp"
#include "Mesh.hpp"

using namespace std;

class JStatisticsDialog : public QDialog, public Ui::StatisticsDialog
{
    Q_OBJECT

    /*
        class MyThread : public QThread {
        public:
            MyThread();

            int     entityDim;
            Color   redColor, greenColor, blueColor;
            JMeshViewer *meshViewer;
            vector<double> data;

            void run();
            void execute();
        private:
        };
        boost::scoped_ptr<MyThread>  thread;
    */

public:
    JStatisticsDialog( QWidget *parent = 0);
    ~JStatisticsDialog() {}

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);
    void setEntityDimension(int entity)  {
        entityDim = entity;
    }
    void setQualityName( const string &s);
    void setValues();

private slots:
    void setSliderValue();
//    void getUpperEntities();
    void editLowerCutoff();
    void editUpperCutoff();
    void resetLowerCutoff();
    void resetUpperCutoff();
    void getCount();

    void keyPressEvent( QKeyEvent *e);

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr    mesh;

    int     entityDim;
    string  qualityName;

    double  minVal, maxVal;
    double  lowerCutoff, upperCutoff;
    bool    displayLower, displayUpper, displayMiddle;
    JColor  greenColor, redColor, blueColor;
    vector<double> data;

    void init();
    void makeConnections();
    void display(const JEdgePtr &e);
    void display(const JFacePtr &f);
    void display(const JCellPtr &c);
    void display_edges_quality();
    void display_faces_quality();
    void display_cells_quality();
};

