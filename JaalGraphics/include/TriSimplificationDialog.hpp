#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_TriSimplificationDialog.hpp"
#include "MeshViewer.hpp"

//#include "qslim.h"

class JTriSimplificationDialog : public QDialog, public Ui::TriSimplificationDialog {
    Q_OBJECT

public:
    JTriSimplificationDialog( QWidget *parent = 0);
    ~JTriSimplificationDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:

    void accept();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    vector<JMeshPtr>  levelmesh;
    /*
         boost::scoped_ptr<MxSMFReader>  smfReader;
         boost::scoped_ptr<MxStdModel>   inmodel;
    */

    int    placement_policy;
    int    weighting_policy;
    int    boundary_weight;
    int    face_target;
    bool   use_fslim;
    bool   record_history;
    bool   join_only;
    double meshing_penalty;
    double compactness_ratio;

    void init();
    void makeConnections();
    void startup_and_init();
    void write_smf();
    void read_smf();
};
