#pragma once

#include <QDialog>

#include "Ui_MorseAnalysisDialog.hpp"
#include "MeshViewer.hpp"
#include "MorseAnalysis.hpp"
#include "MarchingTriangles.hpp"

class JMorseAnalysisDialog : public QDialog, public Ui::MorseAnalysisDialog {
    Q_OBJECT

public:
    JMorseAnalysisDialog( QWidget *parent = 0);
    ~JMorseAnalysisDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void genHeight();
    void genContours();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    vector<JEdgeSequence> contours;

    void init();
    void makeConnections();
};

