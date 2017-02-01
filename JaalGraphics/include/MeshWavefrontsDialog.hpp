#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshWavefrontsDialog.hpp"
#include "MeshViewer.hpp"


class JMeshWavefrontsDialog : public QDialog, public Ui::MeshWavefrontsDialog {
    Q_OBJECT

public:
    JMeshWavefrontsDialog( QWidget *parent = 0);
    ~JMeshWavefrontsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:

    void nextWave();
    void waveAnimation();
    void closeDialog();
    void setNewWave();
    void keyPressEvent( QKeyEvent *e);

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    int currWaveID, maxWaveID, waveSteps;
    map<int, JCellSequence>  cellwaves;

    void setWave(int wid);

    void init();
    void makeConnections();
};

