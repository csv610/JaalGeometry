#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_SceneFloorDialog.hpp"
#include "MeshViewer.hpp"
#include "SceneFloor.hpp"

class JSceneFloorDialog : public QDialog, public Ui::SceneFloorDialog {
    Q_OBJECT

public:
    JSceneFloorDialog( QWidget *parent = 0);
    ~JSceneFloorDialog();

    void setViewManager( JaalViewer *v);

private slots:
    void setPattern();
    void setColor();
    void checkDisplay();
    void setDirection();
    void setDistance();
    void setLength();
    void setLines();
    void keyPressEvent( QKeyEvent *e);

private:
    JaalViewer *viewManager;
    boost::shared_ptr<JSceneFloor> sceneFloor;
    int pattern;

    void init();
    void makeConnections();

};

