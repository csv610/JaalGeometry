#pragma once

#include <QDialog>

#include "JaalViewer.hpp"
#include "Ui_CameraDialog.hpp"

class JCameraDialog : public QDialog, public Ui::CameraDialog {
    Q_OBJECT

public:
    JCameraDialog( QWidget *parent = 0);
    ~JCameraDialog();

    void setViewManager( JaalViewer *v);

protected:
    virtual void keyPressEvent( QKeyEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);
    virtual void wheelEvent(QWheelEvent *e);


private slots:
    void setViewSide();
    void setViewPoint();
    void setCamera();
    void setProjection();
    void entireScene();
    void setSceneRadius();
    void setRotationConstraint();
    void setTranslationConstraint();
    void setPosition();
    void setWindowSize();
    void closeDialog();

private:
    JaalViewer *viewManager;

    void setValues();
    void makeConnections();

};

