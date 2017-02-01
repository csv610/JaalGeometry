#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_LightsDialog.hpp"
#include "JaalViewer.hpp"
#include "Lights.hpp"

class JLightsDialog : public QDialog, public Ui::LightsDialog {
    Q_OBJECT

public:
    JLightsDialog( QWidget *parent = 0);
    ~JLightsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void init();

    void ambientModel();
    void lightModel();

    void getLightNum();
    void setLightType();
    void switchLight();
    void ambientColor();
    void diffuseColor();
    void specularColor();

    void setPosition();

    void spotDirection();
    void spotExponent();
    void spotCutoff();

    void lightBulbs();

    void  keyPressEvent( QKeyEvent *e);

private:
    JaalViewer *viewManager;
    boost::shared_ptr<JLights> lights;
    vector<GLenum> lightID;
    void makeConnections();
};
