#pragma once

#include "Ui_GlobalSettingsDialog.hpp"

#include <QDialog>
#include <QColorDialog>

#include "CameraDialog.hpp"
#include "FontsDialog.hpp"
#include "LightsDialog.hpp"
#include "SceneFloorDialog.hpp"
#include "ScreenShotDialog.hpp"

class JGlobalSettingsDialog : public QDialog, public Ui::GlobalSettingsDialog {
    Q_OBJECT

public:
    JGlobalSettingsDialog( QWidget *parent = 0);
    ~JGlobalSettingsDialog();

    void setViewManager(JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void  openFontsDialog();
    void  openLightsDialog();
    void  openFloorDialog();
    void  openCameraDialog();
    void  openScreenShotDialog();

    void  closeDialog();
    void  keyPressEvent( QKeyEvent *e);

    void setBackgroundColor();

    void checkBoundingBox();
    void checkAxis();

private:
    JaalViewer *viewManager;

    boost::scoped_ptr<JFontsDialog> fontsDialog;
    boost::scoped_ptr<JCameraDialog> cameraDialog;
    boost::scoped_ptr<JSceneFloorDialog> floorDialog;
    boost::scoped_ptr<JLightsDialog> lightsDialog;
    boost::scoped_ptr<JScreenShotDialog> screenShotDialog;

    void makeConnections();
    void init();
};
