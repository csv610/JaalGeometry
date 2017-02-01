#pragma once

#include <QDialog>
#include "Ui_ScreenShotDialog.hpp"

#include "JaalViewer.hpp"
#include "gl2ps.hpp"

class JScreenShotDialog : public QDialog, public Ui::ScreenShotDialog {
    Q_OBJECT

public:
    JScreenShotDialog( QWidget *parent = 0);
    ~JScreenShotDialog();

    void setViewManager( JaalViewer *v);

private slots:
    void getShot();

private:
    JaalViewer *viewManager;
    int  prevWidth, prevHeight;
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
