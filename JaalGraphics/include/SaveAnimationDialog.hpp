#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_SaveAnimationDialog.hpp"
#include "MeshViewer.hpp"


class JSaveAnimationDialog : public QDialog, public Ui::SaveAnimationDialog {
    Q_OBJECT

public:
    JSaveAnimationDialog( QWidget *parent = 0);
    ~JSaveAnimationDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
    }

    /*
    private slots:
         void init();
         void setColor();
         void removeAll();
         void accept();
         void reject();
    */

private:
    JaalViewer *viewManager;
    void makeConnections();
};
