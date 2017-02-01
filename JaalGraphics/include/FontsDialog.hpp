#pragma once

#include <QDialog>

#include "Ui_FontsDialog.hpp"
#include "MeshViewer.hpp"

class JFontsDialog : public QDialog, public Ui::FontsDialog {
    Q_OBJECT

public:
    JFontsDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v) {
        viewManager = v;
    }

private slots:
    void xRotate();
    void yRotate();
    void zRotate();
    void setFontSize();
    void setColor();

private:
    JaalViewer *viewManager;
    void makeConnections();
};
