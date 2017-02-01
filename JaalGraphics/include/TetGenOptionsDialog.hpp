#pragma once

#include "Ui_TetGenOptionsDialog.hpp"
#include "MeshViewer.hpp"
#include <sstream>

class JTetGenOptionsDialog : public QDialog, public Ui::TetgenOptionsDialog {
    Q_OBJECT

public:
    JTetGenOptionsDialog( QWidget *parent = 0);
    ~JTetGenOptionsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    string getOptions() const {
        return options;
    }

private slots:
    void setOptions();

private:
    JaalViewer *viewManager;
    string   options;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
