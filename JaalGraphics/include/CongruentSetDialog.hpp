#pragma once

#include <QDialog>

#include "Ui_CongruentSetDialog.hpp"
#include "MeshViewer.hpp"
#include "CongruentSet.hpp"

class JCongruentSetDialog : public QDialog, public Ui::CongruentSetDialog {
    Q_OBJECT

public:
    JCongruentSetDialog( QWidget *parent = 0);
    ~JCongruentSetDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void generateSet();
    void displayGroup();
    void getCongruent();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    CongruentSet cset;

    void init();
    void makeConnections();
};
