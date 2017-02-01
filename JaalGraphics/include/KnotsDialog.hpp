#pragma once

#include <QDialog>
#include <QFileDialog>

#include "Ui_KnotsDialog.hpp"
#include "MeshViewer.hpp"
#include "Curve.hpp"

class JKnotsDialog : public QDialog, public Ui::KnotsDialog
{
    Q_OBJECT

public:
    JKnotsDialog( QWidget *parent = 0);
    ~JKnotsDialog();

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }

    JCurve*  getKnot() const
    {
        return newCurve;
    }

private slots:
    void readFromFile();
    void closeDialog();


private:
    JMeshPtr oldmesh, newmesh;
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    JCurve *newCurve;
    QString lastSelectedDirectory;

    void init();
    void makeConnections();
};
