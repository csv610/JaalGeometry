#pragma once

#include <QDialog>
#include "Ui_QuadVerdictDialog.hpp"
#include "MeshViewer.hpp"

class JQuadVerdictDialog : public QDialog, public Ui::QuadVerdictDialog {
    Q_OBJECT

public:
    JQuadVerdictDialog( QWidget *parent = 0);
    ~JQuadVerdictDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void getData();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;

    void init();
    void makeConnections();
    void setData( const vector<double> &quality, double minAccept, QLineEdit *minLineEdit,
                  double maxAccept,  QLineEdit *maxLineEdit);
    double getAcceptable( const vector<double> &quality, double minval, double maxval);

};
////////////////////////////////////////////////////////////////////////////////
