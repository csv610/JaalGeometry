#pragma once

#include <QDialog>

#include "Ui_TetMesherDialog.hpp"
#include "MeshViewer.hpp"
#include "AllTetMeshGenerator.hpp"
#include "TetGenOptionsDialog.hpp"
#include "MeshOptDialog.hpp"
#include "ImplicitMeshCutterDialog.hpp"

class JMeshOptDialog;

class JTetMesherDialog : public QDialog, public Ui::TetMesherDialog {
    Q_OBJECT

    struct MyThread : public QThread {
        void run();
        string cmd;
        JMeshPtr mesh;
        JTetMesherDialog *dialog;
    };

    friend class MyThread;

public:
    JTetMesherDialog( QWidget *parent = 0);
    ~JTetMesherDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

protected:
    virtual void  keyPressEvent( QKeyEvent *e);

private slots:
    void genNewMesh();
    void fromHexMesh();
    void openTetGenOptionsDialog();
    void getBoundedDistortion();
    void getMesquiteOpt();
    void getStellarOpt();
    void openTetViewerDialog();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;
    JMeshPtr newTetMesh;

    boost::scoped_ptr<MyThread>  thread;
    boost::scoped_ptr<JTetGenOptionsDialog>  tetgenOptionsDialog;
    boost::scoped_ptr<JMeshOptDialog>  meshOptDialog;
    boost::scoped_ptr<JImplicitMeshCutterDialog>  meshCutterDialog;

    void init();
    void makeConnections();
};
