#pragma once

#include <QDialog>
#include <sstream>
#include "MeshViewer.hpp"
#include "ImageViewer.hpp"
#include "Ui_LoadNewDataDialog.hpp"
#include <thread>

class JLoadNewDataDialog : public QDialog, public Ui::LoadNewDataDialog {
    Q_OBJECT

    struct ThreadWork1
    {
        ThreadWork1() {
            finished = -1;
        }
        int     finished;
        std::string  name;
        JMeshViewerPtr meshViewer;

        void  run()
        {
            finished = 0;
            JMeshPtr msh = JMeshIO::readFile(name);
            if( msh) {
                msh->setFileName(name);
                meshViewer->addObject(msh);
                JMeshPtr texmesh;
                msh->getAttribute("TextureMesh", texmesh);
                if( texmesh ) meshViewer->addObject(texmesh);
            }
            finished = 1;
        }

        void  operator() ()
        {
            run();
        }
    };

    JaalViewer   *viewManager;
public:
    JLoadNewDataDialog( QWidget *parent = 0);
    ~JLoadNewDataDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void loadNewData();
    void closeDialog();

private:
    JMeshViewerPtr  meshViewer;
    JImageViewerPtr imageViewer;
    boost::scoped_ptr<std::thread> readThread;

    QString lastSelectedDirectory;

    void init();
    void makeConnections();

    void loadNewMesh();
    void loadNewImage();
};
