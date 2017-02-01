#pragma once

#include <boost/smart_ptr.hpp>

#include <QMainWindow>
#include <QMessageBox>
#include <QFileDialog>
#include <QAction>

#include "Ui_JaalMainWindow.hpp"

#include "JaalViewer.hpp"
#include "ImageViewer.hpp"

#include "LoadNewDataDialog.hpp"
#include "GlobalSettingsDialog.hpp"
#include "MeshToolsDialog.hpp"

class JaalMainWindow  : public QMainWindow, public Ui::JaalMainWindow {
    Q_OBJECT

public:
    explicit JaalMainWindow( QWidget *parent = NULL );
    ~JaalMainWindow();

private slots:
    void resizeEvent( QResizeEvent *e);
    void openNewDataDialog();
    void openMeshToolsDialog();
    void openGlobalSettingsDialog();
    void openObjectsListDialog();
    void Quit();


private:
    boost::shared_ptr<JMeshViewer>      meshViewer;
    boost::shared_ptr<JImageViewer>     imageViewer;

    boost::scoped_ptr<JLoadNewDataDialog>   loadNewDataDialog;
    boost::scoped_ptr<JGlobalSettingsDialog> globalSettingsDialog;
    boost::scoped_ptr<JMeshToolsDialog> meshToolsDialog;
    boost::scoped_ptr<JObjectsListDialog>      objectsListDialog;

    void makeConnections();
    void clearup();
};

