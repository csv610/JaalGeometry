#pragma once

#include "Ui_MaterialDialog.hpp"

#include <QDialog>
#include <QColorDialog>
#include "MeshViewer.hpp"

class JMaterialDialog : public QDialog, public Ui::MaterialDialog {
    Q_OBJECT

public:
    JMaterialDialog( QWidget *parent = 0);
    ~JMaterialDialog();

    void setMeshViewer( JMeshViewer *v) {
        meshViewer = v;
    }

    void update(JMaterial *m)  {
        material = m;
    }

private slots:
    void setAmbient();
    void setDiffuse();
    void setEmission();
    void setSpecular();
    void setStdMaterial();

private:
    Array3F rgb;
    Array4F rgba;
    JMaterial *material;

    JMeshViewer *meshViewer;
    void makeConnections();
    void setColor();
    JMaterial  newmaterial;
};
