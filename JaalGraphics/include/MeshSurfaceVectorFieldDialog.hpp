#pragma once

#include <QDialog>

#include "Ui_MeshSurfaceVectorFieldDialog.hpp"
#include "MeshViewer.hpp"

#include "EdgeAttributesDialog.hpp"
#include "NodeAttributesDialog.hpp"

#include "NRoSyField.hpp"
#include "PolyVectorsField.hpp"
#include "ConjugateField.hpp"
#include "IntegrableField.hpp"
#include "MeshCurvature.hpp"

class JVectorFieldViewer : public JViewComponent
{
public:

    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }
    void setField( const JMeshPtr &f) {
        vecField = f;
    }
    void setNumSamplesPerFace( int n) {
        numSamples = n;
    }
    void genField();

    void draw();
private:
    JMeshPtr mesh;
    JMeshPtr vecField;
    int numSamples = 3;

    struct FieldVec {
        float  head[3];
        float  tail[3];
    };
    vector<FieldVec>  fieldVec;
    void sampleField( const JFacePtr &t);
};


class JMeshSurfaceVectorFieldDialog : public QDialog, public Ui::MeshSurfaceVectorFieldDialog {
    Q_OBJECT

    const static int   NROSY_FIELD       = 0;
    const static int   POLYVECTORS_FIELD = 1;
    const static int   CONJUGATE_FIELD   = 2;

public:
    JMeshSurfaceVectorFieldDialog( QWidget *parent = 0);
    ~JMeshSurfaceVectorFieldDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

protected:
    void  mouseReleaseEvent( QMouseEvent *e);
    void  keyPressEvent( QKeyEvent *e);

private slots:
    void getVecField();
    void clearField();
    void setVecLength();
    void openEdgeAttribsDialog();
    void readFixedFacesID();
    void readFixedFacesVectors();
    void setConstraints();
    void closeDialog();
    void showAllVectors();
    void showVecComponent();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
//   boost::shared_ptr<JVectorFieldViewer>    vecFieldViewer;

    JMeshEntityPickerPtr picker;
    JFaceSequence pickedFaces;

//  JNodeSequence singularNodes;
    JFaceSequence constrainedFaces;
    JFaceSequence singularFaces;

    int currFieldGenerator;
    string constraintFiles[2];

//    QString lastSelectedDirectory;

    double length;

    JFacePtr currPickedFace;
    JMeshPtr mesh;
    vector<JMeshPtr> vecFields;
// constrainedMesh;

    boost::scoped_ptr<JMeshSurfaceVectorField> surfVecFieldPtr;
/*
    boost::scoped_ptr<JNRoSyField>       nRoSyFieldPtr;
    boost::scoped_ptr<JPolyVectorsField> polyVecFieldPtr;
    boost::scoped_ptr<JConjugateField>   conjugateFieldPtr;
    boost::scoped_ptr<JIntegrableField>  integrableFieldPtr;
*/
    boost::scoped_ptr<JEdgeAttributesDialog> edgeAttribsDialog;
    boost::scoped_ptr<JNodeAttributesDialog> nodeAttribsDialog;

    void init();
    void makeConnections();
    void getNRoSyField();
    void getPolyVecField();
    void getConjugateField();
    void getIntegrableField();
    void displayField(JMeshSurfaceVectorField *p);
};
////////////////////////////////////////////////////////////////////////////////
