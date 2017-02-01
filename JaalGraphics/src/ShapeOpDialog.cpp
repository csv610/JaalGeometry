#include "ShapeOpDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JShapeOpDialog :: JShapeOpDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    shapeOpt.reset( new JShapeOptimizer);
}

///////////////////////////////////////////////////////////////////////////////

JShapeOpDialog :: ~JShapeOpDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    mesh->buildRelations(0,2);
    shapeOpt->setMesh(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyAreaConstraints()
{
    cout << "Area Constraint " << endl;
    if( areaCheckBox->isChecked() )
        shapeOpt->addAreaConstraints();
    else
        shapeOpt->removeAreaConstraints();
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyCocircularConstraints()
{
    shapeOpt->addCircleConstraints();
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyCoplanarConstraints()
{
    shapeOpt->addPlaneConstraints();
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyEdgeStrainConstraints()
{
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyFixedLengthConstraints()
{
//   shapeOpt->addLengthConstraint();
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyLaplaceConstraints()
{
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyParallelogramConstraints()
{
    shapeOpt->addParallelogramConstraints();
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyTriangleStrainConstraints()
{
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyRectangleConstraints()
{
    shapeOpt->addRectangleConstraints();
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyClosenessConstraints()
{
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyFixedBoundaryConstraints()
{
    shapeOpt->addBoundaryConstraints();
}
///////////////////////////////////////////////////////////////////////////////
void JShapeOpDialog :: applyCrossFieldConstraints()
{
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applySimilarityConstraints()
{
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: applyRigidConstraints()
{
}
///////////////////////////////////////////////////////////////////////////////
void JShapeOpDialog :: displayConstraints()
{
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: solve()
{
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: closeDialog()
{
}
///////////////////////////////////////////////////////////////////////////////

void JShapeOpDialog :: makeConnections()
{
    CheckBox( areaCheckBox, [=] {applyAreaConstraints();});
    CheckBox( cocircularCheckBox, [=] {applyCocircularConstraints();});
    CheckBox( coplanarCheckBox, [=] {applyCoplanarConstraints();});
    CheckBox( crossFieldCheckBox, [=] {applyCrossFieldConstraints(); });
    CheckBox( edgeStrainCheckBox, [=] {applyEdgeStrainConstraints() ;});
    CheckBox( laplaceFairCheckBox, [=] {applyLaplaceConstraints() ;});
    CheckBox( parallelogramCheckBox, [=] {applyParallelogramConstraints() ;});
    CheckBox( triangleStrainCheckBox, [=] {applyTriangleStrainConstraints();});
    CheckBox( rectangleCheckBox, [=] {applyRectangleConstraints();});
    CheckBox( vertexClosenessCheckBox, [=] {applyClosenessConstraints(); });
    CheckBox( similarityCheckBox, [=] {applySimilarityConstraints();});
    CheckBox( rigidCheckBox, [=] {applyRigidConstraints();});

    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
