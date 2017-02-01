#include "MeshStackDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshStackDialog :: JMeshStackDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    oldmesh = nullptr;
    hexmesh = nullptr;
    trajCurve = nullptr;
    baseQuadmesh = nullptr;
    rotAngle = 0.0;

    rotateLineEdit->setText( QString::number(rotAngle) );
}

///////////////////////////////////////////////////////////////////////////////

JMeshStackDialog :: ~JMeshStackDialog()
{ }

///////////////////////////////////////////////////////////////////////////////
void JMeshStackDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);


    /*
        if( meshViewer == nullptr ) return;
        oldmesh = meshViewer->getMesh();
        baseQuadmesh = Mesh::newObject();
    */

}
///////////////////////////////////////////////////////////////////////////////

void JMeshStackDialog :: refineQuads()
{
    if( refineDialog == nullptr )
        refineDialog.reset(new JMeshRefine2DDialog(this));
    refineDialog->setViewManager(viewManager);

    baseQuadmesh->deleteAll();

    JQuadrilateralPtr quad = JQuadrilateral::getCanonical();
    JNodeSequence qnodes = quad->getNodes();
    baseQuadmesh->addObjects( qnodes);
    baseQuadmesh->addObject(quad);
    meshViewer->addObject(baseQuadmesh);
    refineDialog->show();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshStackDialog :: gen3DMesh()
{
    /*
        if( hexmesh ) hexmesh->deleteAll();

        if( hexmesh == nullptr ) hexmesh = JMesh::newObject();

        QString qstr = rotateLineEdit->text();
        double rotAngle = qstr.toDouble();

        JNodeSequence curvenodes, meshnodes, qnodes;
        JFaceSequence meshfaces;
        if( curveDialog )
        {
            JCurve *newCurve = curveDialog->getCurve();
            if( newCurve )
            {
                newCurve->getNodes( curvenodes );
                int numPlanes = curvenodes.size();
                cout << "#Planes " << numPlanes << endl;
                meshPlanes.resize( numPlanes );
                Vec3D srcVec, dstVec;
                srcVec[0] = 0.0;
                srcVec[1] = 0.0;
                srcVec[2] = 1.0;
                for( int i = 0; i < numPlanes; i++)
                {
                    meshPlanes[i] = baseQuadmesh->deepCopy();
                    newCurve->getTangentAt( curvenodes[i], dstVec);
                    meshViewer->alignAlong( meshPlanes[i], srcVec, dstVec);
                    AffineTransform af( meshPlanes[i] );
                    Point3D  p3d = curvenodes[i]->getXYZCoords();
                    af.translate( p3d[0], p3d[1], p3d[2] );
                    meshPlanes[i]->getNodes(qnodes);
    //                    af.rotate( qnodes, rotAngle, 2);
                    meshPlanes[i]->getNodes( meshnodes );
                    hexmesh->addObjects( meshnodes );

                    meshPlanes[i]->getFaces( meshfaces );
                    hexmesh->addObjects( meshfaces );
                }
                meshViewer->addObject(hexmesh);
            }
        }

        hexmesh = AllHexMeshGenerator::stackQuadMesh(meshPlanes);
        meshViewer->addObject(hexmesh);
        meshViewer->refreshDisplay();
    */
}
///////////////////////////////////////////////////////////////////////////////

void JMeshStackDialog :: openCurveDialog()
{
    if( hexmesh )  hexmesh->deleteAll();

    if( curveDialog == nullptr)
        curveDialog.reset(new JCurveGenDialog(this));

    curveDialog->setViewManager( viewManager );
    curveDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshStackDialog :: makeConnections()
{
    connect( trajectoryCurvePushButton,  SIGNAL( clicked() ), this, SLOT( openCurveDialog() ));
    connect( refinePushButton, SIGNAL( clicked() ), this, SLOT( refineQuads() ));
    connect( applyPushButton, SIGNAL( clicked() ), this, SLOT(  gen3DMesh() ));

    connect( closePushButton,  SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
