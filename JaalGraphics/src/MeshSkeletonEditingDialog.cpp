#include "MeshSkeletonEditingDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshSkeletonEditingDialog :: JMeshSkeletonEditingDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshSkeletonEditingDialog :: ~JMeshSkeletonEditingDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonEditingDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c1 = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c1);

    JViewComponentPtr c2 = viewManager->getComponent("ShapeViewer");
    if( c2 == nullptr)
        c2 = JShapeViewer::registerComponent(viewManager);

    shapesViewer = dynamic_pointer_cast<JShapeViewer>(c2);
}

///////////////////////////////////////////////////////////////////////////////


void JMeshSkeletonEditingDialog :: setSkeleton( const JMeshSkeletonPtr &s) {
    meshSkel = s;
    if( meshSkel == nullptr)return;
    skelGraph = meshSkel->getSkeleton();
    surfMesh  = meshSkel->getMesh();
    int nbr   = meshSkel->getNumBranches();
    selectBranchSpinBox->setMaximum(nbr);
    skelContourPtr.reset( new JMeshSkeletonContours);
    skelContourPtr->setSkeleton(s);

}

///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonEditingDialog :: displayBranch()
{
    int branchID = selectBranchSpinBox->value();

    JMeshSkeletonBranchPtr branch = meshSkel->getBranch(branchID);
    if( branch == nullptr) return;

    JEdgeSequence branchEdges = branch->skelEdges;
    int numEdges = branchEdges.size();
    numBranchEdgesLineEdit->setText( QString::number(numEdges));

    double len = branch->getLength();
    branchLengthLineEdit->setText( QString::number(len));

    if( changeCenterCheckBox->isChecked() ) {
        Point3D startPoint = branchEdges[0]->getNodeAt(0)->getXYZCoords();
        JMeshAffineTransform  affine;
        affine.setMesh(surfMesh);
        affine.translate( -startPoint[0], -startPoint[1], -startPoint[2] );

        affine.setMesh(skelGraph);
        affine.translate( -startPoint[0], -startPoint[1], -startPoint[2] );

        Point3D p0 = branchEdges.front()->getNodeAt(0)->getXYZCoords();
        Point3D p1 = branchEdges.back()->getNodeAt(1)->getXYZCoords();

        Vec3D src;
        src[0] = p1[0] - p0[0];
        src[1] = p1[1] - p0[1];
        src[2] = p1[2] - p0[2];

        Vec3D dst;
        dst[0] = 0.0;
        dst[1] = 1.0;
        dst[2] = 0.0;

        JNodeSequence  anodes = surfMesh->getNodes();
        JNodeSequence  bnodes = skelGraph->getNodes();
        for( const JNodePtr &vtx : bnodes)
            anodes.push_back(vtx);
        affine.alignAlong(anodes, src,dst);
    }

    JMeshRenderPtr mrender;
    surfMesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 1;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 1;

    JNodeRenderPtr vAttrib;
    JColor color = JEntityColor::getColor("LightGray");
    size_t numNodes = surfMesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &node = surfMesh->getNodeAt(i);
        node->getAttribute("Render", vAttrib);
        vAttrib->display = 1;
        vAttrib->color   = color;
    }

    JFaceRenderPtr fAttrib;
    int partID = -1;
    size_t numFaces = surfMesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = surfMesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display = 0;
        face->getAttribute("Partition", partID);
        if( partID == branch->id ) fAttrib->display = 1;
    }

    if( displaySurfMeshCheckBox->isChecked() ) {
        for( size_t i = 0; i < numFaces; i++) {
            const JFacePtr &face = surfMesh->getFaceAt(i);
            face->getAttribute("Render", fAttrib);
            fAttrib->display = 0;
            face->getAttribute("Partition", partID);
            fAttrib->display = 1;
        }
    }

    if( displayAllEndSlicesCheckBox->isChecked() ) {
        if( otherBranchesMesh == nullptr) {
            int numBranches = meshSkel->getNumBranches();
            otherBranchesMesh = JMesh::newObject();
            for( int i = 0; i < numBranches; i++) {
                 if( i != branchID) {
                   JSlicePtr slice1  = skelContourPtr->getFirstSlice(i);
                   otherBranchesMesh->addObject( slice1->mesh);
                   JSlicePtr slice2 =  skelContourPtr->getLastSlice(i);
                   otherBranchesMesh->addObject( slice2->mesh);
                 }
            }
        }
        meshViewer->addObject(otherBranchesMesh); // For default display.
        displayContours(otherBranchesMesh);       // Now you can change the properties.
    } else {
       meshViewer->removeObject(otherBranchesMesh);
    }

    meshViewer->updateBuffers(surfMesh);
    meshViewer->updateBuffers(skelGraph);

}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonEditingDialog :: selectBranch()
{
    if( meshSkel == nullptr) return;
    if( selectedBranchMesh ) {
       meshViewer->removeObject( selectedBranchMesh);
       selectedBranchMesh.reset();
    }

    if( otherBranchesMesh )  {
        meshViewer->removeObject( otherBranchesMesh);
        otherBranchesMesh.reset();
    }

    displayBranch();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonEditingDialog :: getCap()
{
    cout << "Get Cap " << endl;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonEditingDialog :: getPoles()
{
    cout << "Get Poles " << endl;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonEditingDialog :: removeBranch()
{
    cout << "Remove branch " << endl;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonEditingDialog :: reparameterize()
{
    if( meshSkel == nullptr) return;

    int branchID = selectBranchSpinBox->value();
    JMeshSkeletonBranchPtr branch = meshSkel->getBranch(branchID);
    if( branch == nullptr) return;
    branch->reparameterize();
    meshViewer->updateBuffers(skelGraph);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonEditingDialog :: splitBranch()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonEditingDialog :: getContours()
{
    if( skelGraph == nullptr) return;

    JWaitCursor wCursor;
    wCursor.start();

    if( selectedBranchMesh ) meshViewer->removeObject( selectedBranchMesh);
    int branchID = selectBranchSpinBox->value();

    if( allSlicesRadioButton->isChecked() )
        getAllBranchContours(branchID);
    else if( bothEndSlicesRadioButton->isChecked() ) 
        getTwoBranchContours(branchID);
    else
        getOneBranchContour(branchID);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonEditingDialog :: displayContours( const JMeshPtr &cmesh)
{
    cmesh->getTopology()->searchBoundary();
    JEdgeRenderPtr eAttrib;
    JColor blue = JEntityColor::getColor("Blue");

    int numEdges = cmesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = cmesh->getEdgeAt(i);
        edge->getAttribute("Render", eAttrib);
        eAttrib->display   = 0;
        if( edge->getNumRelations(2) == 1) {
            eAttrib->lineWidth = 3;
            eAttrib->display   = 1;
            eAttrib->color     = blue;
        }
    }

    JColor red = JEntityColor::getColor("Red");
    JFaceRenderPtr fAttrib;
    int numFaces = cmesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = cmesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->color  = red;
    }

    meshViewer->updateBuffers(cmesh);
}
//////////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonEditingDialog :: getOneBranchContour(int branchID)
{
    currSlice.reset();
    if( firstSliceRadioButton->isChecked() )
        currSlice  = skelContourPtr->getFirstSlice(branchID);

    if( lastSliceRadioButton->isChecked() )
        currSlice = skelContourPtr->getLastSlice(branchID);

    if( currSlice == nullptr) return;

    selectedBranchMesh = currSlice->mesh;
    selectedBranchMesh->setName("MeshContours");

    meshViewer->addObject( selectedBranchMesh);
    displayContours(selectedBranchMesh);
}
///////////////////////////////////////////////////////////////////////////////KE
void JMeshSkeletonEditingDialog :: getTwoBranchContours(int branchID)
{
    selectedBranchMesh = JMesh::newObject();
    selectedBranchMesh->setName("MeshContours");

    JSlicePtr slice1  = skelContourPtr->getFirstSlice(branchID);
    selectedBranchMesh->addObject( slice1->mesh);

    JSlicePtr slice2 =  skelContourPtr->getLastSlice(branchID);
    selectedBranchMesh->addObject( slice2->mesh);

    meshViewer->addObject( selectedBranchMesh);
    displayContours(selectedBranchMesh);
}
///////////////////////////////////////////////////////////////////////////////KE

void JMeshSkeletonEditingDialog :: getAllBranchContours(int branchID)
{
    vector<JSlicePtr> allSlices = skelContourPtr->getSlices(branchID);

    selectedBranchMesh = JMesh::newObject();
    JNodeSequence nodes;

    for( const JSlicePtr &slice : allSlices) {
        selectedBranchMesh->addObject(slice->mesh);
    }
    selectedBranchMesh->setName("MeshContours");

    meshViewer->addObject( selectedBranchMesh);
    displayContours(selectedBranchMesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonEditingDialog :: deleteContour()
{
    if( currSlice == nullptr) return;
    currSlice->active = 0;
    meshViewer->removeObject( selectedBranchMesh );
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonEditingDialog :: enumNodes()
{
    if( skelGraph == nullptr) return;

    JMeshRenderPtr mrender;
    skelGraph->getAttribute("Render", mrender);
    mrender->displayIDs[0] =  enumNodesCheckBox->isChecked();
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonEditingDialog :: closeDialog()
{
    JFaceRenderPtr fAttrib;
    size_t numFaces = surfMesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = surfMesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display = 1;
    }
    meshViewer->updateBuffers(surfMesh);
    parentWidget()->show();
    this->close();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonEditingDialog :: makeConnections()
{
    SpinBoxi( selectBranchSpinBox, [=] {selectBranch(); });
    CheckBox( enumNodesCheckBox,   [=] {enumNodes(); });
    CheckBox( displaySurfMeshCheckBox,   [=] {displayBranch(); });
    CheckBox( displayAllEndSlicesCheckBox,   [=] {displayBranch(); });

    PushButton( getContourPushButton,  [=] {getContours(); });
    PushButton( deleteContourPushButton,  [=] {deleteContour(); });

    PushButton( getCapPushButton,  [=] {getCap(); });
    PushButton( getPolesPushButton,[=] {getPoles(); });
    PushButton( closePushButton,   [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
