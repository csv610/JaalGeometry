#include "MeshEntityPickerDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshEntityPickerDialog :: JMeshEntityPickerDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;

    numPickedLineEdit->setText( QString::number(0) );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    viewManager->attach(this);
    entityPicker = meshViewer->getEntityPicker();
    entityPicker->setActive(1);
    selectEntity();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    selectEntity();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: mouseReleaseEvent( QMouseEvent *e)
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: setPickStatus()
{
    bool val = pickStatusCheckBox->isChecked();
//  if( entityPicker ) entityPicker->setStatus(val);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: setSelectionMode()
{
    bool singleSelect = singlePickRadioButton->isChecked();

    int mode;
    if( singleSelect )
        mode = 1;
    else
        mode = 2;

    if( entityPicker ) entityPicker->setMode( mode );
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: setPickNodeColor()
{
    JColor color;
    QColor qcolor = QColorDialog::getColor();
    color[0] = qcolor.red()/255.0;
    color[1] = qcolor.green()/255.0;
    color[2] = qcolor.blue()/255.0;
    color[3] = 1.0;
    if( entityPicker ) entityPicker->setHighlightColor(color, 0);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: setPickEdgeColor()
{
    JColor color;
    QColor qcolor = QColorDialog::getColor();
    color[0] = qcolor.red()/255.0;
    color[1] = qcolor.green()/255.0;
    color[2] = qcolor.blue()/255.0;
    color[3] = 1.0;
    if( entityPicker ) entityPicker->setHighlightColor(color, 1);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: setPickFaceColor()
{
    JColor color;
    QColor qcolor = QColorDialog::getColor();
    color[0] = qcolor.red()/255.0;
    color[1] = qcolor.green()/255.0;
    color[2] = qcolor.blue()/255.0;
    color[3] = 1.0;
    if( entityPicker ) entityPicker->setHighlightColor(color, 2);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: setPickCellColor()
{
    JColor color;
    QColor qcolor = QColorDialog::getColor();
    color[0] = qcolor.red()/255.0;
    color[1] = qcolor.green()/255.0;
    color[2] = qcolor.blue()/255.0;
    if( entityPicker ) entityPicker->setHighlightColor(color, 3);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: numPicked()
{
    if( entityPicker == nullptr) return;
    int nSize = entityPicker->getSize(pickable);
    numPickedLineEdit->setText( QString::number(nSize) );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: selectEntity()
{
    if( mesh == nullptr) return;

    pickable = -1;
    if( pickNodeRadioButton->isChecked() ) pickable = 0;
    if( pickEdgeRadioButton->isChecked() ) pickable = 1;
    if( pickFaceRadioButton->isChecked() ) pickable = 2;
    if( pickCellRadioButton->isChecked() ) pickable = 3;

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = pickable;

    numPicked();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: setMark()
{
    /*
        bool val = ordinaryMarkRadioButton->isChecked();
        if( val ) return;

        int entity = entityPicker->getPickableEntity();

        size_t  nSize;
        if( entity == 0) {
            JNodeSequence nodes = entityPicker->getPickedNodes();
            val   = featureMarkRadioButton->isChecked();
            if( val ) {
                nSize = nodes.size();
                for( size_t i = 0; i < nSize; i++)
                    nodes[i]->setAttribute("Picked", 1.0);
            }
            return;
        }

        if( entity == 1) {
            JEdgeSequence edges = entityPicker->getPickedEdges();
            val   = featureMarkRadioButton->isChecked();
            if( val ) {
                nSize = edges.size();
                for( size_t i = 0; i < nSize; i++)
                    edges[i]->setAttribute("Picked", 1.0);
            }
            return;
        }

        if( entity == 2) {
            JFaceSequence faces = entityPicker->getPickedFaces();
            val   = featureMarkRadioButton->isChecked();
            if( val ) {
                nSize = faces.size();
                for( size_t i = 0; i < nSize; i++)
                    faces[i]->setAttribute("Picked", 1.0);
            }
            return;
        }
    */
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: clearAll()
{
    if( entityPicker == nullptr) return;
    entityPicker->clearAll();
    int nSize = entityPicker->getSize(pickable);
    numPickedLineEdit->setText( QString::number(nSize) );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: displaySelected()
{
    bool val;
    size_t nSize = 0;
    JNodeSequence nseq;
    JEdgeSequence eseq;
    JFaceSequence fseq;
    JCellSequence cseq;

    JNodeRenderPtr nAttrib;
    JEdgeRenderPtr eAttrib;
    JFaceRenderPtr fAttrib;

    if( selectOp == SELECT_SHOW ) {
        nSize = mesh->getSize(0);
        val = 0;
        for( size_t i = 0; i < nSize; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            v->getAttribute("Render", nAttrib);
            nAttrib->display = val;
            nAttrib->glyph   = 0;
        }
        val = 1;
        nseq = entityPicker->getPickedNodes();
        for( size_t i = 0; i < nseq.size(); i++) {
            const JNodePtr &v = nseq[i];
            v->getAttribute("Render", nAttrib);
            nAttrib->display = val;
            nAttrib->glyph   = 1;
        }

        nSize = mesh->getSize(1);
        val = 0;
        for( size_t i = 0; i < nSize; i++) {
            const JEdgePtr &e = mesh->getEdgeAt(i);
            e->getAttribute("Render",eAttrib);
            eAttrib->display =  val;
        }
        val = 1;
        eseq = entityPicker->getPickedEdges();
        for( size_t i = 0; i < eseq.size(); i++) {
            const JEdgePtr &e = eseq[i];
            e->getAttribute("Render",eAttrib);
            eAttrib->display =  val;
        }

        nSize = mesh->getSize(2);
        val = 0;
        for( size_t i = 0; i < nSize; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            f->getAttribute("Render",fAttrib);
            fAttrib->display =  val;
        }
        val = 1;
        fseq = entityPicker->getPickedFaces();
        for( size_t i = 0; i < fseq.size(); i++) {
            const JFacePtr &f = fseq[i];
            f->getAttribute("Render",fAttrib);
            fAttrib->display =  val;
        }

        nSize = mesh->getSize(3);
        val = 0;
        for( size_t i = 0; i < nSize; i++) {
            const JCellPtr &c = mesh->getCellAt(i);
            c->setAttribute("Display", val);
        }
        val = 1;
        cseq = entityPicker->getPickedCells();
        for( size_t i = 0; i < nseq.size(); i++) {
            const JCellPtr &c = cseq[i];
            c->setAttribute("Display", val);
        }
    }
    /*

        if( displayAllRadioButton->isChecked() ) {
            val = 1;

            nSize = mesh->getSize(0);
            for( size_t i = 0; i < nSize; i++) {
                const JNodePtr &v = mesh->getNodeAt(i);
                v->getAttribute("Render", nAttrib);
                nAttrib->display = val;
            }

            nSize = mesh->getSize(1);
            for( size_t i = 0; i < nSize; i++) {
                const JEdgePtr &e = mesh->getEdgeAt(i);
                e->getAttribute("Render", eAttrib);
                eAttrib->display = val;
            }

            nSize = mesh->getSize(2);
            for( size_t i = 0; i < nSize; i++) {
                const JFacePtr &f = mesh->getFaceAt(i);
                f->getAttribute("Render", fAttrib);
                fAttrib->display = val;
            }

            nSize = mesh->getSize(3);
            for( size_t i = 0; i < nSize; i++) {
                const JCellPtr &c = mesh->getCellAt(i);
                c->setAttribute("Display", val);
            }
        }

        if( hidePickedRadioButton->isChecked() ) {
            nSize = mesh->getSize(0);
            val = 1;
            for( size_t i = 0; i < nSize; i++) {
                const JNodePtr &v = mesh->getNodeAt(i);
                v->getAttribute("Render", nAttrib);
                nAttrib->display = val;
            }

            val = 0;
            nseq = entityPicker->getPickedNodes();
            for( size_t i = 0; i < nseq.size(); i++) {
                const JNodePtr &v = nseq[i];
                v->getAttribute("Render", nAttrib);
                nAttrib->display = val;
            }

            nSize = mesh->getSize(1);
            val = 1;
            for( size_t i = 0; i < nSize; i++) {
                const JEdgePtr &e = mesh->getEdgeAt(i);
                e->getAttribute("Render", eAttrib);
                eAttrib->display = val;
            }
            val = 0;
            eseq = entityPicker->getPickedEdges();
            for( size_t i = 0; i < eseq.size(); i++) {
                const JEdgePtr &e = eseq[i];
                e->getAttribute("Render", eAttrib);
                eAttrib->display = val;
            }

            nSize = mesh->getSize(2);
            val = 1;
            for( size_t i = 0; i < nSize; i++) {
                const JFacePtr &f = mesh->getFaceAt(i);
                f->getAttribute("Render", fAttrib);
                fAttrib->display = val;
            }
            val = 0;
            fseq = entityPicker->getPickedFaces();
            for( size_t i = 0; i < fseq.size(); i++) {
                const JFacePtr &f = fseq[i];
                f->getAttribute("Render", fAttrib);
                fAttrib->display = val;
            }

            nSize = mesh->getSize(3);
            val = 1;
            for( size_t i = 0; i < nSize; i++) {
                const JCellPtr &c = mesh->getCellAt(i);
                c->setAttribute("Display", val);
            }
            val = 0;
            cseq = entityPicker->getPickedCells();
            for( size_t i = 0; i < nseq.size(); i++) {
                const JCellPtr &c = cseq[i];
                c->setAttribute("Display", val);
            }
        }

        if( meshViewer ) meshViewer->refreshDisplay();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: deleteSelected()
{
    /*
        if( entityPicker == nullptr ) return;
        int entity = entityPicker->getPickableEntity();

        size_t  nSize;
        if( entity == 0) {
            JNodeSequence nodes = entityPicker->getPickedNodes();
            nSize = nodes.size();
            for( size_t i = 0; i < nSize; i++)
                nodes[i]->setStatus(JMeshEntity::REMOVE);
        }

        if( entity == 1) {
            JEdgeSequence edges = entityPicker->getPickedEdges();
            nSize = edges.size();
            for( size_t i = 0; i < nSize; i++)
                edges[i]->setStatus( JMeshEntity::REMOVE);
        }

        if( entity == 2) {
            JFaceSequence faces = entityPicker->getPickedFaces();
            nSize = faces.size();
            for( size_t i = 0; i < nSize; i++)
                faces[i]->setStatus( JMeshEntity::REMOVE );
        }

        if( entity == 3) {
            JCellSequence cells = entityPicker->getPickedCells();
            nSize = cells.size();
            for( size_t i = 0; i < nSize; i++)
                cells[i]->setStatus( JMeshEntity::REMOVE );
        }
        if( meshViewer) meshViewer->refreshDisplay();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: reject()
{
    /*
         if( meshViewer ) {
              entityPicker->setStatus(0);
              meshViewer->displayAll(3,1);
              clearAll();
         }
    */
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: setPickWidth()
{
    /*
        if( viewManager == nullptr ) return;
        int w = pickWidthSpinBox->value();
        viewManager->setSelectRegionWidth( w );
    */
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: setPickHeight()
{

    /*
        if( viewManager == nullptr ) return;
        int h = pickHeightSpinBox->value();
        viewManager->setSelectRegionHeight( h );
    */
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: addManually()
{
    /*
        QString qstr;
        qstr = manualLineEdit->text();
        string str = qstr.toUtf8().constData();

        StringTokenizer  token(str);
        StringTokenizer::iterator iter;
        for(iter = token.begin(); iter != token.end(); ++iter) {
            string s = *iter;
            int id = atoi(s.c_str());
            if( entityPicker) entityPicker->append(id);
        }
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: selectAll()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: hideSelected()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: growSelected()
{
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEntityPickerDialog :: shrinkSelected()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: setSensitivity()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: setLensRadius()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: invertSelected()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: closeDialog()
{
    if( entityPicker ) entityPicker->setActive(0);
    viewManager->detach(this);
    this->close();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEntityPickerDialog :: makeConnections()
{
    PushButton( pickNodeColorPushButton, [=] { setPickNodeColor();});
    PushButton( pickEdgeColorPushButton, [=] { setPickEdgeColor();});
    PushButton( pickFaceColorPushButton, [=] { setPickFaceColor();});
    PushButton( pickCellColorPushButton, [=] { setPickCellColor();});
    PushButton( clearAllPickPushButton,  [=] { clearAll();});
    PushButton( numPickedPushButton,     [=] { numPicked();});

    RadioButton( pickNodeRadioButton,   [=] {selectEntity();});
    RadioButton( pickEdgeRadioButton,   [=] {selectEntity();});
    RadioButton( pickFaceRadioButton,   [=] {selectEntity();});
    RadioButton( pickCellRadioButton,   [=] {selectEntity();});
    RadioButton( singlePickRadioButton, [=] {setSelectionMode();});
    RadioButton( multiPickRadioButton,  [=] {setSelectionMode();});

    CheckBox( pickStatusCheckBox, [=] {setPickStatus();});
    PushButton( closePushButton,  [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
