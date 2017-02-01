#include "SuggestiveContoursDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
JSuggestiveContoursDialog :: JSuggestiveContoursDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JSuggestiveContoursDialog :: ~JSuggestiveContoursDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JSuggestiveContoursDialog :: init()
{
    if( viewManager == nullptr ) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    scViewer.reset( new JSuggestiveContoursViewer(viewManager) );
    scViewer->setName("SuggestiveContoursViewer");
    viewManager->attach( scViewer );
    setMesh( meshViewer->getCurrentMesh() );
}
///////////////////////////////////////////////////////////////////////////////

void JSuggestiveContoursDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    JMeshOFFExporter mexp;
    mexp.writeFile(mesh, "scmodel.off");
    scViewer->loadNewMesh( "scmodel.off");

    // Let Suggestive Contour display takeover ...
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->display = 0;

    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JSuggestiveContoursDialog :: closeDialog()
{
    if( scViewer ) viewManager->detach(scViewer);

    // Let standard display takeover ...
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->display = 1;

    this->close();
}

///////////////////////////////////////////////////////////////////////////////
void JSuggestiveContoursDialog :: setParam()
{
    if( scViewer == nullptr) return;

    bool val;
    val = drawEdgesCheckBox->isChecked();
    scViewer->draw_edges = val;

    val = exteriorSilhouettesCheckBox->isChecked();
    scViewer->draw_extsil = val;

    val = occludingContoursCheckBox->isChecked();
    scViewer->draw_c = val;

    val = suggestiveContoursCheckBox->isChecked();
    scViewer->draw_sc = val;

    val = suggestiveHighlightsCheckBox->isChecked();
    scViewer->draw_sh = val;

    val = principalHighlightRCheckBox->isChecked();
    scViewer->draw_phridges = val;

    val = principalHighlightVCheckBox->isChecked();
    scViewer->draw_phvalleys = val;

    val = ridgesCheckBox->isChecked();
    scViewer->draw_ridges = val;

    val = apparentRidgesCheckBox->isChecked();
    scViewer->draw_apparent = val;

    val = valleysCheckBox->isChecked();
    scViewer->draw_valleys = val;

    val = drawKCheckBox->isChecked();
    scViewer->draw_K = val;

    val = drawHCheckBox->isChecked();
    scViewer->draw_H = val;

    val = drawDwKrCheckBox->isChecked();
    scViewer->draw_DwKr = val;

    val = drawBoundariesCheckBox->isChecked();
    scViewer->draw_bdy  = val;

    val = drawIsotopesCheckBox->isChecked();
    scViewer->draw_isoph  = val;

    val = drawHiddenCheckBox->isChecked();
    scViewer->draw_hidden  = val;

    val = drawFadedCheckBox->isChecked();
    scViewer->draw_faded  = val;

    val = useTextureCheckBox->isChecked();
    scViewer->use_texture  = val;

    val = drawColorsCheckBox->isChecked();
    scViewer->draw_colors  = val;

    val = hermiteCheckBox->isChecked();
    scViewer->use_hermite  = val;

    int id = lightStyleComboBox->currentIndex();
    scViewer->lighting_style = id;

    viewManager->refreshDisplay();
}

void JSuggestiveContoursDialog :: smoothNormals()
{
    if( scViewer == nullptr) return;
    scViewer->smoothNormals();
    viewManager->refreshDisplay();
}

void JSuggestiveContoursDialog :: smoothCurvature()
{
    if( scViewer == nullptr) return;
    scViewer->smoothCurvature();
    viewManager->refreshDisplay();
}

void JSuggestiveContoursDialog :: smoothCurvatureDeriv()
{
    if( scViewer == nullptr) return;
    scViewer->smoothCurvatureDeriv();
    viewManager->refreshDisplay();
}

void JSuggestiveContoursDialog :: subdivisionSurface()
{
    if( scViewer == nullptr) return;
    scViewer->subdivisionSurface();
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JSuggestiveContoursDialog :: makeConnections()
{
    CheckBox( drawEdgesCheckBox, [=] {setParam(); });
//    CheckBox( ridgesCheckBox,    [=]{setParam();));
    CheckBox( drawKCheckBox,   [=] {setParam(); });
    CheckBox( drawHCheckBox,   [=] {setParam();});
    CheckBox( valleysCheckBox,  [=] {setParam(); });
    CheckBox( drawDwKrCheckBox,   [=] {setParam(); });
    CheckBox( drawIsotopesCheckBox, [=] {setParam(); });
    CheckBox( drawBoundariesCheckBox,  [=] {setParam(); });
    CheckBox( apparentRidgesCheckBox,   [=] {setParam(); });
    CheckBox( occludingContoursCheckBox,   [=] {setParam(); });
    CheckBox( suggestiveContoursCheckBox,   [=] {setParam(); });
    CheckBox( suggestiveHighlightsCheckBox,   [=] {setParam(); });
    CheckBox( principalHighlightRCheckBox,   [=] {setParam(); });
    CheckBox( principalHighlightVCheckBox,   [=] {setParam(); });
    CheckBox( exteriorSilhouettesCheckBox,  [=] {setParam(); });

    ComboBox( lightStyleComboBox, [=] {setParam(); });

    PushButton( smoothNormalsPushButton,  [=] {smoothNormals();});
    PushButton( smoothCurvaturePushButton, [=] {smoothCurvature(); });
    PushButton( smoothCurvatureDerivPushButton, [=] {smoothCurvatureDeriv(); });
    PushButton( subdividePushButton, [=] {subdivisionSurface(); });

    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
