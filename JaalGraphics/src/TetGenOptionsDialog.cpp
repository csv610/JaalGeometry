#include "TetGenOptionsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JTetGenOptionsDialog :: JTetGenOptionsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
    options = "qzCVY";
    double maxV = std::numeric_limits<double>::max();
    maxVolLineEdit->setText( QString::number(maxV));
}

///////////////////////////////////////////////////////////////////////////////

JTetGenOptionsDialog :: ~JTetGenOptionsDialog()
{

}
///////////////////////////////////////////////////////////////////////////////

void JTetGenOptionsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
}

///////////////////////////////////////////////////////////////////////////////
void JTetGenOptionsDialog :: setOptions()
{
    options = "qzC";

    ostringstream oss;
    if( opt_a_checkBox->isChecked()) {
        QString qstr = maxVolLineEdit->text();
        oss << "a" << qstr.toDouble();
    }

    if( opt_r_checkBox->isChecked() )
        oss << "r";
    /*
           oss << "q";
           oss << radiusEdgeRatioSpinBox->value();
           oss << "/" << minDihedralAngleSpinBox->value();
    */

    if( opt_Y_checkBox->isChecked() )
        oss << "Y";

    /*
       if( tetgenOptCheckBox->isChecked() )
                oss << "O9/7";
    */
    options += oss.str();
}
///////////////////////////////////////////////////////////////////////////////
void JTetGenOptionsDialog :: makeConnections()
{
    PushButton( applyPushButton, [=] {setOptions();});

    CheckBox( opt_a_checkBox, [=] {setOptions();});
    CheckBox( opt_r_checkBox, [=] {setOptions();});
    CheckBox( opt_w_checkBox, [=] {setOptions();});
    CheckBox( opt_Y_checkBox, [=] {setOptions();});
    CheckBox( opt_q_checkBox, [=] {setOptions();});
    CheckBox( optCheckBox,    [=] {setOptions();});

    PushButton(closePushButton, [=] { close(); });
}

///////////////////////////////////////////////////////////////////////////////
