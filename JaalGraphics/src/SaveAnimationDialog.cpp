#include "SaveAnimationDialog.hpp"
///////////////////////////////////////////////////////////////////////////////

JSaveAnimationDialog :: JSaveAnimationDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JSaveAnimationDialog :: ~JSaveAnimationDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JSaveAnimationDialog :: makeConnections()
{
    PushButton( closePushButton,  [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
