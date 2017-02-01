#include "LocallyInjectiveMapParamsDialog.hpp"
///////////////////////////////////////////////////////////////////////////////

JLocallyInjectiveMapParamsDialog :: JLocallyInjectiveMapParamsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    limDeformer = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JLocallyInjectiveMapParamsDialog :: ~JLocallyInjectiveMapParamsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JLocallyInjectiveMapParamsDialog :: init()
{
    if( limDeformer == nullptr ) return;

    QString qs = energyComboBox->currentText();
    string str = qs.toUtf8().constData();

    /*
        // What parameters the solver holding before the solver starts ...
        double alpha = limDeformer->getAlpha();
        alphaLineEdit->setText( QString::number(alpha) );
        double alphaRatio = limDeformer->getAlphaRatio();
        alphaRatioLineEdit->setText( QString::number(alphaRatio) );
    */

    double beta = limDeformer->getBeta();
    betaLineEdit->setText( QString::number(beta) );

    double gamma = limDeformer->getGamma();
    gammaLineEdit->setText( QString::number(gamma) );


    int maxIter  = limDeformer->getMaxIterations();
    maxIterationsLineEdit->setText( QString::number(maxIter) );

}
///////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMapParamsDialog :: setParams()
{
    if( limDeformer == nullptr ) return;
    ///////////////////////////////////////////////////////////////////////////
    // With parameters, the solver will run ? set all the parameters. These
    // parameters can not modified in the middle of the some run..
    ///////////////////////////////////////////////////////////////////////////
    QString qstr = energyComboBox->currentText();
    string str = qstr.toUtf8().constData();

    if( str == "AsRigidAsPossible")
        limDeformer->setEnergyType( JLocallyInjectiveMap::AS_RIGID_AS_POSSIBLE_ENERGY);

    if( str == "Laplacian")
        limDeformer->setEnergyType( JLocallyInjectiveMap::LAPLACIAN_ENERGY);

    if( str == "Dirichlet")
        limDeformer->setEnergyType( JLocallyInjectiveMap::DIRICHLET_ENERGY);

    if( str == "GreenStrain")
        limDeformer->setEnergyType( JLocallyInjectiveMap::GREEN_STRAIN_ENERGY);

    if( str == "Identity")
        limDeformer->setEnergyType( JLocallyInjectiveMap::IDENTITY_ENERGY);

    if( str == "LeastSquareConformal")
        limDeformer->setEnergyType( JLocallyInjectiveMap::LEAST_SQUARE_CONFORMAL_ENERGY);

    if( str == "Poisson2D")
        limDeformer->setEnergyType( JLocallyInjectiveMap::POISSON_ENERGY);

    if( str == "UniformLaplacian")
        limDeformer->setEnergyType( JLocallyInjectiveMap::UNIFORM_LAPLACE_ENERGY);

    qstr = maxIterationsLineEdit->text() ;
    limDeformer->setMaxIterations(qstr.toInt());

    bool val;
    barriersCheckBox->isChecked();
    limDeformer->setBarriers( val );

    val = logBarrierCheckBox->isChecked();
    limDeformer->setBarrierType(JLocallyInjectiveMap::LOG_BARRIER, val );

    val = neoHookeanBarrierCheckBox->isChecked();
    limDeformer->setBarrierType( JLocallyInjectiveMap::NEOHOOKEAN_BARRIER, val );

    val = barrierCompensationCheckBox->isChecked();
    limDeformer->setBarrierCompensation( val );

    val = alphaUpdateCheckBox->isChecked();
    limDeformer->setAlphaUpdate( val );

    val = substeppingCheckBox->isChecked();
    limDeformer->setSubstepping( val );

    /*
        qstr = alphaLineEdit->text() ;
        limDeformer->setAlpha( qstr.toDouble() );
        qstr = alphaRatioLineEdit->text() ;
        limDeformer->setAlphaRatio( qstr.toDouble() );
    */
    qstr = betaLineEdit->text() ;
    limDeformer->setBeta( qstr.toDouble() );

    qstr = gammaLineEdit->text() ;
    limDeformer->setGamma( qstr.toDouble() );
}


///////////////////////////////////////////////////////////////////////////////

void JLocallyInjectiveMapParamsDialog :: makeConnections()
{
    connect( applyPushButton, SIGNAL( clicked() ), this, SLOT( setParams() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
