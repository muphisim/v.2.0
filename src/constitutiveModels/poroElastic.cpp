//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#include "poroElastic.h"
#include "maths.h"
#include "classGPs.h"

DarcyLaw::DarcyLaw(int ndim, string name, string name_cons, double kappa, double mu, double Q,  double alpha ):
	constitutiveModels(ndim, name, name_cons, 0.), _kappa(kappa), _mu(mu), _Q(Q), _alpha(alpha), _intVarsSize(0)
	{
		INFO("create a DarcyLaw with _kappa(%g), _mu(%g), _Q(%g), _alpha(%g)",_kappa, _mu, _Q, _alpha);
	}
	
void DarcyLaw::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    _term = new ExtraDofMechanicalFullCouplingTerm(_extraFieldIndexes[0]);
};

double DarcyLaw::soundSpeed() const {
    return 0;
}

void DarcyLaw::initIntVars(vector<double> &intVars) 
{
  _intVarsSize = intVars.size();
  intVars.push_back(0);
};


/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void DarcyLaw::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &FCurr = currState.getDeformationGradient();
    
    double p = currState.getExtraFields()[_extraFieldIndexes[0]];
    double pprev = prevState.getExtraFields()[_extraFieldIndexes[0]];
    const classTensor1& gradp = currState.getExtraFieldGradients()[_extraFieldIndexes[0]];
    
    // temperature -dependent of conductivity
    
    int ndim = getDim();

    classTensor1& flux = currState.getExtraFieldFluxes()[_extraFieldIndexes[0]];
    double& source = currState.getExtraFieldSources()[_extraFieldIndexes[0]];
    double& fieldDensity = currState.getExtraFieldDensities()[_extraFieldIndexes[0]];
    
    fieldDensity = 1./_Q;
    flux = gradp;
    flux.scale(-_kappa/_mu); // Q=-kappa/mu*gradp
    
    // source must be added up to existing one
    source += fieldDensity*(p-pprev)/dt;
    
    // update 
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    // liquid volume fraction
    double trE=0;
    for (int i=0; i< ndim; i++)
    {
        trE += FCurr[i * ndim + i]-1.;
    }
    intVarsCurr[_intVarsSize] = _alpha*trE + p/_Q;
    
    if (flagTanMod)
    {
        double& DsourceDp = (*currState.getDextraFieldSourcesDextraFields())[_extraFieldIndexes[0]][_extraFieldIndexes[0]];
        classTensor1& DfluxDp = (*currState.getDextraFieldFluxesDextraFields())[_extraFieldIndexes[0]][_extraFieldIndexes[0]];
        classTensor2& DfluxDgradp = (*currState.getDextraFieldFluxesDextraFieldGrads())[_extraFieldIndexes[0]][_extraFieldIndexes[0]];
        
        // add-up
        DsourceDp += fieldDensity/dt;
        DfluxDp.setAll(0.);
        DfluxDgradp.setAll(0.);
        DfluxDgradp.addDiagonal(-_kappa/_mu);
    }
};


LinearPoroelastic::LinearPoroelastic(int ndim, string name, string name_cons, double rho, double young, double nu, double alpha, int stressState):
	linearElastic(ndim, name, name_cons, rho, young, nu, stressState), _alpha(alpha){
		
		INFO("create a LinearPoroelastic with young %g  nu %g alpha %g",young, nu, alpha);
	}
	
void LinearPoroelastic::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term; 
    _term = new MechanicalExtraDofFullCouplingTermPoroElastic(std::vector<int>(1,0), _alpha);
};

void LinearPoroelastic::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
	linearElastic::updateConstitutive(GaussPoint,dt,timeRun,flagTanMod);
	// update mechanical source
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    
    int ndim = getDim();
    double& mecaSource = currState.getExtraFieldSources()[0];
    for (int i=0; i<ndim; i++)
    {
      mecaSource += _alpha* (FCurr[i * ndim + i]-FPrev[i * ndim + i])/dt; //
    }
    
    if (flagTanMod)
    {
        classTensor2& DsourceDF = (*currState.getDextraFieldSourcesDdeformationGradient())[0];        
        for (int i = 0; i < ndim; i++) 
        {
            DsourceDF(i,i) = _alpha/dt;
        }
    }
};