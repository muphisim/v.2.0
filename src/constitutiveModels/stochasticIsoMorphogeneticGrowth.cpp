//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file stochastiIsoMorphogeneticGrowth.cpp
  \brief This file contains all functions related to isotropic growth constitutive model. There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/
#include "stochasticIsoMorphogeneticGrowth.h"
#include "maths.h"
#include "classGPs.h"
/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young Young modulus
  @param[in] nu Poisson ratio
  @param[in] Gc Growth multiplier
*/
classStochasticIsoMorphogeneticGrowth::classStochasticIsoMorphogeneticGrowth(int ndim, string name, string name_cons, double rho,
                                                         double Young, double nu, double Gc, vector<double> Stochastic_numbers,
                                                     vector<double> Stochastic_parameters,
                                                     vector<double> Stochastic_function, string approximation, int order, int resolution, bool stochasticMapping, string distribution) : 
                                                         
            constitutiveModels(ndim, name,  name_cons, rho){


    int number_func = 0;
    for (int j = 0; j < Stochastic_numbers.size(); j++) {
        number_func = number_func + Stochastic_numbers[j];
    }
    int number_stochastic = 0;
    if (approximation == "Haar") {
        number_stochastic = binomialCoefficients(pow(2, order + 1) + number_func - 1, number_func);
    } else {
        if (approximation == "Mixed") {
        number_stochastic = binomialCoefficients((order + 1)*(2*pow(2,resolution)) + number_func - 1, number_func);
        } else {
            number_stochastic = binomialCoefficients(order + number_func, order);
        }
    }
    this->_stochasticMapping = stochasticMapping;
    vector<double> zeros(number_stochastic, 0);
    this->elasticStiffTensor = new classTensor6(ndim, number_stochastic);
    this->_Young = zeros;
    this->_nu = zeros;
    this->_Gc = zeros;
    this->_Young[0] = Young;
    this->_nu[0] = nu;
    this->_Gc[0] = Gc;
    int offset = 1;
    for (int j = 0; j < Stochastic_parameters.size(); j++) {
        if (Stochastic_parameters[j] == 1) {
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->_Young[offset + i] = Young * Stochastic_function[i];
                //count++;
            }
            offset = offset + Stochastic_numbers[j];
        }
        if (Stochastic_parameters[j] == 2) {
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->_nu[offset + i] = nu * Stochastic_function[i + offset - 1];
            }
            offset = offset + Stochastic_numbers[j];
        }
        if (Stochastic_parameters[j] == 3) {
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->_Gc[offset + i] = Gc * Stochastic_function[i + offset - 1];
            }
            offset = offset + Stochastic_numbers[j];
        }
    }

    
    vector<double> C(number_stochastic * number_stochastic * number_stochastic, 0);
    if (approximation == "Haar") {
        Build_C_Haar(order, number_func, C);
        ConvertToHaar(this->_Young, this->_nu, order, Stochastic_parameters, Stochastic_numbers, distribution);
    } else {
        if (approximation == "Mixed"){
            Build_C_Mixed(order, resolution, number_func, C, distribution);
        } else {
 
            build_ThirdOrderTensor_C(order, number_func, C, distribution);
        }
    }
    this->C_Stoch = C;
    this->nstoch = number_stochastic;

        vector<double> lambda(number_stochastic, 0);
    vector<double> mu(number_stochastic, 0);
    vector<double> Young_times_nu(number_stochastic, 0);
    vector<double> Young_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> oneplusnu_times_oneminustwonu(number_stochastic, 0);
    vector<double> oneminustwicenu(number_stochastic, 0);
    vector<double> factor(number_stochastic, 0);
    vector<double> factor_times_nu(number_stochastic, 0);
    vector<double> young_over_oneminustwonu(number_stochastic, 0);
    oneminustwicenu = this->_nu;
    vector<double> oneplusnu = this->_nu;
    vector<double> __nu = this->_nu;
    vector<double> __Young = this->_Young; 
    linearCombinationVariableStochastic(oneplusnu, 1, 1);
    linearCombinationVariableStochastic(oneminustwicenu, -2, 1);


    DivisionRandomVariableRandomVariableStochastic(__Young, oneplusnu, C_Stoch, mu);
    vector<double> Young_over_oneplusnu = mu;
    MultiplicationRandomVariableScalar(mu, 0.5);

   
    DivisionRandomVariableRandomVariableStochastic(Young_over_oneplusnu, oneminustwicenu, C_Stoch, factor);
    multiRandomVariableRandomVariableStochastic(factor, __nu, C_Stoch, factor_times_nu);
    lambda = factor_times_nu;
    DivisionRandomVariableRandomVariableStochastic(__Young, oneminustwicenu, C_Stoch, young_over_oneminustwonu);
    this->_lambda = lambda;
    this->_mu = mu;

    vector<double> dmufdf(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    
    for(int alpha =0; alpha < nstoch; alpha++){
        for(int gamma =0; gamma < nstoch; gamma++){
            for(int delta =0; delta < nstoch; delta++){
                for(int a =0; a < ndim; a++){
                    for(int b =0; b < ndim; b++){
                        dmufdf [gamma*ndim*ndim*ndim*ndim*nstoch + delta*ndim*ndim*ndim*ndim + b*ndim*ndim*ndim + a*ndim*ndim + b*ndim + a] += C[gamma*nstoch*nstoch + delta*nstoch + alpha]*mu[alpha];
                    }
                }
            }
        }
    }
    this->dmuFdF = dmufdf;

};

classStochasticIsoMorphogeneticGrowth::~classStochasticIsoMorphogeneticGrowth()
{
 
}

void classStochasticIsoMorphogeneticGrowth::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
     _term = new MechanicalTermStochastic();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classStochasticIsoMorphogeneticGrowth::soundSpeed() const 
{
    double rho = this->getRho();
    double sound = sqrt(_Young[0] / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classStochasticIsoMorphogeneticGrowth::initIntVars(vector<double> &intVars) {
    intVars.resize(nstoch,0);
    intVars[0] = 1.0;

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classStochasticIsoMorphogeneticGrowth::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    classStates& initState = GaussPoint->getInitialState();
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &FIni = initState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    vector<double> PK2(ndim *ndim *nstoch, 0);
    vector<double>FIniInv(nstoch*ndim*ndim,0);
    vector<double>FCurrMecha(nstoch*ndim*ndim,0);
    vector<double>det_FIni(nstoch,0);
    if(_stochasticMapping) {
        InverseTensorStochastic(FIni, ndim, nstoch, C_Stoch, FIniInv);
        multTensorTensor3Stochastic(FCurr, ndim, ndim, FIniInv, ndim, FCurrMecha, nstoch, C_Stoch);
        DeterminantStochasticTensor(FIni, ndim, nstoch, det_FIni, C_Stoch);
    }else {
            FCurrMecha = FCurr;
    }

    computeGLstoch(FCurrMecha, ndim, nstoch, C_Stoch, E);
    //computeEquivalentStrainStochastic(E, ndim, nstoch, C_Stoch, equivalentStrain);
    computeVolumetricStrainStochastic(E, ndim, nstoch, volumetricStrain);
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();

    vector<double> theta = intVarsPrev;
    if (timeRun > 1E-16) {
        for(int i=0;i<nstoch;i++){
        theta[i] += _Gc[i] * dt;
        }
    } else {
        theta = intVarsPrev;
    }
    intVarsCurr = theta;

    vector<double> Jacobian(nstoch,0);
    vector<double> Je(nstoch,0);
    DeterminantStochasticTensor(FCurr, ndim, nstoch, Jacobian, C_Stoch);
    
    stress(PK2, PK1, FCurr, theta);

    vector<double>FIniInvT(nstoch*ndim*ndim,0);
    vector<double>PK2_FIniInvT(nstoch*ndim*ndim,0);
    vector<double>FIniInv_PK2_FIniInvT(nstoch*ndim*ndim,0);
    if(_stochasticMapping){
        transposeStochastic(FIniInv, ndim, nstoch, FIniInvT);
        multTensorTensor3Stochastic(PK2, ndim, ndim, FIniInvT, ndim, PK2_FIniInvT, nstoch, C_Stoch);
        multTensorTensor3Stochastic(FIniInv, ndim, ndim, PK2_FIniInvT, ndim, FIniInv_PK2_FIniInvT, nstoch, C_Stoch);
        multiSTensorRandomVariableStochastic(FIniInv_PK2_FIniInvT, ndim, ndim, det_FIni, nstoch, C_Stoch, PK2);
    }

    vector<double> dmuFinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dlambdalogJFinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dmuthetapowerdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);   
    double tol = 0.0000001;
    build_derivatives(FCurr, dmuFinvTdF, dlambdalogJFinvTdF, dmuthetapowerdF, theta, C_Stoch, tol, nstoch, ndim);
    
    double val = 0;
    computeCauchyStoch(FCurr, ndim, nstoch, C_Stoch, PK1, Cauchy);
    //computeStochasticVMS(Cauchy, ndim, nstoch, VMS, C_Stoch);
    if (flagTanMod) {
        build_derivatives(FCurr, dmuFinvTdF, dlambdalogJFinvTdF, dmuthetapowerdF, theta, C_Stoch, tol, nstoch, ndim);
        classTensor6* tanModuli = currState.getTangentStochastic();
        for (int iota = 0; iota < nstoch; iota++) {
            for (int i = 0; i < ndim; i++) {
                for (int J = 0; J < ndim; J++) {
                    for (int eps = 0; eps < nstoch; eps++) {
                        for (int k = 0; k < ndim; k++) {
                            for (int M = 0; M < ndim; M++) {
                                val = dmuthetapowerdF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J] - dmuFinvTdF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J] + dlambdalogJFinvTdF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J];
                                tanModuli->setValues(eps, iota, k, M, i, J, val);
                            }
                        }
                    }
                }
            }
        }
    }

}


/*! \brief This function provides the growth tensor. It would used for external calls, as the growth part is unstressed in the deformation gradient decomposition.
  @param[in] FCurr Current deformation gradient tensor
  @param[in] FPrev Current deformation gradient tensor
  @param[in] intVarsPrev Array of the previous state of the internal variables
  @param[out] intVarsCurr Array of the current state of the internal variables
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
  @param[in] PK1 First Piola Kirchhoff stress tensor
  @param[in] tanModuli Tangent moduli
  @param[out] Fg Growth tensor
*/
void classStochasticIsoMorphogeneticGrowth::getFg(const vector<double> &FCurr, const vector<double> &FPrev, const vector<double> &intVarsPrev,
                                        vector<double> &intVarsCurr, double dt, double timeRun, bool flagTanMod,
                                        const vector<double> &PK1, classTensor4 *tanModuli, vector<double> &Fg) const {
    int ndim = getDim();
    vector<double> theta_power(nstoch,0);
    vector<double> theta = intVarsPrev;
    if (timeRun > 1E-16) {
        for(int i=0;i<nstoch;i++){
        theta[i] += _Gc[i] * dt;
        }
    } else {
        theta = intVarsPrev;
    }
    intVarsCurr = theta; // update internal var

    double ndime = double(ndim + 1E-10);
    nthPowerIntegral(theta, theta_power, C_Stoch, (1.0/ndime));
    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }
    for (int alpha = 0; alpha < ndim; alpha++) {
        for (int i = 0; i < ndim; i++) {
            for (int j = 0; j < ndim; j++) {
                Fg[alpha*ndim*ndim + i * ndim + j] = theta_power[alpha] * delta[i * ndim + j];
            }
        }
    }
}

/*! \brief This function provides the elastic Green-Lagrangian strain tensor. Please not that this E is not the same than the one provided in the results of MuPhiSim (that one is calculated with the whole deformation gradient) 
  @param[in] ndim Dimension of the domain
  @param[out] E Green-Lagrangian strain tensor
  @param[in] FCurr Current deformation gradient tensor
  @param[in] FPrev Previous deformation gradient tensor
  @param[in] intVarsPrev Array of the internal variables
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] nStep Current number of time steps
*/
void classStochasticIsoMorphogeneticGrowth::getE(int ndim, vector<double> &E, const vector<double> &FCurr, const vector<double> &FPrev,
                                       const vector<double> &intVarsPrev, double dt, double timeRun, int nStep) const{
    vector<double> theta = intVarsPrev;
    vector<double> theta_power(nstoch,0);
    vector<double> Fe(nstoch*ndim * ndim, 0);
    double ndime = double(ndim + 1E-10);
    nthPowerIntegral(theta, theta_power, C_Stoch, (-1.0/ndime));

    multiSTensorRandomVariableStochastic(FCurr, ndim, ndim, theta_power, nstoch, C_Stoch,
                                                 Fe);
    computeGLstoch(Fe, ndim, nstoch, C_Stoch, E);

};

void classStochasticIsoMorphogeneticGrowth::stress(vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr, const vector<double> &Theta) {
    int number_stochastic = nstoch;
    int ndim = getDim();
    setAll(PK1,0);
    vector<double> FInv(ndim *ndim *number_stochastic, 0);
    vector<double> FInvT(ndim *ndim *number_stochastic, 0);
    vector<double> Jacobian(number_stochastic,0);
    vector<double> Je(number_stochastic,0);
    vector<double> theta_power(number_stochastic,0);
    vector<double> logJe(nstoch,0);
    vector<double> Jacobian_over_theta(nstoch,0);
    vector<double> lambda_logJe(nstoch,0);
    vector<double> lambda_logJ_over_theta(nstoch,0);
    vector<double> factor1(ndim *ndim *number_stochastic,0);
    vector<double> factor2(ndim *ndim *number_stochastic,0);
    vector<double> factor3(ndim *ndim *number_stochastic,0);
    vector<double> firstPart(ndim *ndim *number_stochastic,0);
    vector<double> Jacobian_power(number_stochastic,0);
    vector<double> mu_over_theta_power(number_stochastic,0);
    vector<double> ones(number_stochastic,0);
    vector<double> one_over_theta(number_stochastic,0);
    ones[0] =1;
    DivisionRandomVariableRandomVariableStochastic(ones, Theta, C_Stoch, one_over_theta);
    
    double power = -2.0/(double(ndim));
     

    DeterminantStochasticTensor(FCurr, ndim, number_stochastic, Jacobian, C_Stoch);
    multiRandomVariableRandomVariableStochastic(Jacobian, one_over_theta, C_Stoch,
                                                 Je);
    logIntegral(Je, logJe, C_Stoch);
    
    multiRandomVariableRandomVariableStochastic(_lambda, logJe, C_Stoch,
                                                 lambda_logJe); 
    
    InverseTensorStochastic(FCurr, ndim, number_stochastic, C_Stoch, FInv);
    transposeStochastic(FInv, ndim, number_stochastic, FInvT);
    
    nthPowerIntegral(Theta, theta_power, C_Stoch, power);

    multiRandomVariableRandomVariableStochastic(_mu, theta_power, C_Stoch,
                                                 mu_over_theta_power);

    multiSTensorRandomVariableStochastic(FCurr, ndim, ndim, mu_over_theta_power, number_stochastic, C_Stoch,
                                                 factor1);                                        
    multiSTensorRandomVariableStochastic(FInvT, ndim, ndim, lambda_logJe, number_stochastic, C_Stoch,
                                                 factor2);
    multiSTensorRandomVariableStochastic(FInvT, ndim, ndim, _mu, number_stochastic, C_Stoch,
                                                 factor3);
                                                                                                                                     
    

    
    
    for (int i = 0; i < number_stochastic; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                PK1[i*ndim*ndim + j*ndim +k] = factor1[i*ndim*ndim + j*ndim +k] + factor2[i*ndim*ndim + j*ndim +k] - factor3[i*ndim*ndim + j*ndim +k];
            }
        }
    }

    multTensorTensor3Stochastic(FInv, ndim, ndim, PK1, ndim, PK2, nstoch, C_Stoch);
    

}

void classStochasticIsoMorphogeneticGrowth::build_derivatives(const vector<double> &FCurr, vector<double> &dmuFinvTdF, vector<double> &dlambdalogJeFinvTdF, vector<double> &dmuthetapowerdF, vector<double> &theta, vector<double> &C_Stoch, const double tol, const int nstoch, const int ndim){

    setAll(dmuFinvTdF,0);
    setAll(dlambdalogJeFinvTdF,0);

    vector<double> lambda_logJe(nstoch,0);
    vector<double> Je(nstoch,0);
    vector<double> Jacobian(nstoch,0);
    vector<double> Jacobian_Perturbed(nstoch,0);
    vector<double> FPerturbed(ndim*ndim*nstoch,0);
    vector<double> FInv(ndim *ndim *nstoch, 0);
    vector<double> FInvT(ndim *ndim *nstoch, 0);
    vector<double> FInvPerturbed(ndim *ndim *nstoch, 0);
    vector<double> FInvTPerturbed(ndim *ndim *nstoch, 0);
    vector<double> logJ(nstoch,0);
    vector<double> logJe(nstoch,0);
    vector<double> dlogJedF(nstoch*nstoch*ndim*ndim,0);
    vector<double> dlambdalogJedF(nstoch*nstoch*ndim*ndim,0);
    vector<double> one_overJe(nstoch,0);
    vector<double> ones(nstoch,0);
    vector<double> Cmu(nstoch*nstoch,0);
    ones[0] = 1;
    vector<double> dJdF(nstoch*nstoch*ndim*ndim,0);
    vector<double> dFinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dFdFFinv(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> CFinv(nstoch*nstoch*ndim*ndim,0);
    vector<double> theta_power(nstoch,0);
    vector<double> dtheta_powerdF(nstoch*nstoch*ndim*ndim,0);
    
    vector<double> muF(nstoch*ndim*ndim,0);
    vector<double> one_over_theta(nstoch,0);
    DivisionRandomVariableRandomVariableStochastic(ones, theta, C_Stoch, one_over_theta);

    double power = -2.0/(double(ndim));

    // Compute F-T
    InverseTensorStochastic(FCurr, ndim, nstoch, C_Stoch, FInv);
    transposeStochastic(FInv, ndim, nstoch, FInvT);

    FPerturbed = FCurr;
    DeterminantStochasticTensor(FCurr, ndim, nstoch, Jacobian, C_Stoch);
    multiRandomVariableRandomVariableStochastic(Jacobian, one_over_theta, C_Stoch,
                                                 Je);
    logIntegral(Je, logJe, C_Stoch);
    multiRandomVariableRandomVariableStochastic(_lambda, logJe, C_Stoch, lambda_logJe);
    nthPowerIntegral(theta, theta_power, C_Stoch, power);
    multiSTensorRandomVariableStochastic(FCurr, ndim, ndim, _mu, nstoch, C_Stoch,
                                                 muF);


    for(int beta =0; beta < nstoch; beta++){
        for(int b =0; b < ndim; b++){
            for(int a =0; a < ndim; a++){
                FPerturbed = FCurr;
                FPerturbed[beta*ndim*ndim + b*ndim +a ] += tol;
                DeterminantStochasticTensor(FPerturbed, ndim, nstoch, Jacobian_Perturbed, C_Stoch);
                for(int alpha =0; alpha < nstoch; alpha++){
                    dJdF[alpha*ndim*ndim*nstoch + beta*ndim*ndim + b*ndim +a] = (Jacobian_Perturbed[alpha] - Jacobian[alpha])/tol;
                }
            }
        }
    }

    for(int gamma =0; gamma < nstoch; gamma++){
        for(int delta =0; delta < nstoch;  delta++){
            for(int c =0; c < ndim; c++){
                for(int d =0; d < ndim; d++){
                    for(int beta =0; beta < nstoch; beta++){
                        for(int a =0; a < ndim; a++){
                            dFdFFinv[gamma*ndim*ndim*ndim*ndim*nstoch + delta*ndim*ndim*ndim*ndim + a*ndim*ndim*ndim + c *ndim*ndim + a*ndim +d] -= C_Stoch[delta + beta*nstoch + gamma*nstoch*nstoch]*FInv[beta*ndim*ndim + d*ndim + c];
                        }
                    }
                }
            }
        }
    }

    for(int alpha =0; alpha<nstoch; alpha++){
        for(int beta =0; beta < nstoch; beta++){
            for(int gamma =0; gamma < nstoch; gamma++){
                for(int a=0 ; a < ndim; a++){
                    for(int b=0 ; b < ndim; b++){
                        CFinv[gamma*nstoch*ndim*ndim + beta*ndim*ndim + b*ndim + a] += C_Stoch[gamma*nstoch*nstoch + beta*nstoch + alpha]*FInv[alpha*ndim*ndim + b*ndim + a];
                    }
                }
            }
        }
    }

    for(int alpha =0; alpha<nstoch; alpha++){
        for(int beta =0; beta < nstoch; beta++){
            for(int gamma =0; gamma < nstoch; gamma++){
                for(int a=0 ; a < ndim; a++){
                    for(int b=0 ; b < ndim; b++){
                        for(int c=0 ; c < ndim; c++){
                            for(int d=0 ; d < ndim; d++){
                                for(int K=0 ; K < ndim; K++){
                                dFinvTdF[gamma*nstoch*ndim*ndim*ndim*ndim + alpha*ndim*ndim*ndim*ndim + b*ndim*ndim*ndim + a*ndim*ndim + c*ndim + d] += CFinv[gamma*nstoch*ndim*ndim + beta*ndim*ndim + a*ndim + K]*dFdFFinv[beta*ndim*ndim*ndim*ndim*nstoch + alpha*ndim*ndim*ndim*ndim + K*ndim*ndim*ndim + b *ndim*ndim + c*ndim +d];
                                }
                            }
                        }
                    }
                }
            }
        }
    }



    for(int alpha =0; alpha<nstoch; alpha++){
        for(int beta =0; beta < nstoch; beta++){
            for(int gamma =0; gamma < nstoch; gamma++){
                Cmu[gamma*nstoch + beta] += C_Stoch[gamma*nstoch*nstoch + beta*nstoch + alpha]*_mu[alpha];
            }
        }
    }

    for(int alpha =0; alpha<nstoch; alpha++){
        for(int beta =0; beta < nstoch; beta++){
            for(int gamma =0; gamma < nstoch; gamma++){
                for(int c =0; c < ndim; c++){
                    for(int d =0; d < ndim; d++){
                        for(int a =0; a < ndim; a++){
                            for(int b =0; b < ndim; b++){
                             dmuFinvTdF[gamma*nstoch*ndim*ndim*ndim*ndim + alpha*ndim*ndim*ndim*ndim + a*ndim*ndim*ndim + b*ndim*ndim + c*ndim +d] += Cmu[gamma*nstoch + beta]*dFinvTdF[beta*ndim*ndim*ndim*ndim*nstoch + alpha*ndim*ndim*ndim*ndim + a*ndim*ndim*ndim + b *ndim*ndim + c*ndim +d];
                            }
                        }
                    }
                }
            }
        }
    }
    



    DivisionRandomVariableRandomVariableStochastic(ones, Jacobian, C_Stoch, one_overJe);
    derivativeProductRandomConstantRandomVariable(one_overJe,dJdF,dlogJedF,C_Stoch,ndim,nstoch);
    derivativeProductRandomConstantRandomVariable(_lambda,dlogJedF,dlambdalogJedF,C_Stoch,ndim,nstoch);
    derivativeProductRandomVariableTensorRandomVariable(lambda_logJe, FInvT, dlambdalogJedF, dFinvTdF, dlambdalogJeFinvTdF, C_Stoch, ndim, nstoch);
    derivativeProductRandomVariableTensorRandomVariable(theta_power, muF, dtheta_powerdF, dmuFdF, dmuthetapowerdF, C_Stoch, ndim, nstoch);  



}