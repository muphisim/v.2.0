//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file stochastiAreaGrowth.cpp
  \brief This file contains all functions related to isotropic growth constitutive model. There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/
#include "stochasticAreaGrowth.h"
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
classStochasticAreaGrowth::classStochasticAreaGrowth(int ndim, string name, string name_cons, double rho,
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
    this->_n0.resize(ndim*number_stochastic,0);
    _n0[0] = 0;
    _n0[1] = 1;
    if (ndim == 3) {
        _n0[0] = 0;
        _n0[1] = 0;
        _n0[2] = 1;
    }
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
        if (Stochastic_parameters[j] == 4) {
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                for (int k =0; k < ndim; k++){
                this->_n0[ndim*(offset + i) + k] = Stochastic_function[i + offset - 1];
                }
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
    stochasticDyadicProduct(_n0, ndim, _n0, nstoch, C_Stoch, _n0Xn0);

    vector<double> dfdf(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    
    for(int alpha =0; alpha < nstoch; alpha++){
        for(int a =0; a < ndim; a++){
            for(int b =0; b < ndim; b++){
                    dfdf [alpha*ndim*ndim*ndim*ndim*nstoch + alpha*ndim*ndim*ndim*ndim + b*ndim*ndim*ndim + a*ndim*ndim + b*ndim + a] += 1;               
            }
        }
    }
    this->dFdF = dfdf;
    
    vector<double>CN(nstoch*nstoch*ndim,0);

    for(int gamma = 0; gamma <nstoch; gamma++){
        for(int alpha = 0; alpha < nstoch; alpha++){
            for(int j = 0; j < ndim; j++){
                for(int beta = 0; beta < nstoch; beta ++){
                    CN[gamma*nstoch*ndim + alpha*ndim + j] += C[gamma*nstoch*nstoch + beta*nstoch + alpha]*_n0[beta*ndim +j];
                }
            }
        }
    }

    vector<double>dndF(nstoch*nstoch*ndim*ndim*ndim,0);
    for(int gamma =0; gamma<nstoch; gamma++){
        for(int delta =0; delta<nstoch; delta++){
            for(int beta =0; beta<nstoch; beta++){
                for(int i=0; i <ndim; i++){
                    for(int l = 0; l < ndim; l++){
                        dndF[gamma*nstoch*ndim*ndim*ndim + delta*ndim*ndim*ndim + i *ndim*ndim + i*ndim + l] += C[gamma*nstoch*nstoch + beta*nstoch + delta]*_n0[beta*ndim +l];
                    }
                }
            }
        }
    }
    vector<double>dnXn0df(nstoch*nstoch*ndim*ndim*ndim*ndim,0);

    for(int gamma =0; gamma<nstoch; gamma++){
        for(int delta =0; delta<nstoch; delta++){
            for(int alpha =0; alpha<nstoch; alpha++){
                for(int i=0; i <ndim; i++){
                    for(int j = 0; j < ndim; j++){
                        for(int k=0; k <ndim; k++){
                            for(int l = 0; l < ndim; l++){
                                dnXn0df[gamma*nstoch*ndim*ndim*ndim*ndim +delta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l] += dndF[alpha*nstoch*ndim*ndim*ndim + delta*ndim*ndim*ndim + i *ndim*ndim + k*ndim + l]*CN[gamma*nstoch*ndim + alpha*ndim +j];
                                
                            }
                        }
                    }
                }
            }
        }
    }
    this->dnXn0dF = dnXn0df;


};

classStochasticAreaGrowth::~classStochasticAreaGrowth()
{
 
}

void classStochasticAreaGrowth::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
     _term = new MechanicalTermStochastic();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classStochasticAreaGrowth::soundSpeed() const 
{
    double rho = this->getRho();
    double sound = sqrt(_Young[0] / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classStochasticAreaGrowth::initIntVars(vector<double> &intVars) {
    intVars.resize(nstoch,0);
    intVars[0] = 1.0;

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classStochasticAreaGrowth::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
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

    vector<double> sqrt_theta(nstoch, 0);
    vector<double> sqrt_theta_n0Xn0(ndim*ndim*nstoch, 0);
    vector<double> PK1e(nstoch *ndim *ndim, 0);
    vector<double> Fe(nstoch *ndim *ndim, 0);
    vector<double> delta(nstoch*ndim*ndim,0);
    vector<double> Fg(nstoch *ndim *ndim, 0);
    vector<double> Fginv(nstoch *ndim *ndim, 0);
    vector<double> FginvT(nstoch *ndim *ndim, 0);
    // Deformed surface normal#
    nthPowerIntegral(theta, sqrt_theta, C_Stoch, 0.5);
    multiSTensorRandomVariableStochastic(_n0Xn0, ndim, ndim, sqrt_theta, nstoch, C_Stoch, sqrt_theta_n0Xn0);

    for(int i=0;i<ndim;i++){
        delta[i*ndim+i] = 1;
    }

    for (int alpha = 0; alpha < nstoch; alpha++){
        for (int i = 0; i < ndim; i++) {
            for (int j = 0; j < ndim; j++) {
            Fg[alpha*ndim*ndim + i * ndim + j] = sqrt_theta[alpha] * delta[i * ndim + j] + _n0Xn0[alpha*ndim*ndim + i * ndim + j] - sqrt_theta_n0Xn0[alpha*ndim*ndim + i * ndim + j];
            }
        }
    }

    InverseTensorStochastic(Fg, ndim, nstoch, C_Stoch, Fginv);
    transposeStochastic(Fginv, ndim, nstoch, FginvT);

    stress(PK2, PK1, FCurr, sqrt_theta, FginvT);
    vector<double>FIniInvT(nstoch*ndim*ndim,0);
    vector<double>PK2_FIniInvT(nstoch*ndim*ndim,0);
    vector<double>FIniInv_PK2_FIniInvT(nstoch*ndim*ndim,0);
    if(_stochasticMapping){
        transposeStochastic(FIniInv, ndim, nstoch, FIniInvT);
        multTensorTensor3Stochastic(PK2, ndim, ndim, FIniInvT, ndim, PK2_FIniInvT, nstoch, C_Stoch);
        multTensorTensor3Stochastic(FIniInv, ndim, ndim, PK2_FIniInvT, ndim, FIniInv_PK2_FIniInvT, nstoch, C_Stoch);
        multiSTensorRandomVariableStochastic(FIniInv_PK2_FIniInvT, ndim, ndim, det_FIni, nstoch, C_Stoch, PK2);
    }


    vector<double> dPK1dF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> PK1_var(ndim *ndim *nstoch, 0);
    vector<double> PK2_var(ndim *ndim *nstoch, 0);
    vector<double> FCurr_var(ndim *ndim *nstoch, 0);
    double tol = 0.0000001;
    
    double val = 0;
    if (flagTanMod) {
        build_derivatives(FCurr, dPK1dF, FginvT, theta, C_Stoch, tol, nstoch, ndim);
        classTensor6* tanModuli = currState.getTangentStochastic();
        for (int iota = 0; iota < nstoch; iota++) {
            for (int i = 0; i < ndim; i++) {
                for (int J = 0; J < ndim; J++) {

                    for (int eps = 0; eps < nstoch; eps++) {
                        for (int k = 0; k < ndim; k++) {
                            for (int M = 0; M < ndim; M++) {
                                val = dPK1dF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J] ;
                                tanModuli->setValues(eps, iota, k, M, i, J, val);
                            }
                        }
                    }
                }
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
void classStochasticAreaGrowth::computeFe(const vector<double> &FCurr, const vector<double> &sqrt_theta, vector<double> &Fe) {
    int number_stochastic = nstoch;
    int ndim = getDim();
    setAll(Fe,0);
    vector<double>n(ndim*nstoch,0);
    vector<double>nXn0(ndim*ndim*nstoch,0);
    vector<double>ones(nstoch,0);
    vector<double>inv_sqrt_theta(nstoch,0);
    vector<double>inv_sqrt_thetaF(ndim*ndim*nstoch,0);
    vector<double>inv_sqrt_theta_nXn0(ndim*ndim*nstoch,0);
    ones[0] = 1;

    DivisionRandomVariableRandomVariableStochastic(ones, sqrt_theta, C_Stoch, inv_sqrt_theta);

    multTensorVectorStochastic(FCurr, ndim, nstoch, _n0, C_Stoch, n);
    stochasticDyadicProduct(n, ndim, _n0, nstoch, C_Stoch, nXn0);
    multiSTensorRandomVariableStochastic(nXn0, ndim, ndim, inv_sqrt_theta, nstoch, C_Stoch,inv_sqrt_theta_nXn0);

    multiSTensorRandomVariableStochastic(FCurr, ndim, ndim, inv_sqrt_theta, nstoch, C_Stoch,
                                     inv_sqrt_thetaF);

    for (int alpha = 0; alpha < nstoch; alpha++){
        for (int i = 0; i < ndim; i++) {
            for (int j = 0; j < ndim; j++) {
            Fe[alpha*ndim*ndim + i * ndim + j] = inv_sqrt_thetaF[alpha*ndim*ndim + i * ndim + j] + nXn0[alpha*ndim*ndim + i * ndim + j] - inv_sqrt_theta_nXn0[alpha*ndim*ndim + i * ndim + j];
            }
        }
    }

};

void classStochasticAreaGrowth::stress(vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr, const vector<double> &sqrtTheta, const vector<double> &FginvT) {
    int number_stochastic = nstoch;
    int ndim = getDim();
    setAll(PK1,0);
    vector<double>Fe(ndim*ndim*nstoch,0);
    computeFe(FCurr, sqrtTheta, Fe);                                                                                                                                
    vector<double>Je(nstoch,0);
    vector<double>logJe(nstoch,0);
    vector<double>FeInv(nstoch*ndim*ndim,0);
    vector<double>FInv(nstoch*ndim*ndim,0);
    vector<double>FeInvT(nstoch*ndim*ndim,0);
    vector<double>firstPart(nstoch*ndim*ndim,0);
    vector<double>secondPart(nstoch*ndim*ndim,0);
    vector<double>thirdPart(nstoch*ndim*ndim,0);
    vector<double>lambda_logJe(nstoch,0);
    vector<double>PK1e(nstoch*ndim*ndim,0);
    
    DeterminantStochasticTensor(Fe, ndim, nstoch, Je, C_Stoch);
    logIntegral(Je, logJe, C_Stoch);
    InverseTensorStochastic(FCurr, ndim, nstoch, C_Stoch, FInv);
    InverseTensorStochastic(Fe, ndim, nstoch, C_Stoch, FeInv);
    transposeStochastic(FeInv, ndim, nstoch, FeInvT);

    multiSTensorRandomVariableStochastic(Fe, ndim, ndim, _mu, nstoch, C_Stoch,
                                                firstPart);

    multiSTensorRandomVariableStochastic(FeInvT, ndim, ndim, _mu, nstoch, C_Stoch,
                                                thirdPart);                                            
    multiRandomVariableRandomVariableStochastic(_lambda, logJe, C_Stoch,
                                                 lambda_logJe);                                            
    multiSTensorRandomVariableStochastic(FeInvT, ndim, ndim, lambda_logJe, nstoch, C_Stoch,
                                                secondPart);

    for (int i = 0; i < nstoch; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                PK1e[i*ndim*ndim + j*ndim +k] = firstPart[i*ndim*ndim + j*ndim +k] - thirdPart[i*ndim*ndim + j*ndim +k] + secondPart[i*ndim*ndim + j*ndim +k];
            }
        }
    }
    multTensorTensor3Stochastic(PK1e, ndim, ndim, FginvT, ndim, PK1, nstoch, C_Stoch);
    multTensorTensor3Stochastic(FInv, ndim, ndim, PK1, ndim, PK2, nstoch, C_Stoch);
    

}

void classStochasticAreaGrowth::build_derivatives(const vector<double> &FCurr, vector<double> &dPK1dF, const vector<double> &FgInvT, const vector<double> &theta, const vector<double> &C_Stoch, const double tol, const int nstoch, const int ndim){

    setAll(dPK1dF,0);

    vector<double> lambda_logJe(nstoch,0);
    vector<double> Je(nstoch,0);
    vector<double> Je_Perturbed(nstoch,0);
    vector<double> FPerturbed(ndim*ndim*nstoch,0);
    vector<double> FePerturbed(ndim*ndim*nstoch,0);
    vector<double> Fe(ndim *ndim *nstoch, 0);
    vector<double> FeInv(ndim *ndim *nstoch, 0);
    vector<double> FeInvT(ndim *ndim *nstoch, 0);
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
    vector<double> dFdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dFedF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dJedF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dFedF_Feinv(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dFeinvdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dFeinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dmuFeinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dJdF(nstoch*nstoch*ndim*ndim,0);
    vector<double> dFinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dFdFFinv(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dmuFedF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dlambdalogJeFeinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> sqrt_theta(nstoch,0);
    vector<double> inv_sqrt_theta(nstoch,0);
    vector<double> inv_sqrt_thetadFdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> inv_sqrt_thetadndF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dPK1edF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);

    
    vector<double> muFe(nstoch*ndim*ndim,0);

    double power = -0.5;
    nthPowerIntegral(theta, inv_sqrt_theta, C_Stoch, power);
    nthPowerIntegral(theta, sqrt_theta, C_Stoch, 0.5);
    // Compute F-T
    computeFe(FCurr, sqrt_theta, Fe);
    InverseTensorStochastic(Fe, ndim, nstoch, C_Stoch, FeInv);
    transposeStochastic(FeInv, ndim, nstoch, FeInvT);

    DeterminantStochasticTensor(Fe, ndim, nstoch, Je, C_Stoch);
    logIntegral(Je, logJe, C_Stoch);
    multiRandomVariableRandomVariableStochastic(_lambda, logJe, C_Stoch, lambda_logJe);
    multiSTensorRandomVariableStochastic(Fe, ndim, ndim, _mu, nstoch, C_Stoch,
                                                 muFe);


    for(int alpha =0; alpha < nstoch; alpha++){
        for(int b =0; b < ndim; b++){
            for(int a =0; a < ndim; a++){
                dFdF[alpha*ndim*ndim*ndim*ndim*nstoch + alpha*ndim*ndim*ndim*ndim + b*ndim*ndim*ndim + a*ndim*ndim + b*ndim +a] = 1;            
            }
        }
    }
    derivativeProductRandomConstantTensorRandomVariable(inv_sqrt_theta, dFdF, inv_sqrt_thetadFdF, C_Stoch, ndim, nstoch);
    derivativeProductRandomConstantTensorRandomVariable(inv_sqrt_theta, dnXn0dF, inv_sqrt_thetadndF, C_Stoch, ndim, nstoch);

    for(int alpha =0; alpha < nstoch; alpha++){
        for(int beta =0; beta < nstoch; beta++){
            for(int i =0; i < ndim; i++){
                for(int j =0; j < ndim; j++){
                    for(int k =0; k < ndim; k++){
                        for(int l =0; l < ndim; l++){
                            dFedF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l] = inv_sqrt_thetadFdF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l] + dnXn0dF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l] -inv_sqrt_thetadndF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l]; 
                        }
                    }
                }
            }
        }
    }

    for(int beta =0; beta < nstoch; beta++){
        for(int b =0; b < ndim; b++){
            for(int a =0; a < ndim; a++){
                FPerturbed = FCurr;
                FPerturbed[beta*ndim*ndim + b*ndim +a ] += tol;
                computeFe(FPerturbed, sqrt_theta, FePerturbed);
                DeterminantStochasticTensor(FePerturbed, ndim, nstoch, Je_Perturbed, C_Stoch);
                for(int alpha =0; alpha < nstoch; alpha++){
                    dJedF[alpha*ndim*ndim*nstoch + beta*ndim*ndim + b*ndim +a] = (Je_Perturbed[alpha] - Je[alpha])/tol;
                }
            }
        }
    }

    derivativeProductTensorConstantTensorRandomVariable(dFedF, FeInv, "RIGHT", dFedF_Feinv, C_Stoch, ndim, nstoch);
    derivativeProductTensorConstantTensorRandomVariable(FeInv, dFedF_Feinv, "LEFT", dFeinvdF, C_Stoch, ndim, nstoch);

    for(int alpha =0; alpha < nstoch; alpha++){
        for(int beta =0; beta < nstoch; beta++){
            for(int i =0; i < ndim; i++){
                for(int j =0; j < ndim; j++){
                    for(int k =0; k < ndim; k++){
                        for(int l =0; l < ndim; l++){
                            dFeinvTdF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l] = -dFeinvdF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + j*ndim*ndim*ndim + i*ndim*ndim + k*ndim +l]; 
                        }
                    }
                }
            }
        }
    }

    derivativeProductRandomConstantTensorRandomVariable(_mu, dFedF, dmuFedF, C_Stoch, ndim, nstoch);
    derivativeProductRandomConstantTensorRandomVariable(_mu, dFeinvTdF, dmuFeinvTdF, C_Stoch, ndim, nstoch);

    



    DivisionRandomVariableRandomVariableStochastic(ones, Je, C_Stoch, one_overJe);
    derivativeProductRandomConstantRandomVariable(one_overJe,dJedF,dlogJedF,C_Stoch,ndim,nstoch);
    derivativeProductRandomConstantRandomVariable(_lambda,dlogJedF,dlambdalogJedF,C_Stoch,ndim,nstoch);
    derivativeProductRandomVariableTensorRandomVariable(lambda_logJe, FeInvT, dlambdalogJedF, dFeinvTdF, dlambdalogJeFeinvTdF, C_Stoch, ndim, nstoch);

    for(int alpha =0; alpha < nstoch; alpha++){
        for(int beta =0; beta < nstoch; beta++){
            for(int i =0; i < ndim; i++){
                for(int j =0; j < ndim; j++){
                    for(int k =0; k < ndim; k++){
                        for(int l =0; l < ndim; l++){
                            dPK1edF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l] = dmuFedF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l] + dlambdalogJeFeinvTdF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l] - dmuFeinvTdF[alpha*ndim*ndim*ndim*ndim*nstoch + beta*ndim*ndim*ndim*ndim + i*ndim*ndim*ndim + j*ndim*ndim + k*ndim +l]; 
                        }
                    }
                }
            }
        }
    }

    derivativeProductTensorConstantTensorRandomVariable(dPK1edF, FgInvT, "RIGHT", dPK1dF, C_Stoch, ndim, nstoch);    



}