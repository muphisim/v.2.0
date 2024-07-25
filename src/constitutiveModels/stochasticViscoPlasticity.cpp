//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file stochastiViscoPlasticity.cpp
  \brief This file contains all functions related to isotropic growth constitutive model. There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/
#include "stochasticViscoPlasticity.h"
#include "stochasticHyperElasticStVenantKirchhoff.h"
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
classStochasticViscoPlasticity::classStochasticViscoPlasticity(int ndim, string name, string name_cons, double rho,
                                                         double Young, double nu, double eps_0, double s, double n, int stressState, vector<double> Stochastic_numbers,
                                                     vector<double> Stochastic_parameters,
                                                     vector<double> Stochastic_function, string approximation, int order, int resolution, bool stochasticMapping, string distribution) : 
                                                         
            constitutiveModels(ndim, name,  name_cons, rho){

    
    this->HyperElasticStVenantKirchhoff = new  classStochasticHyperElasticStVenantKirchhoff(ndim, name, name_cons, rho,
                                                                                         Young, nu, stressState,
                                                                                         Stochastic_numbers,
                                                                                         Stochastic_parameters,
                                                                                         Stochastic_function,
                                                                                         approximation, order,
                                                                                         resolution, stochasticMapping, distribution);
    int number_func = 0;
    for (int j = 0; j < Stochastic_numbers.size(); j++) {
        number_func = number_func + Stochastic_numbers[j];
    }
    int number_stochastic = 0;

    if(order == 0){
        number_stochastic = 1;
    } else {
        if (approximation == "Haar") {
            number_stochastic = binomialCoefficients(pow(2, order + 1) + number_func - 1, number_func);
        } else {
            if (approximation == "Mixed") {
            number_stochastic = binomialCoefficients((order + 1)*(2*pow(2,resolution)) + number_func - 1, number_func);
            } else {
                number_stochastic = binomialCoefficients(order + number_func, order);
            }
        }
    }
    this->_stochasticMapping = stochasticMapping;
    vector<double> zeros(number_stochastic, 0);
    this->elasticStiffTensor = new classTensor6(ndim, number_stochastic);

    this->_Young = Young;
    this->_eps_0 = zeros;
    this->_s = zeros;
    this->_n = 0;
    this->_eps_0[0] = eps_0;
    this->_s[0] = s;
    this->_n = n;
    int offset = 1;

    for (int j = 0; j < Stochastic_parameters.size(); j++) {
        if (Stochastic_parameters[j] == 1) {
            offset = offset + Stochastic_numbers[j];
        }
        if (Stochastic_parameters[j] == 2) {
            offset = offset + Stochastic_numbers[j];
        }
        if (Stochastic_parameters[j] == 3) {
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->_eps_0[offset + i] = eps_0 * Stochastic_function[i + offset - 1];
            }
            offset = offset + Stochastic_numbers[j];
        }
        if (Stochastic_parameters[j] == 4) {
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->_s[offset + i] = s * Stochastic_function[i + offset - 1];
            }
            offset = offset + Stochastic_numbers[j];
        }
    }

    
    vector<double> C(number_stochastic * number_stochastic * number_stochastic, 0);
    if(order == 0){
        C[0]=1;
    } else {
        if (approximation == "Haar") {
            Build_C_Haar(order, number_func, C);
            //ConvertToHaar(this->_Young, this->_nu, order, Stochastic_parameters, Stochastic_numbers, distribution);
        } else {
            if (approximation == "Mixed"){
                Build_C_Mixed(order, resolution, number_func, C, distribution);
            } else {
    
                build_ThirdOrderTensor_C(order, number_func, C, distribution);
            }
        }
    }
    this->C_Stoch = C;
    this->nstoch = number_stochastic;


};

classStochasticViscoPlasticity::~classStochasticViscoPlasticity()
{
 
}

void classStochasticViscoPlasticity::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
     _term = new MechanicalTermStochastic();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classStochasticViscoPlasticity::soundSpeed() const 
{
    double rho = this->getRho();
    double sound = sqrt(_Young / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classStochasticViscoPlasticity::initIntVars(vector<double> &intVars) {
    int ndim = getDim();
    intVars.resize(ndim*ndim*nstoch,0);
    for(int i=0; i<ndim; i++){
    intVars[i+ndim*i] = 1;
    }

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classStochasticViscoPlasticity::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
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
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    intVarsCurr = intVarsPrev;

    if(_stochasticMapping) {
        if((timeRun>1+0.0001)&&(timeRun<1+dt+0.0001)){
            initIntVars(intVarsCurr);
            INFO("Reinitilisation");
        }
        InverseTensorStochastic(FIni, ndim, nstoch, C_Stoch, FIniInv);
        multTensorTensor3Stochastic(FCurr, ndim, ndim, FIniInv, ndim, FCurrMecha, nstoch, C_Stoch);
        DeterminantStochasticTensor(FIni, ndim, nstoch, det_FIni, C_Stoch);
    }else {
            FCurrMecha = FCurr;
    }

    computeGLstoch(FCurrMecha, ndim, nstoch, C_Stoch, E);
    //computeEquivalentStrainStochastic(E, ndim, nstoch, C_Stoch, equivalentStrain);
    //computeVolumetricStrainStochastic(E, ndim, nstoch, volumetricStrain);



    computePlasticIncrement(intVarsCurr, PK2, PK1, FCurrMecha, dt);


    //stress(PK2, PK1, FCurr, sqrt_theta, FginvT);
    vector<double>FIniInvT(nstoch*ndim*ndim,0);
    vector<double>PK2_FIniInvT(nstoch*ndim*ndim,0);
    vector<double>FIniInv_PK2_FIniInvT(nstoch*ndim*ndim,0);
    if(_stochasticMapping){
        transposeStochastic(FIniInv, ndim, nstoch, FIniInvT);
        multTensorTensor3Stochastic(PK2, ndim, ndim, FIniInvT, ndim, PK2_FIniInvT, nstoch, C_Stoch);
        multTensorTensor3Stochastic(FIniInv, ndim, ndim, PK2_FIniInvT, ndim, FIniInv_PK2_FIniInvT, nstoch, C_Stoch);
        multiSTensorRandomVariableStochastic(FIniInv_PK2_FIniInvT, ndim, ndim, det_FIni, nstoch, C_Stoch, PK2);
    }
    setAll(PK1,0);
    multTensorTensor3Stochastic(FCurr, ndim, ndim, PK2, ndim, PK1, nstoch, C_Stoch);
    computeCauchyStoch(FCurr, ndim, nstoch, C_Stoch, PK1, Cauchy);
    //computeStochasticVMS(Cauchy, ndim, nstoch, VMS, C_Stoch);
    double val = 0;
    if (flagTanMod) {
        classTensor6* tanModuli = currState.getTangentStochastic();
        static_cast<classStochasticHyperElasticStVenantKirchhoff *>(this->HyperElasticStVenantKirchhoff)->getTanModuli(tanModuli, FCurr, PK2, flagTanMod);
    }

}


/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
void classStochasticViscoPlasticity::computePlasticIncrement(vector<double> &FPCurr, vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr, const double dt) 
{
    double tol_1 = 0.0001;
    double tol_2 = 0.0001;
    double res = 1000;
    int ite = 0;
    int ndim = getDim();
    const vector<double> FPPrev = FPCurr;
    vector<double> S(nstoch*ndim*ndim,0);
    vector<double> SDev(nstoch*ndim*ndim,0);
    vector<double> SigEv(nstoch,0);
    vector<double> E(nstoch*ndim*ndim,0);
    vector<double> C(nstoch*ndim*ndim,0);
    vector<double> Lp(nstoch*ndim*ndim,0);
    vector<double> R(nstoch*ndim*ndim,0);
    vector<double> Eps_P_R(nstoch*ndim*ndim,0);
    vector<double> EpsP(nstoch,0);
    vector<double> dResdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> delta(nstoch*ndim*ndim,0);
    vector<double> Tang(nstoch*ndim*ndim,0);
    vector<double> ResVec(nstoch*ndim*ndim,0);
    vector<double> ResVecPerturbed(nstoch*ndim*ndim,0);
    double Res = 10;
    double ResPerturbed = 0;
    vector<double> FPCurrPerturbed(nstoch*ndim*ndim,0);

    vector<double> FeCurr(nstoch*ndim*ndim,0);
    vector<double> FPInvCurr(nstoch*ndim*ndim,0);
    vector<double> FPInvCurrT(nstoch*ndim*ndim,0);
    vector<double> DetFPCurr(nstoch,0);

    vector<double> PK2_FPInvCurrT(nstoch*ndim*ndim,0);
    vector<double> FPInvCurr_PK2_FPInvCurrT(nstoch*ndim*ndim,0);

    InverseTensorStochastic(FPCurr, ndim, nstoch, C_Stoch, FPInvCurr);
    multTensorTensor3Stochastic(FCurr, ndim, ndim, FPInvCurr, ndim, FeCurr, nstoch, C_Stoch);

    while ((Res>tol_2) || (ite>100)){

        Res = 0;

        computeGLstoch(FeCurr, ndim, nstoch, C_Stoch, E);

        static_cast<classStochasticHyperElasticStVenantKirchhoff *>(this->HyperElasticStVenantKirchhoff)->stress(S,E);
        computeSDev(SDev, C, FeCurr, S);     
        computeSigEq(SDev, C, SigEv);
        computeR(SDev, C, SigEv, R);
        computeEpsP(EpsP, SigEv);
 
        multiSTensorRandomVariableStochastic(R, ndim, ndim, EpsP, nstoch, C_Stoch,
                                                 Eps_P_R);
                                        
        computeRes(FPCurr, FPInvCurr, FPPrev, Eps_P_R, ResVec, Res, dt);                              
        if(Res>tol_2){
            FPCurrPerturbed.resize(nstoch*ndim*ndim,0);
            dResdF.resize(nstoch*nstoch*ndim*ndim*ndim*ndim,0);

            for(int alpha =0; alpha < nstoch; alpha++){
                for(int i =0; i < ndim; i++){
                    for(int j =0; j < ndim; j++){
                    ResPerturbed = 0;
                    FPCurrPerturbed = FPCurr;
                    FPCurrPerturbed[alpha*ndim*ndim + i*ndim +j ] += tol_1;
                    R.resize(nstoch*ndim*ndim,0);
                    C.resize(nstoch*ndim*ndim,0);
                    S.resize(nstoch*ndim*ndim,0);
                    E.resize(nstoch*ndim*ndim,0);
                    SDev.resize(nstoch*ndim*ndim,0);
                    FeCurr.resize(nstoch*ndim*ndim,0);
                    SigEv.resize(nstoch,0);
                    EpsP.resize(nstoch,0);
                    Eps_P_R.resize(nstoch*ndim*ndim,0);
                    ResVecPerturbed.resize(nstoch*ndim*ndim,0);
                    FPInvCurr.resize(ndim*ndim*nstoch,0);
                    FeCurr.resize(ndim*ndim*nstoch,0);
                    InverseTensorStochastic(FPCurrPerturbed, ndim, nstoch, C_Stoch, FPInvCurr);
                    multTensorTensor3Stochastic(FCurr, ndim, ndim, FPInvCurr, ndim, FeCurr, nstoch, C_Stoch);
                    computeGLstoch(FeCurr, ndim, nstoch, C_Stoch, E);
                    static_cast<classStochasticHyperElasticStVenantKirchhoff *>(this->HyperElasticStVenantKirchhoff)->stress(S,E);
                    computeSDev(SDev, C, FeCurr, S);
                    computeSigEq(SDev, C, SigEv);
                    computeR(SDev, C, SigEv, R);
                    computeEpsP(EpsP, SigEv);
                    
                    multiSTensorRandomVariableStochastic(R, ndim, ndim, EpsP, nstoch, C_Stoch,
                                                            Eps_P_R);
                    computeRes(FPCurrPerturbed, FPInvCurr, FPPrev, Eps_P_R, ResVecPerturbed, ResPerturbed, dt); 
                        for(int beta=0; beta < nstoch; beta ++){
                            for(int k =0; k < ndim; k++){
                                for(int l =0; l < ndim; l++){
                                dResdF[beta*ndim*ndim*ndim*ndim*nstoch + alpha*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + l*ndim*ndim + i*ndim +j] = (ResVecPerturbed[beta*ndim*ndim + k*ndim + l]-ResVec[beta*ndim*ndim + k*ndim + l])/tol_1;
                                }
                            }
                        }                
                    }
                }
            }            
            delta.resize(ndim*ndim*nstoch, 0);

            InverseSystem(dResdF, ResVec, ndim, nstoch, delta);

            for(int alpha =0; alpha < nstoch; alpha++){
                for(int i =0; i < ndim; i++){
                    for(int j =0; j < ndim; j++){
                        FPCurr[alpha*ndim*ndim + i*ndim +j] -= delta[alpha*ndim*ndim + i*ndim +j];
                    }
                }
            }
            FPInvCurr.resize(ndim*ndim*nstoch,0);
            FeCurr.resize(ndim*ndim*nstoch,0);
            InverseTensorStochastic(FPCurr, ndim, nstoch, C_Stoch, FPInvCurr);
            multTensorTensor3Stochastic(FCurr, ndim, ndim, FPInvCurr, ndim, FeCurr, nstoch, C_Stoch);
        }
        ite++;
        if(ite>100){
            INFO("RADIAL RETURN DOES NOT CONVERGE");
        }        
    }
    transposeStochastic(FPInvCurr, ndim, nstoch, FPInvCurrT);
    DeterminantStochasticTensor(FPCurr, ndim, nstoch, DetFPCurr, C_Stoch);
    multTensorTensor3Stochastic(S, ndim, ndim, FPInvCurrT, ndim, PK2_FPInvCurrT, nstoch, C_Stoch);
    multTensorTensor3Stochastic(FPInvCurr, ndim, ndim, PK2_FPInvCurrT, ndim, FPInvCurr_PK2_FPInvCurrT, nstoch, C_Stoch);
    multiSTensorRandomVariableStochastic(PK2_FPInvCurrT, ndim, ndim, DetFPCurr, nstoch, C_Stoch,
                                                            PK2);
    multTensorTensor3Stochastic(FCurr, ndim, ndim, PK2, ndim, PK1, nstoch, C_Stoch);                                                        

}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
void classStochasticViscoPlasticity::computeSDev(vector<double> &SDev, vector<double> &C, const vector<double> &FeCurr, const vector<double> &S) 
{
    double tol_1 = 0.1;
    double tol_2 = 0.001;
    double res = 1000;
    int ite = 0;
    int ndim = getDim();
    setAll(SDev,0);
    setAll(C,0);
    vector<double> doubleContractSC(nstoch,0);
    vector<double> SC(nstoch*ndim*ndim,0);
    vector<double> Cinv(nstoch*ndim*ndim,0);
    vector<double> doubleContractSC_times_Cinv(ndim*ndim*nstoch,0);

    multSTensor3FirstTransposeStochastic(FeCurr, ndim, ndim, FeCurr, ndim, C, nstoch, C_Stoch);
    InverseTensorStochastic(C, ndim, nstoch, C_Stoch, Cinv);
    multSTensor3FirstTransposeStochastic(S, ndim, ndim, C, ndim, SC, nstoch, C_Stoch);
    getStochasticTrace(SC,ndim,nstoch,doubleContractSC);
    multiSTensorRandomVariableStochastic(Cinv, ndim, ndim, doubleContractSC, nstoch, C_Stoch,
                                                 doubleContractSC_times_Cinv);

    for(int alpha =0; alpha < nstoch; alpha++){
        for(int i =0; i < ndim; i++){
            for(int j =0; j < ndim; j++){
                SDev[alpha*ndim*ndim + i*ndim + j] = S[alpha*ndim*ndim + i*ndim + j] - (1.0/double(ndim))*doubleContractSC_times_Cinv[alpha*ndim*ndim + i*ndim + j];
            }
        }
    }                                                          
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
void classStochasticViscoPlasticity::computeSigEq(const vector<double> &SDev, const vector<double> &C, vector<double>&SigEv)
{
    int ndim = getDim();
    
    vector<double> CSdev(nstoch*ndim*ndim,0);
    vector<double> SDevC(nstoch*ndim*ndim,0);
    vector<double> CSdev_SdevC(nstoch*ndim*ndim,0);
    vector<double> Trace_CSdev_SdevC(nstoch,0);
    double res = 0;
    setAll(SigEv,0);
    multTensorTensor3Stochastic(C, ndim, ndim, SDev, ndim, CSdev, nstoch, C_Stoch);
    multTensorTensor3Stochastic(SDev, ndim, ndim, C, ndim, SDevC, nstoch, C_Stoch);
    multSTensor3FirstTransposeStochastic(SDevC, ndim, ndim, CSdev, ndim, CSdev_SdevC, nstoch, C_Stoch);
    getStochasticTrace(CSdev_SdevC, ndim, nstoch, Trace_CSdev_SdevC);
    MultiplicationRandomVariableScalar(Trace_CSdev_SdevC, 1.5);
    

    for(int alpha=0; alpha<nstoch;alpha++){
        res += Trace_CSdev_SdevC[alpha];
    }
    if(res<1){
        //INFO("SigEv too small so we suppose it is equal to 0");
        }else{
        double normalisation = 1.0/(Trace_CSdev_SdevC[0]);
        double inv_normalisation = pow(Trace_CSdev_SdevC[0],0.5);
        MultiplicationRandomVariableScalar(Trace_CSdev_SdevC,normalisation);
        //squareRootIntegral(Trace_CSdev_SdevC, SigEv, C_Stoch);
        nthPowerIntegral(Trace_CSdev_SdevC, SigEv, C_Stoch, 0.5);
        MultiplicationRandomVariableScalar(SigEv,inv_normalisation);
    }
                                                
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
void classStochasticViscoPlasticity::computeR(const vector<double> &SDev, const vector<double> &C, const vector<double>&SigEv, vector<double>&R)
{
    int ndim = getDim();
    setAll(R,0);
    vector<double> CSdev(nstoch*ndim*ndim,0);
    vector<double> CSdevC(nstoch*ndim*ndim,0);
    vector<double> InvSigEv(nstoch,0);
    vector<double> ones(nstoch,0);
    ones[0] = 1;
    double res = 0;
    for(int alpha=0; alpha<nstoch;alpha++){
        res += SigEv[alpha];
    }
    if(res<1){
        //INFO("SigEv too small so we suppose R is equal to 0");
        }else{
        DivisionRandomVariableRandomVariableStochastic(ones, SigEv, C_Stoch, InvSigEv);
        MultiplicationRandomVariableScalar(InvSigEv, 1.5);
        multTensorTensor3Stochastic(C, ndim, ndim, SDev, ndim, CSdev, nstoch, C_Stoch);
        multTensorTensor3Stochastic(CSdev, ndim, ndim, C, ndim, CSdevC, nstoch, C_Stoch);
        multiSTensorRandomVariableStochastic(CSdevC, ndim, ndim, InvSigEv, nstoch, C_Stoch,
                                                    R);
    }
      
                                             
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
void classStochasticViscoPlasticity::computeEpsP(vector<double> &EpsP, const vector<double> &SigEv)
{
    int ndim = getDim();
    vector<double> SigEv_over_s(nstoch,0);
    vector<double> SigEv_over_s_powern(nstoch,0);
    vector<double> InvSigEv(nstoch,0);
    vector<double> ones(nstoch,0);
    ones[0] = 1;
    double res = 0;
    setAll(EpsP,0);
    for(int alpha=0; alpha<nstoch;alpha++){
        res += SigEv[alpha];
    }
    if(res<1){
        //INFO("SigEv too small so we suppose EPS_P is equal to 0");
        }else{
        //printVector(SigEv,"SigEv");
        DivisionRandomVariableRandomVariableStochastic(SigEv, _s, C_Stoch, SigEv_over_s);
        double normalize = double(1.0)/double(SigEv_over_s[0]);
        double inv_normalize = pow(SigEv_over_s[0],_n);
        MultiplicationRandomVariableScalar(SigEv_over_s , normalize);
        nthPowerIntegral(SigEv_over_s, SigEv_over_s_powern, C_Stoch, _n);
        MultiplicationRandomVariableScalar(SigEv_over_s_powern , inv_normalize);
        multiRandomVariableRandomVariableStochastic(_eps_0, SigEv_over_s_powern, C_Stoch, EpsP);
    }
                                              
                                             
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/

void classStochasticViscoPlasticity::computeRes(const vector<double> &FPCurr, const vector<double> &FPInvCurr, const vector<double> &FPPrev, const vector<double> &Eps_P_R, vector<double> &ResVec, double &Res, double dt)
{
    int ndim = getDim();
    setAll(ResVec,0);
    Res = 0;
    vector<double> FPdot(nstoch*ndim*ndim,0);
    vector<double> FPdot_FPinv(nstoch*ndim*ndim,0);

    for(int alpha =0; alpha < nstoch; alpha++){
        for(int i =0; i < ndim; i++){
            for(int j =0; j < ndim; j++){
                FPdot[alpha*ndim*ndim + i*ndim + j] = (FPCurr[alpha*ndim*ndim + i*ndim + j] - FPPrev[alpha*ndim*ndim + i*ndim + j])/dt;
            }
        }
    }
    multTensorTensor3Stochastic(FPdot, ndim, ndim, FPInvCurr, ndim, FPdot_FPinv, nstoch, C_Stoch);
    


    for(int alpha =0; alpha < nstoch; alpha++){
        for(int i =0; i < ndim; i++){
            for(int j =0; j < ndim; j++){
                ResVec[alpha*ndim*ndim + i*ndim + j] = FPdot_FPinv[alpha*ndim*ndim + i*ndim + j] - Eps_P_R[alpha*ndim*ndim + i*ndim + j];
            }
        }
    }

    for(int alpha =0; alpha < nstoch; alpha++){
        for(int i =0; i < ndim; i++){
            for(int j =0; j < ndim; j++){
                Res += ResVec[alpha*ndim*ndim + i*ndim + j] * ResVec[alpha*ndim*ndim + i*ndim + j];
            }
        }
    }



                                             
}
