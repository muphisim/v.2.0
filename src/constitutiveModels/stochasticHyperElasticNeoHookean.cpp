//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file stochasitcHyperElasticNeoHookean.cpp
  \brief 
*/
#include "stochasticHyperElasticNeoHookean.h"
#include "maths.h"
#include "classGPs.h"


/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young Young modulus
  @param[in] nu Poisson ratio
*/
classStochasticHyperElasticNeoHookean::classStochasticHyperElasticNeoHookean(int ndim, string name,
                                                                                           string name_cons, double rho,
                                                                                           double Young, double nu,
                                                                                           vector<double> Stochastic_numbers,
                                                                                           vector<double> Stochastic_parameters,
                                                                                           vector<double> Stochastic_function,
                                                                                           string approximation,
                                                                                           int order, int resolution, bool stochasticMapping, string distribution)
        : constitutiveModels(ndim, name, name_cons, rho) {
    
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
    this->Young = zeros;
    this->nu = zeros;
    this->Young[0] = Young;
    this->nu[0] = nu;
    int offset = 1;
    for (int j = 0; j < Stochastic_parameters.size(); j++) {
        if (Stochastic_parameters[j] == 1) {
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->Young[offset + i] = Young * Stochastic_function[i];
                //count++;
            }
        }
        if (Stochastic_parameters[j] == 2) {
            if(Stochastic_parameters.size()==2){
                offset = offset + Stochastic_numbers[0];
            }
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->nu[offset + i] = nu * Stochastic_function[i + offset - 1];
                offset++;
            }
        }
    }
    
    vector<double> C(number_stochastic * number_stochastic * number_stochastic, 0);
    if (approximation == "Haar") {
        Build_C_Haar(order, number_func, C);
        ConvertToHaar(this->Young, this->nu, order, Stochastic_parameters, Stochastic_numbers, distribution);
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
    vector<double> K1(number_stochastic, 0);
    vector<double> twicemu(number_stochastic, 0);
    vector<double> oneminusnu_carre(number_stochastic, 0);
    vector<double> Young_times_nu(number_stochastic, 0);
    vector<double> fac(number_stochastic, 0);
    vector<double> Young_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> nutimesYoung_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> oneplusnu_times_oneminustwonu(number_stochastic, 0);
    vector<double> oneminustwicenu(number_stochastic, 0);
    vector<double> factor(number_stochastic, 0);
    vector<double> factor_times_nu(number_stochastic, 0);
    vector<double> young_over_oneminustwonu(number_stochastic, 0);
    oneminustwicenu = this->nu;
    vector<double> oneplusnu = this->nu;
    vector<double> _nu = this->nu;
    vector<double> _Young = this->Young; 
    linearCombinationVariableStochastic(oneplusnu, 1, 1);
    linearCombinationVariableStochastic(oneminustwicenu, -2, 1);


    DivisionRandomVariableRandomVariableStochastic(_Young, oneplusnu, C_Stoch, mu);
    vector<double> Young_over_oneplusnu = mu;
    MultiplicationRandomVariableScalar(mu, 0.5);

   
    DivisionRandomVariableRandomVariableStochastic(Young_over_oneplusnu, oneminustwicenu, C_Stoch, factor);
    multiRandomVariableRandomVariableStochastic(factor, _nu, C_Stoch, factor_times_nu);
    lambda = factor_times_nu;
    DivisionRandomVariableRandomVariableStochastic(_Young, oneminustwicenu, C_Stoch, young_over_oneminustwonu);
    K1 = young_over_oneminustwonu;
    MultiplicationRandomVariableScalar(K1, double(1./3));
    this->lambda = lambda;
    this->mu = mu;
    this->K1 = K1;
    
    // build dFdF once and for all;
    vector<double> dfdf(nstoch*nstoch*ndim*ndim*ndim*ndim,0);

    for(int alpha =0; alpha < nstoch; alpha++){
        for(int a =0; a < ndim; a++){
            for(int b =0; b < ndim; b++){
                        dfdf [alpha*ndim*ndim*ndim*ndim*nstoch + alpha*ndim*ndim*ndim*ndim + b*ndim*ndim*ndim + a*ndim*ndim + b*ndim + a] = 1;
            }
        }
    }
    this->dFdF = dfdf;


}

classStochasticHyperElasticNeoHookean::~classStochasticHyperElasticNeoHookean()
{
}

void classStochasticHyperElasticNeoHookean::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
     _term = new MechanicalTermStochastic();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classStochasticHyperElasticNeoHookean::soundSpeed() const{
    double rho = this->getRho();
    double Young = this->Young[0];
    //  double factornu = (1.-nu)/((1.+nu)*(1.-2.*nu));
    double sound = sqrt(Young / rho);
    return sound;
}

/*! \brief Compute derivatives.
  @param[in&out] dJdF
  @param[in&out] dFdF
  @param[in&out] dFdinvT
  @param[in] FCurr
  @param[in] tol
  @param[in] nstoch
 */
void classStochasticHyperElasticNeoHookean::build_derivatives_NeoHookean(const vector<double> &FCurr, vector<double> &dJdF, vector<double> &dFinvTdF, vector<double> &dtraceBdF, vector<double> &C_Stoch, const double tol, const int nstoch, const int ndim) {
    setAll(dtraceBdF,0);
    vector<double> Jacobian(nstoch,0);
    vector<double> Jacobian_Perturbed(nstoch,0);
    vector<double> FPerturbed(ndim*ndim*nstoch,0);
    vector<double> FInv(ndim *ndim *nstoch, 0);
    vector<double> FInvT(ndim *ndim *nstoch, 0);
    vector<double> FInvPerturbed(ndim *ndim *nstoch, 0);
    vector<double> FInvTPerturbed(ndim *ndim *nstoch, 0);
    InverseTensorStochastic(FCurr, ndim, nstoch, C_Stoch, FInv);
    transposeStochastic(FInv, ndim, nstoch, FInvT);
    FPerturbed = FCurr;
    DeterminantStochasticTensor(FCurr, ndim, nstoch, Jacobian, C_Stoch);
    vector<double> dFdFFinv(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> CFinv(nstoch*nstoch*ndim*ndim,0);

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
        for(int eps =0; eps < nstoch; eps++){
                for(int a =0; a < ndim; a++){
                    for(int b =0; b < ndim; b++){
                        for(int alpha =0; alpha < nstoch; alpha++){
                        dtraceBdF[gamma * ndim * ndim * nstoch + eps *ndim *ndim + a *ndim + b] += 2*C_Stoch[gamma*nstoch*nstoch + eps*nstoch + alpha]*FCurr[alpha*ndim*ndim + a *ndim + b];
                    }
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


}



/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classStochasticHyperElasticNeoHookean::initIntVars(vector<double> &intVars) {
    // This constitutive model does not have any internal variable. The arrays is not modified
    // The last position is always for the number of remeshings. MuPhiSim Internal variable
    intVars.push_back(0.0); // This internal variable is initialised to 0
}


void classStochasticHyperElasticNeoHookean::stress(vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr) {
    int number_stochastic = nstoch;
    int ndim = getDim();
    setAll(PK1,0);
    vector<double> FCurrT(ndim *ndim *number_stochastic,0);
    vector<double> FInv(ndim *ndim *number_stochastic, 0);
    vector<double> FInvT(ndim *ndim *number_stochastic, 0);
    vector<double> Jacobian(number_stochastic,0);
    vector<double> Jacobian_minusone(number_stochastic,0);
    vector<double> JJ_minusone(number_stochastic,0);
    vector<double> K1_JJ_minusone(number_stochastic,0);
    vector<double> traceB(number_stochastic,0);
    vector<double> B(ndim *ndim *number_stochastic,0);
    vector<double> factor1(ndim *ndim *number_stochastic,0);
    vector<double> factor2(ndim *ndim *number_stochastic,0);
    vector<double> factor3(ndim *ndim *number_stochastic,0);
    vector<double> firstPart(ndim *ndim *number_stochastic,0);
    vector<double> Jacobian_power(number_stochastic,0);
    vector<double> mu_over_jacobian_power(number_stochastic,0);
    
     

    DeterminantStochasticTensor(FCurr, ndim, number_stochastic, Jacobian, C_Stoch);
    Jacobian_minusone = Jacobian;
    Jacobian_minusone[0] = Jacobian_minusone[0] - 1.0;
    
    InverseTensorStochastic(FCurr, ndim, number_stochastic, C_Stoch, FInv);
    transposeStochastic(FInv, ndim, number_stochastic, FInvT);
    transposeStochastic(FCurr, ndim, number_stochastic, FCurrT);
    
    multSTensor3SecondTransposeStochastic(FCurr, ndim, ndim, FCurr, ndim, B, number_stochastic, C_Stoch);

    for (int m = 0; m < number_stochastic; m++) {
        for (int i = 0; i < ndim; i++) {
            traceB[m] += B[m * ndim * ndim + i * ndim + i];
        }
    }

    multiRandomVariableRandomVariableStochastic(Jacobian, Jacobian_minusone, C_Stoch,
                                                 JJ_minusone);
    multiRandomVariableRandomVariableStochastic(K1, JJ_minusone, C_Stoch,
                                                 K1_JJ_minusone);                                            
    multiSTensorRandomVariableStochastic(FInvT, ndim, ndim, K1_JJ_minusone, number_stochastic, C_Stoch,
                                                 factor3);
    multiSTensorRandomVariableStochastic(FInvT, ndim, ndim, traceB, number_stochastic, C_Stoch,
                                                 factor2);
    double power = ((5.0 - double(ndim))/ double(ndim));                                                                                                                                  
    nthPowerIntegral(Jacobian, Jacobian_power, C_Stoch, power);
    for (int i = 0; i < number_stochastic; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                factor1[i*ndim*ndim + j*ndim +k] = FCurr[i*ndim*ndim + j*ndim +k] -1*(1/(double(ndim)))*factor2[i*ndim*ndim + j*ndim +k];
            }
        }
    }
    
    DivisionRandomVariableRandomVariableStochastic(mu, Jacobian_power, C_Stoch,
                                                   mu_over_jacobian_power);
                                                  
    multiSTensorRandomVariableStochastic(factor1, ndim, ndim, mu_over_jacobian_power, number_stochastic, C_Stoch,
                                                firstPart);
    
    
    for (int i = 0; i < number_stochastic; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                PK1[i*ndim*ndim + j*ndim +k] = firstPart[i*ndim*ndim + j*ndim +k] + factor3[i*ndim*ndim + j*ndim +k];
            }
        }
    }

    multTensorTensor3Stochastic(FInv, ndim, ndim, PK1, ndim, PK2, nstoch, C_Stoch);
    

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classStochasticHyperElasticNeoHookean::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    //if(timeRun > 1){
    //    this->Young[0] = 1e7;
    //}
    
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    classStates& initState = GaussPoint->getInitialState();
    vector<double> &FIni = initState.getDeformationGradient();
    vector<double> &FCurr = currState.getDeformationGradient();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    int ndim = getDim();
    vector<double>FIniInv(nstoch*ndim*ndim,0);
    vector<double>FCurrMecha(nstoch*ndim*ndim,0);
    vector<double>det_FIni(nstoch,0);
    double power = -((5.0 - double(ndim))/ double(ndim)); 
    double power_minusone = power - 1.0 ;
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
    vector<double> PK2(ndim *ndim *nstoch, 0);
    vector<double> FInv(ndim *ndim *nstoch, 0);
    vector<double> FInvT(ndim *ndim *nstoch, 0);

    InverseTensorStochastic(FCurrMecha, ndim, nstoch, C_Stoch, FInv);
    transposeStochastic(FInv, ndim, nstoch, FInvT);

    stress(PK2, PK1, FCurrMecha);
    vector<double> dJdF(ndim*ndim*nstoch*nstoch,0);
    vector<double> dJpowerdF(ndim*ndim*nstoch*nstoch,0);
    vector<double> dTraceBdF(ndim*ndim*nstoch*nstoch,0);
    vector<double> dFInvTdF(ndim*ndim*ndim*ndim*nstoch*nstoch,0);
    vector<double> Jacobian(nstoch,0);
    vector<double> Jacobian_power(nstoch,0);
    vector<double> Jacobian_power_minus_one(nstoch,0);
    vector<double> mu_Jpower(nstoch,0);
    vector<double> dK1JdF(ndim*ndim*nstoch*nstoch,0);
    vector<double> dJ_squaredF(ndim*ndim*nstoch*nstoch,0);
    vector<double> dK1J_squaredF(ndim*ndim*nstoch*nstoch,0);
    vector<double> K1J(nstoch,0);
    vector<double> K1J_square(nstoch,0);
    vector<double> mu_jpower_traceB(nstoch,0);
    vector<double> dmuJpowerdF(ndim*ndim*nstoch*nstoch,0);
    vector<double> dfactor1dF(ndim*ndim*ndim*ndim*nstoch*nstoch,0);
    vector<double> dfactor2dF(ndim*ndim*ndim*ndim*nstoch*nstoch,0);
    vector<double> dfactor3dF(ndim*ndim*ndim*ndim*nstoch*nstoch,0);
    vector<double> dfactor4dF(ndim*ndim*ndim*ndim*nstoch*nstoch,0);


    DeterminantStochasticTensor(FCurrMecha, ndim, nstoch, Jacobian, C_Stoch);
    nthPowerIntegral(Jacobian, Jacobian_power, C_Stoch, power);
    nthPowerIntegral(Jacobian, Jacobian_power_minus_one, C_Stoch, power_minusone);
    double tol = 0.000001;
    build_derivatives_NeoHookean(FCurrMecha, dJdF, dFInvTdF, dTraceBdF, C_Stoch, tol, nstoch, ndim);
    derivativeProductRandomConstantRandomVariable(Jacobian,dJdF,dJ_squaredF,C_Stoch,ndim,nstoch);
    // PK1 = mu * j-power * (F - 1/ndim trace(B)*F^-t) + K1J^2*F-t - K1J*F-t
    // PK1 = factor_1 + factor_2 + factor_3 + factor_4

    // building dFactor1dF
    derivativeProductRandomConstantRandomVariable(Jacobian_power_minus_one,dJdF,dJpowerdF,C_Stoch,ndim,nstoch);
    derivativeProductRandomConstantRandomVariable(mu,dJpowerdF,dmuJpowerdF,C_Stoch,ndim,nstoch);
    MultiplicationRandomVariableScalar(dmuJpowerdF,power);
    multiRandomVariableRandomVariableStochastic(mu, Jacobian_power, C_Stoch,
                                                 mu_Jpower);
    derivativeProductRandomVariableTensorRandomVariable(mu_Jpower, FCurr, dmuJpowerdF,dFdF, dfactor1dF, C_Stoch, ndim, nstoch); 

    //building dFactor4dF
    multiRandomVariableRandomVariableStochastic(K1, Jacobian, C_Stoch,
                                                 K1J);
    derivativeProductRandomConstantRandomVariable(K1,dJdF,dK1JdF,C_Stoch,ndim,nstoch);
    derivativeProductRandomVariableTensorRandomVariable(K1J, FInvT, dK1JdF,dFInvTdF, dfactor4dF, C_Stoch, ndim, nstoch); 

    // building dFactor3dF                                                                                             
                                             
    multiRandomVariableRandomVariableStochastic(K1J, Jacobian, C_Stoch,
                                                 K1J_square);
    MultiplicationRandomVariableScalar(dJ_squaredF,2.0);
    derivativeProductRandomConstantRandomVariable(K1,dJ_squaredF,dK1J_squaredF,C_Stoch,ndim,nstoch);
    derivativeProductRandomVariableTensorRandomVariable(K1J_square, FInvT, dK1J_squaredF,dFInvTdF, dfactor3dF, C_Stoch, ndim, nstoch); 

    // building dFactor2dF

    vector<double>traceB(nstoch,0);
    vector<double>derivative(ndim*ndim*nstoch*nstoch,0);
    vector<double>derivative_tot(ndim*ndim*nstoch*nstoch,0);
    vector<double>dmu_jpower_traceBdF(ndim*ndim*nstoch*nstoch,0);
    vector<double>B(nstoch*ndim*ndim,0);
    vector<double>dmu_power(nstoch*nstoch*ndim*ndim,0);
    vector<double>mu_jpower_trB(nstoch,0);

    multSTensor3SecondTransposeStochastic(FCurrMecha, ndim, ndim, FCurrMecha, ndim, B, nstoch, C_Stoch);

    for (int m = 0; m < nstoch; m++) {
        for (int i = 0; i < ndim; i++) {
        traceB[m] += B[m * ndim * ndim + i * ndim + i];
        }
    }

    multiRandomVariableRandomVariableStochastic(mu_Jpower, traceB, C_Stoch,
                                                 mu_jpower_trB);



    multiRandomVariableRandomVariableStochastic(mu_Jpower, traceB, C_Stoch,
                                                 mu_jpower_traceB);
    derivativeProductRandomVariableRandomVariable(mu_Jpower, traceB, dmuJpowerdF,dTraceBdF, dmu_jpower_traceBdF, C_Stoch, ndim, nstoch);

    derivativeProductRandomVariableTensorRandomVariable(mu_jpower_traceB, FInvT, dmu_jpower_traceBdF,dFInvTdF, dfactor2dF, C_Stoch, ndim, nstoch);  



    vector<double>FIniInvT(nstoch*ndim*ndim,0);
    vector<double>PK2_FIniInvT(nstoch*ndim*ndim,0);
    vector<double>FIniInv_PK2_FIniInvT(nstoch*ndim*ndim,0);

    if(_stochasticMapping){
        transposeStochastic(FIniInv, ndim, nstoch, FIniInvT);
        multTensorTensor3Stochastic(PK2, ndim, ndim, FIniInvT, ndim, PK2_FIniInvT, nstoch, C_Stoch);
        multTensorTensor3Stochastic(FIniInv, ndim, ndim, PK2_FIniInvT, ndim, FIniInv_PK2_FIniInvT, nstoch, C_Stoch);
        multiSTensorRandomVariableStochastic(FIniInv_PK2_FIniInvT, ndim, ndim, det_FIni, nstoch, C_Stoch, PK2);
        setAll(PK1,0);
        multTensorTensor3Stochastic(FCurr, ndim, ndim, PK2, ndim, PK1, nstoch, C_Stoch);
    }
        
    computeCauchyStoch(FCurr, ndim, nstoch, C_Stoch, PK1, Cauchy);
    //computeStochasticVMS(Cauchy, ndim, nstoch, VMS, C_Stoch);
    double val =0;
    if (flagTanMod) {
        classTensor6* tanModuli = currState.getTangentStochastic();
        for (int iota = 0; iota < nstoch; iota++) {
            for (int i = 0; i < ndim; i++) {
                for (int J = 0; J < ndim; J++) {
                    for (int eps = 0; eps < nstoch; eps++) {
                        for (int k = 0; k < ndim; k++) {
                            for (int M = 0; M < ndim; M++) {
                                val = dfactor1dF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J] - (1/(double(ndim)))*dfactor2dF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J] + dfactor3dF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J] - dfactor4dF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J];
                                tanModuli->setValues(eps, iota, k, M, i, J, val);
                            }
                        }
                    }
                }
            }
        }
    }

}

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young Young modulus
  @param[in] nu Poisson ratio
*/
classStochasticHyperElastic3DCompressibleNeoHookean::classStochasticHyperElastic3DCompressibleNeoHookean(int ndim, string name,
                                                                                           string name_cons, double rho,
                                                                                           double Young, double nu,
                                                                                           vector<double> Stochastic_numbers,
                                                                                           vector<double> Stochastic_parameters,
                                                                                           vector<double> Stochastic_function,
                                                                                           string approximation,
                                                                                           int order, int resolution, bool stochasticMapping, string distribution) : 
                        constitutiveModels(ndim, name, name_cons, rho) {
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

    this->Young = zeros;
    this->nu = zeros;
    this->Young[0] = Young;
    this->nu[0] = nu;
    int offset = 1;
    for (int j = 0; j < Stochastic_parameters.size(); j++) {
        if (Stochastic_parameters[j] == 1) {
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->Young[offset + i] = Young * Stochastic_function[i];
                //count++;
            }
        }
        if (Stochastic_parameters[j] == 2) {
            if(Stochastic_parameters.size()==2){
                offset = offset + Stochastic_numbers[0];
            }
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->nu[offset + i] = nu * Stochastic_function[i + offset - 1];
                offset++;
            }
        }
    }
    vector<double> C(number_stochastic * number_stochastic * number_stochastic, 0);
    if (approximation == "Haar") {
        Build_C_Haar(order, number_func, C);
        ConvertToHaar(this->Young, this->nu, order, Stochastic_parameters, Stochastic_numbers, distribution);
    } else {
        if (approximation == "Mixed"){
            Build_C_Mixed(order, resolution, number_func, C, distribution);
        } else {
            build_ThirdOrderTensor_C(order, number_func, C, distribution);
        }
    }
    C_Stoch = C;
    this->nstoch = number_stochastic;

    vector<double> lambda(number_stochastic, 0);
    vector<double> mu(number_stochastic, 0);
    vector<double> K1(number_stochastic, 0);
    vector<double> twicemu(number_stochastic, 0);
    vector<double> oneminusnu_carre(number_stochastic, 0);
    vector<double> Young_times_nu(number_stochastic, 0);
    vector<double> fac(number_stochastic, 0);
    vector<double> Young_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> nutimesYoung_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> oneplusnu_times_oneminustwonu(number_stochastic, 0);
    vector<double> factor(number_stochastic, 0);
    vector<double> factor_times_nu(number_stochastic, 0);
    vector<double> young_over_oneminustwonu(number_stochastic, 0);
    vector<double> oneminustwicenu = this->nu;
    vector<double> oneplusnu = this->nu;


    linearCombinationVariableStochastic(oneplusnu, 1, 1);
    linearCombinationVariableStochastic(oneminustwicenu, -2, 1);
    DivisionRandomVariableRandomVariableStochastic(this->Young, oneplusnu, C_Stoch, mu);
    vector<double> Young_over_oneplusnu = mu;
    MultiplicationRandomVariableScalar(mu, 0.5);
    
    DivisionRandomVariableRandomVariableStochastic(Young_over_oneplusnu, oneminustwicenu, C_Stoch, factor);
    multiRandomVariableRandomVariableStochastic(factor, this->nu, C_Stoch, factor_times_nu);
    lambda = factor_times_nu;

    DivisionRandomVariableRandomVariableStochastic(this->Young, oneminustwicenu, C_Stoch, young_over_oneminustwonu);
    K1 = young_over_oneminustwonu;
    MultiplicationRandomVariableScalar(K1, 0.33);
    this->lambda = lambda;
    this->mu = mu;
    this->K1 = K1; 

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

}

classStochasticHyperElastic3DCompressibleNeoHookean::~classStochasticHyperElastic3DCompressibleNeoHookean()
{
   
}

void classStochasticHyperElastic3DCompressibleNeoHookean::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term; 
    _term = new MechanicalTermStochastic();;
};

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classStochasticHyperElastic3DCompressibleNeoHookean::soundSpeed() const
{
    return sqrt((lambda[0] + 2 * mu[0]) / _rho);
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classStochasticHyperElastic3DCompressibleNeoHookean::initIntVars(vector<double> &intVars) {
    // This constitutive model does not have internal variables

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classStochasticHyperElastic3DCompressibleNeoHookean::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    classStates& initState = GaussPoint->getInitialState();
    vector<double> &FIni = initState.getDeformationGradient();
    vector<double> &FCurr = currState.getDeformationGradient();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    int ndim = getDim();
    vector<double>FIniInv(nstoch*ndim*ndim,0);
    vector<double>FCurrMecha(nstoch*ndim*ndim,0);
    vector<double>det_FIni(nstoch,0);
    vector<double>Jacobian(nstoch,0);
    vector<double>logJacobian(nstoch,0);
    vector<double>invC(nstoch*ndim*ndim,0);
    vector<double>I_minus_invC(nstoch*ndim*ndim,0);
    vector<double>rightCauchy(nstoch*ndim*ndim,0);
    vector<double>firstPart(nstoch*ndim*ndim,0);
    vector<double>secondPart(nstoch*ndim*ndim,0);
    vector<double>lambda_logJ(nstoch,0);
    
    DeterminantStochasticTensor(FCurr, ndim, nstoch, Jacobian, C_Stoch);
    logIntegral(Jacobian, logJacobian, C_Stoch);
    multiRandomVariableRandomVariableStochastic(lambda, logJacobian, C_Stoch, lambda_logJ);
    multSTensor3FirstTranspose(FCurr, ndim, ndim, FCurr, ndim, rightCauchy);
    InverseTensorStochastic(rightCauchy, ndim, nstoch, C_Stoch, invC);

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
    vector<double> PK2(ndim *ndim *nstoch, 0);
    vector<double> Finv(ndim *ndim *nstoch, 0);

    stress(PK2,PK1,FCurr);


    
    vector<double>FIniInvT(nstoch*ndim*ndim,0);
    vector<double>PK2_FIniInvT(nstoch*ndim*ndim,0);
    vector<double>FIniInv_PK2_FIniInvT(nstoch*ndim*ndim,0);
    if(_stochasticMapping){
        transposeStochastic(FIniInv, ndim, nstoch, FIniInvT);
        multTensorTensor3Stochastic(PK2, ndim, ndim, FIniInvT, ndim, PK2_FIniInvT, nstoch, C_Stoch);
        multTensorTensor3Stochastic(FIniInv, ndim, ndim, PK2_FIniInvT, ndim, FIniInv_PK2_FIniInvT, nstoch, C_Stoch);
        multiSTensorRandomVariableStochastic(FIniInv_PK2_FIniInvT, ndim, ndim, det_FIni, nstoch, C_Stoch, PK2);
        setAll(PK1,0);
        multTensorTensor3Stochastic(FCurr, ndim, ndim, PK2, ndim, PK1, nstoch, C_Stoch);
    }
    
    vector<double>dmuFinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double>dlambdalogJFinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    double tol = 0.00000001;

    

                                    
    computeCauchyStoch(FCurr, ndim, nstoch, C_Stoch, PK1, Cauchy);
    //computeStochasticVMS(Cauchy, ndim, nstoch, VMS, C_Stoch);
    double val = 0;
    if (flagTanMod) {
        build_derivatives_CompressiveNeoHookean(FCurr,dmuFinvTdF,dlambdalogJFinvTdF,lambda_logJ,C_Stoch, tol, nstoch,ndim);
        classTensor6* tanModuli = currState.getTangentStochastic();
        for (int iota = 0; iota < nstoch; iota++) {
            for (int i = 0; i < ndim; i++) {
                for (int J = 0; J < ndim; J++) {
                    for (int eps = 0; eps < nstoch; eps++) {
                        for (int k = 0; k < ndim; k++) {
                            for (int M = 0; M < ndim; M++) {
                                val = dmuFdF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J] - dmuFinvTdF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J] + dlambdalogJFinvTdF[eps*ndim*ndim*ndim*ndim*nstoch + iota*ndim*ndim*ndim*ndim + k*ndim*ndim*ndim + M*ndim*ndim + i*ndim + J];
                                tanModuli->setValues(eps, iota, k, M, i, J, val);
                            }
                        }
                    }
                }
            }
        }
    }

};

/*! \brief Compute derivatives.
  @param[in&out] dJdF
  @param[in&out] dFdF
  @param[in&out] dFdinvT
  @param[in] FCurr
  @param[in] tol
  @param[in] nstoch
 */
void classStochasticHyperElastic3DCompressibleNeoHookean::build_derivatives_CompressiveNeoHookean(const vector<double> &FCurr, vector<double> &dmuFinvTdF, vector<double> &dlambdalogJFinvTdF, vector<double> &lambda_logJ,  vector<double> &C_Stoch, const double tol, const int nstoch, const int ndim) {
    setAll(dmuFinvTdF,0);
    setAll(dlambdalogJFinvTdF,0);

    vector<double> Jacobian(nstoch,0);
    vector<double> Jacobian_Perturbed(nstoch,0);
    vector<double> FPerturbed(ndim*ndim*nstoch,0);
    vector<double> FInv(ndim *ndim *nstoch, 0);
    vector<double> FInvT(ndim *ndim *nstoch, 0);
    vector<double> FInvPerturbed(ndim *ndim *nstoch, 0);
    vector<double> FInvTPerturbed(ndim *ndim *nstoch, 0);
    vector<double> logJ(nstoch,0);
    vector<double> dlogJdF(nstoch*nstoch*ndim*ndim,0);
    vector<double> dlambdalogJdF(nstoch*nstoch*ndim*ndim,0);
    vector<double> one_overJacobian(nstoch,0);
    vector<double> ones(nstoch,0);
    vector<double> Cmu(nstoch*nstoch,0);
    vector<double> logJPerturb(nstoch,0);
    ones[0] = 1;
    vector<double> dJdF(nstoch*nstoch*ndim*ndim,0);
    vector<double> dFinvTdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dFdFFinv(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> CFinv(nstoch*nstoch*ndim*ndim,0);

    // Compute F-T
    InverseTensorStochastic(FCurr, ndim, nstoch, C_Stoch, FInv);
    transposeStochastic(FInv, ndim, nstoch, FInvT);

    FPerturbed = FCurr;
    DeterminantStochasticTensor(FCurr, ndim, nstoch, Jacobian, C_Stoch);
    logIntegral(Jacobian, logJ, C_Stoch);

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
                Cmu[gamma*nstoch + beta] += C_Stoch[gamma*nstoch*nstoch + beta*nstoch + alpha]*mu[alpha];
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
    



    DivisionRandomVariableRandomVariableStochastic(ones, Jacobian, C_Stoch, one_overJacobian);
    derivativeProductRandomConstantRandomVariable(one_overJacobian,dJdF,dlogJdF,C_Stoch,ndim,nstoch);
    derivativeProductRandomConstantRandomVariable(lambda,dlogJdF,dlambdalogJdF,C_Stoch,ndim,nstoch);
    derivativeProductRandomVariableTensorRandomVariable(lambda_logJ, FInvT, dlambdalogJdF, dFinvTdF, dlambdalogJFinvTdF, C_Stoch, ndim, nstoch); 

}

void classStochasticHyperElastic3DCompressibleNeoHookean::stress(vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr) {
    
    int ndim = getDim();
    setAll(PK1,0);
    setAll(PK2,0);
    vector<double>Jacobian(nstoch,0);
    vector<double>logJacobian(nstoch,0);
    vector<double>FInv(nstoch*ndim*ndim,0);
    vector<double>FInvT(nstoch*ndim*ndim,0);
    vector<double>firstPart(nstoch*ndim*ndim,0);
    vector<double>secondPart(nstoch*ndim*ndim,0);
    vector<double>thirdPart(nstoch*ndim*ndim,0);
    vector<double>lambda_logJ(nstoch,0);
    
    DeterminantStochasticTensor(FCurr, ndim, nstoch, Jacobian, C_Stoch);
    logIntegral(Jacobian, logJacobian, C_Stoch);
    InverseTensorStochastic(FCurr, ndim, nstoch, C_Stoch, FInv);
    transposeStochastic(FInv, ndim, nstoch, FInvT);

    multiSTensorRandomVariableStochastic(FCurr, ndim, ndim, mu, nstoch, C_Stoch,
                                                firstPart);

    multiSTensorRandomVariableStochastic(FInvT, ndim, ndim, mu, nstoch, C_Stoch,
                                                thirdPart);                                            
    multiRandomVariableRandomVariableStochastic(lambda, logJacobian, C_Stoch,
                                                 lambda_logJ);                                            
    multiSTensorRandomVariableStochastic(FInvT, ndim, ndim, lambda_logJ, nstoch, C_Stoch,
                                                secondPart);

    for (int i = 0; i < nstoch; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                PK1[i*ndim*ndim + j*ndim +k] = firstPart[i*ndim*ndim + j*ndim +k] - thirdPart[i*ndim*ndim + j*ndim +k] + 1*secondPart[i*ndim*ndim + j*ndim +k];
            }
        }
    }

    multTensorTensor3Stochastic(FInv, ndim, ndim, PK1, ndim, PK2, nstoch, C_Stoch);


}

double classStochasticHyperElastic3DCompressibleNeoHookean::defoEnergy(classGPs *GaussPoint) const
{
    // W = 0.5*lambda*log(J)**2 - mu*log(J)+ 0.5*mu*trace(FT*F)
    //To be computed
    return 0.0;
}

