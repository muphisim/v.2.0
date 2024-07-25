//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _stochasticViscoPlasticity_H_
#define _stochasticViscoPlasticity_H_

#include "constitutiveModels.h"

/*! \brief Ellen Kuhl formulation for the area growth, Journal of Mechanical Behavior of Biomedical Materials, 29 (2014), 529-543.
Decomposition of the deformation gradient into the elastic (Hyperelastic) and growth contribution: F=Fe*Fg
*/
class classStochasticViscoPlasticity : public constitutiveModels {
    protected:
        constitutiveModels *HyperElasticStVenantKirchhoff;
        classTensor6 *elasticStiffTensor;
        vector<double> _s;
        double _n;
        vector<double> _eps_0;
        double _Young;
        int nstoch;
        vector<double> C_Stoch;
        bool _stochasticMapping;
    public:
        classStochasticViscoPlasticity(int ndim, string name, string name_cons, double rho, double Young, double nu, double eps_0, double s, double n, int stressState, vector<double> Stochastic_numbers,
                                                     vector<double> Stochastic_parameters,
                                                     vector<double> Stochastic_function, string approximation, int order, int resolution, bool stochasticMapping, string distribution);
        virtual ~classStochasticViscoPlasticity();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return nstoch*getDim();}; 
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{}
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{}
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);  
        const vector<double>& getC_Stoch() const { return C_Stoch; }


    private:
        void computePlasticIncrement(vector<double> &intVarsCurr, vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr, const double dt);
        void computeSDev(vector<double> &SDev, vector<double> &C, const vector<double> &FeCurr, const vector<double> &S);
        void computeSigEq(const vector<double> &SDev, const vector<double> &C, vector<double>&SigEv);
        void computeR(const vector<double> &SDev, const vector<double> &C, const vector<double>&SigEv, vector<double>&R);
        void computeEpsP(vector<double> &EpsP, const vector<double> &SigEv);
        void computeRes(const vector<double> &FPCurr, const vector<double> &FPInvCurr, const vector<double> &FPPrev, const vector<double> &Eps_P_R, vector<double> &ResVector, double &Res, double dt);

};

#endif
