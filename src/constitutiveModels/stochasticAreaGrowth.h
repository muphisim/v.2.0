//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _stochasticAreaGrowth_H_
#define _stochasticAreaGrowth_H_

#include "constitutiveModels.h"

/*! \brief Ellen Kuhl formulation for the area growth, Journal of Mechanical Behavior of Biomedical Materials, 29 (2014), 529-543.
Decomposition of the deformation gradient into the elastic (Hyperelastic) and growth contribution: F=Fe*Fg
*/
class classStochasticAreaGrowth : public constitutiveModels {
    protected:
        classTensor6 *elasticStiffTensor;
        vector<double> _Young;
        vector<double> _nu;
        vector<double> _lambda;
        vector<double> _mu;
        vector<double> _Gc; // Same notation than Ellen Kuhl
        int nstoch;
        vector<double> C_Stoch;
        vector<double> dFdF;
        vector<double> dnXn0dF;
        bool _stochasticMapping;
        vector<double> _n0;
        vector<double> _n0Xn0;
    public:
        classStochasticAreaGrowth(int ndim, string name, string name_cons, double rho, double Young, double nu,
                                    double Gc,vector<double> Stochastic_numbers,
                                                     vector<double> Stochastic_parameters,
                                                     vector<double> Stochastic_function, string approximation, int order, int resolution, bool stochasticMapping, string distribution);
        virtual ~classStochasticAreaGrowth();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return nstoch*getDim();}; 
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{}
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{}
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);  
        const vector<double>& getC_Stoch() const { return C_Stoch; }


    private:
        void computeFe(const vector<double> &FCurr, const vector<double> &sqrt_theta, vector<double> &Fe);
        void stress(vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr, const vector<double> &sqrt_theta, const vector<double> &FginvT);
        void build_derivatives(const vector<double> &FCurr, vector<double> &dPK1dF, const vector<double> &FgInvT, const vector<double> &theta, const vector<double> &C_Stoch, const double tol, const int nstoch, const int ndim);
};

#endif
