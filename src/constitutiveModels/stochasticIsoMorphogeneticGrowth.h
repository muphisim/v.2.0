//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _stochasticIsoMorphogeneticGrowth_H_
#define _stochasticIsoMorphogeneticGrowth_H_

#include "constitutiveModels.h"

/*! \brief Ellen Kuhl formulation for the area growth, Journal of Mechanical Behavior of Biomedical Materials, 29 (2014), 529-543.
Decomposition of the deformation gradient into the elastic (Hyperelastic) and growth contribution: F=Fe*Fg
*/
class classStochasticIsoMorphogeneticGrowth : public constitutiveModels {
    protected:
        classTensor6 *elasticStiffTensor;
        vector<double> _Young;
        vector<double> _nu;
        vector<double> _lambda;
        vector<double> _mu;
        vector<double> _Gc; // Same notation than Ellen Kuhl
        int nstoch;
        vector<double> C_Stoch;
        vector<double> dmuFdF;
        bool _stochasticMapping;
    public:
        classStochasticIsoMorphogeneticGrowth(int ndim, string name, string name_cons, double rho, double Young, double nu,
                                    double Gc,vector<double> Stochastic_numbers,
                                                     vector<double> Stochastic_parameters,
                                                     vector<double> Stochastic_function, string approximation, int order, int resolution, bool stochasticMapping, string distribution);
        virtual ~classStochasticIsoMorphogeneticGrowth();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return nstoch*getDim();}; 
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{}
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{}
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);  
        const vector<double>& getC_Stoch() const { return C_Stoch; }

        // For external calls
        void getFg(const vector<double> &FCurr, const vector<double> &FPrev, const vector<double> &intVarsPrev,
                    vector<double> &intVarsCurr, double dt, double timeRun, bool flagTanMod,
                    const vector<double> &PK1, classTensor4 *tanModuli, vector<double> &Fg) const ;
        void getE(int ndim, vector<double> &E, const vector<double> &Fcurr, const vector<double> &Fprev, const vector<double> &intVars,
                  double dt, double timeRun, int nStep) const;

    private:
        void stress(vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr, const vector<double> &Theta);
        void build_derivatives(const vector<double> &FCurr, vector<double> &dmuFinvTdF, vector<double> &dlambdalogJFinvTdF, vector<double> &dmuthetapowerdF, vector<double> &theta, vector<double> &C_Stoch, const double tol, const int nstoch, const int ndim);
};

#endif
