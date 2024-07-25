//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _stochasticHyperElasticStVenantKirchhoff_H_
#define _stochasticHyperElasticStVenantKirchhoff_H_

#include "constitutiveModels.h"
#include "maths.h"

/*! \brief Hyper Elastic St-Venant - Kirchhoff constitutive model.
*/
class classStochasticHyperElasticStVenantKirchhoff : public constitutiveModels {
    protected:
        classTensor6 *elasticStiffTensor;
        classTensor6 *_tanModuli;
        vector<double> Young;
        vector<double> C;
        vector<double> nu;
        vector<double> _fac;
        vector<double> _nutimesYoung_over_one_minusnu_carre;
        vector<double> _mu2;
        vector<double> _lam;
        int nstoch;
        vector<double> C_Stoch;
        int _stressState;//0 for plain stress, 1 for plain strain
        bool _stochasticMapping;
        
    public:
        classStochasticHyperElasticStVenantKirchhoff(int ndim, string name, string name_cons, double rho, double young,
                                                     double nu, int stressState, vector<double> Stochastic_numbers,
                                                     vector<double> Stochastic_parameters,
                                                     vector<double> Stochastic_function, string approximation, int order, int resolution, bool stochasticMapping, string distribution);
        virtual ~classStochasticHyperElasticStVenantKirchhoff();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return nstoch*getDim();}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{};
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{};
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
        
        const vector<double>& getC_Stoch() const { return C_Stoch; }
        void stress(vector<double> &Sig, const vector<double> &E);
        void returnTangentModulus(classTensor6 *tangentModulus){tangentModulus = _tanModuli;}
        void getTanModuli(classTensor6* tanModuli, const vector<double> &FCurr, const vector<double> &PK2, const bool flagTanMod);
        
    private:
        void build_CF(vector<double> &CF, const vector<double> &FCurr);
        void build_dEdF(vector<double> &dEdF, const vector<double> &FCurr);
        void build_dSdF(vector<double> &dSdF,vector<double> &dEdF);
        void fillElasticStiffTensor(int ndim, int n_stoch);
        void pushBackdEdF(vector<double> &dEdF, const vector<double> &FIniInv);
        void pushBackdSdF(vector<double> &dSdF, const vector<double> &FIniInv, const vector<double> &FIniInvT, const vector<double> &det_FIni);
       
};

#endif
