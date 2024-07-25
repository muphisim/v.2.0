//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _poroElastic_H_
#define _poroElastic_H_

#include "constitutiveModels.h"
#include "linearElastic.h"

class DarcyLaw : public constitutiveModels 
{
    protected:
      double _kappa, _mu, _Q, _alpha;
      int _intVarsSize;
        
    public:
        DarcyLaw(int ndim, string name, string name_cons, double kappa, double mu, double Q,  double alpha );
        virtual ~DarcyLaw(){};
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return 1;}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{};
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{};
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod); 
};

class LinearPoroelastic : public linearElastic 
{
    protected:
        double _alpha;

    public:
        LinearPoroelastic(int ndim, string name, string name_cons, double rho, double young, double nu, double alpha, int stressState);
        virtual ~LinearPoroelastic(){}
       
        virtual void initLawFromOptions();
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
}; 

#endif
