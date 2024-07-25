//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _boundaryConditions_H_
#define _boundaryConditions_H_

#include "configuration.h"
#include "classGPs.h"
#include <petscsnes.h>
#include <petscksp.h>
#include "GPsDistribution.h"
#include "classNeumannBCs.h"
#include "constitutiveList.h"

extern void BCsNodalForces(vector<classNodes *> &nodes, vector<classElements *> &elements,
                           vector<classNeumann *> &NeumannNodes, int &ndim, vector<double> &force, double timeRun,
                           double dt);

extern void BCsTractionRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classTractionRamp *&BCsTraction, int &ndim, vector<double> &force, double timeRun, double dt);

extern void BCsTractionInst(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classTractionInst *&BCsTraction, int &ndim, vector<double> &force, double timeRun, double dt);

extern void BCsStressRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classStressRamp *&BCsTraction, int &ndim, vector<double> &force, double timeRun, double dt);

extern void BCsStressInst(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classStressInst *&BCsTraction, int &ndim, vector<double> &force, double timeRun, double dt);

extern void BCsStochasticPressureRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classStochasticPressureRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt);

extern void BCsPressureRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classPressureRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt,
                            bool tangentEstimation, map<pair<int,int>,double>& stiffBC);

extern void BCsHertzianRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classHertzianRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt, 
                             bool tangentEstimation, map<pair<int,int>,double>& stiffBC);

extern void BCsStochasticHertzianRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classStochasticHertzianRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt);

extern void BCsPressureInst(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classPressureInst *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt);

extern void BCsGravity(vector<classNodes *> &nodes, vector<classElements *> &elements, int &ndim, vector<double> &force,
                       classGravity *&BCsGrav);
                       
extern void BCsVolHeatFlux(vector<classGPs *> &GPs,  vector<classNodes *> &nodes, classVolHeatFlux *&BCsVolHeatFlux, int &ndim,
                       vector<double>& force, double timeRun, double dt);
                       
extern void BCsVolHeatFluxGaussian(const vector<double> &uK, vector<classGPs *> &GPs,  vector<classNodes *> &nodes, classVolHeatFluxGaussian *&BCsVolHeatFlux, int &ndim,
                       vector<double>& force, double timeRun, double dt);
                       
extern void BCsCurrentInst(vector<classGPs *> &GPs, vector<classNodes *> &nodes, classCurrentInst *&BCsCurrent, int &ndim,
                        vector<double>& force, double timeRun, double dt);
                        
extern void BCsHeatFluxInst(const vector<double> &uK, vector<double>& force,  vector<classNodes *> &nodes,
                            vector<classElements *> &elements, classSurfaceHeatFluxInst *&BCsHeat, int &ndim, double timeRun,
                            double dt) ;
extern void BCsConvection(const vector<double> &uK, vector<double>& force, vector<classNodes *> &nodes,
                          vector<classElements *> &elements, classConvection *&BCsConvection, int &ndim, double timeRun,
                          double dt,  bool tangentEstimation, map<pair<int,int>,double>& stiffBC);
extern void BCsRadiation(const vector<double> &uK, vector<double>& force, vector<classNodes *> &nodes,
                         vector<classElements *> &elements, classRadiation *&BCsRadiation, int &ndim, double timeRun,
                         double dt,  bool tangentEstimation, map<pair<int,int>,double>& stiffBC);

extern void NeumannBCManagement(const vector<double> &uK, vector<classGPs *> &GPs, vector<classNodes *> &nodes, vector<classElements *> &element,
                                vector<classNeumannBCs *> &NeumannBCs, int &ndim, vector<double> &force, double timeRun,
                                double dt, bool tangentEstimation, map<pair<int,int>,double>& stiffPart);

extern void deleteNeumannBCs(vector<classNeumannBCs *> &NeumannBCs);

extern void findActiveDof(vector<classDirichlet *> DirichletNodes, int ndim, int numDof, vector<int> &activeDof);

extern void activeSearching(vector<classDirichlet *> DirichletNodes, int ndim, int numDof, vector<int> &activeDof);

#ifdef PARALLEL
extern void findActiveDof(vector<classDirichlet*> DirichletNodes, int ndim, int numDof, vector<int> &activeDof, vector<int> nod_local, vector<int> &nDof_local, int nNod);
#endif
#endif
