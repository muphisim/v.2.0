//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classElements_H
#define _classElements_H

#include "configuration.h"
#include "constitutiveModels.h"
#include "geometry.h"
#include "classNodes.h"
#include "classGPs.h"
#include "GPsDistribution.h"

/*! \brief  This is the wrapper for all type of elements. This base class implements general and virtual functions. ALL VIRTUAL FUNCTIONS MUST BE IMPLEMENTED BY THE USER FOR ALL NEW TYPE OF ELEMENTS*/
class classElements {

protected:

    int me; /*!< \brief Who are you. Number of the element. It begins from zero */
    int bulkElement; //in case of subelement, number of parent element
    int localNdim;
    vector<constitutiveModels *> constitutiveModel;/*!< \brief Constitutive model for the element*/
    string type; /*!< \brief String with the type of element*/
    vector<int> myNodes;/*!< \brief This array stores the nodes that are part of the element*/
    vector<double> Cauchy; /*!< \brief Cauchy stress of the element*/
    vector<double> E;/*!< \brief Green strain of the element*/
    double charLength;/*!< \brief Characteristic length of the element*/
    vector<int> myGPs;/*!< \brief This array the labels of the GPs that are inside the element */
    vector <double> VMS;
    vector <double> volumetricStrain;
    vector <double> equivalentStrain;
    vector<double> n0;/*!< \brief anisotropic direction of conductivity tensor */
    vector<double> intVars; /*!< \brief Array with the internal variables of the model (to be plotted)*/
    bool _flagExtra;
    bool _flagStochastic;
    vector<int> intVarsPosition;/*!< \brief This array stores the internal variable in which a value is going to be stored*/
    vector<double> intVarsValue;/*!< \brief This array stores the value of a internal variable*/
public:

    classElements();

    classElements(int me, string type, const vector<int>& myNodes, int ndim, int bulkElement);

    inline string getType() const { return type; }; // Common

    inline void setConsModel(vector<constitutiveModels *>& constitutiveModel) {
        this->constitutiveModel = constitutiveModel;
    }; // Common

    inline void setConsModel(constitutiveModels *constitutiveModel) {
        this->constitutiveModel.push_back(constitutiveModel);
    }; // Common

    inline vector<constitutiveModels *>& getConsModel() {
        return constitutiveModel;
    }; // Common
    inline int getBulkElement() const { return bulkElement; }

    inline void setFlagExtra(bool flagExtra) {
        _flagExtra = flagExtra;
    }

    inline void setFlagStochastic(bool flagStochastic) {
        _flagStochastic = flagStochastic;
    }

    inline bool getFlagExtra() const {
        return _flagExtra;
    }

    inline bool getFlagStochastic() const{
        return _flagStochastic;
    }
    
    inline const vector<int>& getMyNodes() const { return myNodes; };// Common
    inline const vector<double>& getCauchy() const { return Cauchy; };// Common
    inline const vector<double>& getE() const { return E; };// Common
    inline int getMe() const { return me; };// Common
    inline void setMe(int me) {
        this->me = me;
    };// Common
    inline void sumCauchy(const vector<double>& Cauchy) {
        int size = Cauchy.size();
        int size_elem = (this->Cauchy).size();
        if(size != size_elem){
           (this->Cauchy).resize(size,0);
        }
        for (int i = 0; i < size; i++) {
            this->Cauchy[i] += Cauchy[i];
        }
    };// Common

    inline const vector<int>& getMyGPs() const{ return myGPs; };// Common
    inline void addMyGPs(int ind) {
        myGPs.push_back(ind);
    };// Common
    inline void sumE(const vector<double>& E){
        int size = E.size();
        int size_elem = (this->E).size();
        if(size != size_elem){
           (this->E).resize(size,0);
        }
        for (int i = 0; i < size; i++) {
            this->E[i] += E[i];
        }
    };// Common
    inline void resetE() {
        fill(this->E.begin(), this->E.end(), 0);
    }

    inline void resetCauchy() {
        fill(this->Cauchy.begin(), this->Cauchy.end(), 0);
    }

    inline double getCharLength() const{
        return this->charLength;
    };

    inline void sumVMS(vector<double>& VMS) {
        int size = VMS.size();
        int size_elem = (this->VMS).size();
        if(size != size_elem){
           (this->VMS).resize(size,0);
        }
        for (int i = 0; i < size; i++) {
            this->VMS[i] += VMS[i];
        }
    };// Common
    inline void resetVMS() {
        fill(this->VMS.begin(), this->VMS.end(), 0);
    };

    inline vector <double> getVMS() const { return VMS; };// Common

    inline void sumVolumetricStrain(vector<double>& volumetricStrain) {
        int size = volumetricStrain.size();
        int size_elem = (this->volumetricStrain).size();
        if(size != size_elem){
           (this->volumetricStrain).resize(size,0);
        }
        for (int i = 0; i < size; i++) {
            this->volumetricStrain[i] += volumetricStrain[i];
        }
    };// Common
    inline void resetVolumetricStrain() {
        fill(this->volumetricStrain.begin(), this->volumetricStrain.end(), 0);
    };

    inline vector <double> getVolumetricStrain() const { return volumetricStrain; };// Common

    inline void sumEquivalentStrain(vector<double>& equivalentStrain) {
        int size = equivalentStrain.size();
        int size_elem = (this->equivalentStrain).size();
        if(size != size_elem){
           (this->equivalentStrain).resize(size,0);
        }
        for (int i = 0; i < size; i++) {
            this->equivalentStrain[i] += equivalentStrain[i];
        }
    };// Common
    inline void resetEquivalentStrain() {
        fill(this->equivalentStrain.begin(), this->equivalentStrain.end(), 0);
    };

    inline vector <double> getEquivalentStrain() const { return equivalentStrain; };// Common

    inline void set_n0(const vector<double>& a0) { this->n0 = a0; };

    inline const vector<double>& get_n0() const { return this->n0; };

    inline void sumIntVars(const vector<double>& intVars2) {
        for (int i = 0; i < intVars2.size(); i++) {
            this->intVars[i] += intVars2[i];
        }
    };// Common
    inline void initIntVars(int size) {
        this->intVars.resize(size, 0);
    };

    inline void resetIntVars() {
        this->intVars.clear();
    };

    inline const vector<double>& getIntVars()  const{ return this->intVars; };// Common

    inline void add_IntVarsPosition(int position) {
        this->intVarsPosition.push_back(position);
    };

    inline const vector<int>& get_IntVarsPosition() const { return this->intVarsPosition; };

    inline void add_IntVarsValue(double value) {
        this->intVarsValue.push_back(value);
    };

    inline const vector<double>& get_IntVarsValue() const { return this->intVarsValue; };
    
    
    virtual void uvwToxyz(const vector<classNodes *> &nodes, int ndim, const vector<double>& uvw, vector<double>& xyz) const;
    virtual void xyzTouvw(const vector<classNodes *> &nodes, int ndim, const vector<double>& xyz, vector<double>& uvw) const;
    virtual void getShapeFunctionsXYZ(const vector<classNodes *> &nodes, int ndim, const vector<double>& uvw, vector<double>& f) const;
    virtual void getGradShapeFunctionsXYZ(const vector<classNodes *> &nodes, int ndim, const vector<double>& uvw, vector<double>& gradf) const;    
    virtual double getJacobian(const vector<classNodes *> &nodes, int ndim, const vector<double>& uvw, mMatrix& jac) const;
    
    virtual void defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim);
    virtual void defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int bulkElement);
    virtual void definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int orderd, vector<classGPs *> &mmGPs);
    virtual void defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim, const vector<double>& Point);
    
    // pure virtual functions
    /*! This function gets the default integration order for FEM
    */
    virtual int getDefaultIntegrationOrder() const = 0;
    /*! This function gets the default integration order for MM
    */
    virtual int getDefaultIntegrationOrderMM() const  = 0;
    /*! This function gets the integration points over elements
      @param[in] order Integration order
      @param[out] npts number of integration points
      @param[out] uvw Isoparametric coordinates
      @param[out] weight all weights
    */
    virtual void getIntegrationPoints(int order, int& npts, vector<vector<double> >& uvw, vector<double>& weight) const =0;
    
    /*! This function gets the shape function in the isoparametric space
      @param[in] uvw Isoparametric coordinates
      @param[out] fuvw Shape functions
    */
    virtual void getShapeFunctionsUVW(const vector<double>& uvw, vector<double>& fuvw) const = 0;
     /*! This function gets the shape function in the isoparametric space
      @param[in] uvw Isoparametric coordinates
      @param[out] gradfuvw Grad of shape functions
    */
    virtual void getGradShapeFunctionsUVW(const vector<double>& uvw, vector<double>& gradfuvw) const = 0;


    /*! This function calculates the normal of the segment based on the two points
      @param[in] nod Array with all nodes in the domain
      @param[in] ndim Dimension of the domain
    */
    virtual vector<double> getn0(vector<classNodes *>& nod, int ndim) = 0;

    /*! This function calculates a subelement from a given element. This function should not be called for segments
      @param[in] me Label of the element
      @param[in] S label of the edge (2D) or surface (3D) that is the subelement
      @param[out] subElement Subelement
      @param[in] nod Array with all nodes in the domain
    */
    virtual void getSubElement(int me, int S, vector<classElements *>& subElements, vector<classNodes *>& nod, int bulkElement) =0;

    /*! This function calculates the characteristic lenght of the element. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
      @param[in] nod Array with all nodes in the domain
      @param[in] ndim Dimension of the domain
    */
    virtual double characteristicLength(int ndim, const vector<vector<double> >& nodesVec) = 0;

    /*! Destructor
     */
    virtual ~classElements() {
    };
};

#include "elementsList.h"

#endif
