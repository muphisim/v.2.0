//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classLinearHexahedra_H
#define _classLinearHexahedra_H

#include "classElements.h"

/*! \brief  This class is for linear hexahedra */
class classLinearHexahedra : public classElements {

public:

    classLinearHexahedra(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
                          int bulkElement);
  
    virtual int getDefaultIntegrationOrder() const {return 2;};
    virtual int getDefaultIntegrationOrderMM() const {return 2;};
    virtual void getIntegrationPoints(int order, int& npts, vector<vector<double> >& uvw, vector<double>& weight) const;
    virtual void getShapeFunctionsUVW(const vector<double>& uvw, vector<double>& fuvw) const;
    virtual void getGradShapeFunctionsUVW(const vector<double>& uvw, vector<double>& gradfuvw) const;
  
    virtual vector<double> getn0(vector<classNodes *>& nod, int ndim);

    virtual void getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod, int bulkElement);

    virtual double characteristicLength(int ndim, const vector<vector<double> >& nodesVec);
    
    ~classLinearHexahedra() {
    };
};


#endif
