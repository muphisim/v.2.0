//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classElements.cpp
  \brief This file contains the definition of the constructors of the classElements
*/

#include "classElements.h"

/*! Constructor for classElements. It should not be called.

*/
classElements::classElements() {
    ERROR("This constructor should not be called");
    exit(-1);
};

/*! Constructor for classElements in case of subelements.
  @param[in] me Number of element
  @param[in] type Type of element
  @param[in] myNodes Nodes that build the element
  @param[in] ndim Dimension of the domain
  @param[in] bulkElement number of parent element to subelement
*/
classElements::classElements(int me, string type, const vector<int>& myNodes, int ndim, int bulkElement) {

    this->bulkElement = bulkElement;
    this->me = me;
    this->type = type;
    this->myNodes = myNodes;
    this->Cauchy.assign(ndim * ndim, 0);
    this->E.assign(ndim * ndim, 0);
    this->VMS.assign(1, 0);
    this->volumetricStrain.assign(1, 0);
    this->equivalentStrain.assign(1, 0);
    this->n0.assign(3, 0);
    this->localNdim = ndim;
};

/*! convert isoparametric coordinates to physical coordinates
  @param[in] nodes nodes list to access their coordinates
  @param[in] ndim Dimension of the domain
  @param[in] uvw Isoparametric coordinates
  @param[out] xyz physical coordinates
*/
void classElements::uvwToxyz(const vector<classNodes *> &nodes, int ndim, const vector<double>& uvw, vector<double>& xyz) const
{
    int numShapeFunctions = myNodes.size();
    xyz.resize(ndim);
    setAll(xyz,0);
    vector<double> sf(numShapeFunctions,0.);
    getShapeFunctionsUVW(uvw,sf);
    
    for (int i=0; i< numShapeFunctions; i++)
    {
        const vector<double>&  xyz_node = nodes[myNodes[i]]->getXYZ();
        for (int j=0; j< ndim; j++)
        {
            xyz[j]+= sf[i]*xyz_node[j];
        }
    }
};

/*! convert physical coordinates to isoparametric coordinates
  @param[in] nodes nodes list to access their coordinates
  @param[in] ndim Dimension of the domain
  @param[in] xyz physical coordinates
  @param[out] uvw Isoparametric coordinates
*/
void classElements::xyzTouvw(const vector<classNodes *> &nodes, int ndim, const vector<double>& xyz, vector<double>& uvw) const
{
    /* NR iterations to solve
     *  xyz - sum_i N(u,v,w)*xyz_i = 0 to find uvw
     * */
     
    if (ndim != localNdim)
    {
        ERROR("classElements::xyzTouvw can be called for bulk elemnts");
        exit(-1);
    };
     
    uvw.resize(ndim,0);
    setAll(uvw,0);
    vector<double> xyz_iter(ndim,0);
    uvwToxyz(nodes,ndim,uvw,xyz_iter);
    mMatrix jac(ndim, ndim);
    mMatrix invjac(ndim, ndim);
    int iter = 1, maxiter = 20;
    double error = 1., tol = 1.e-6;
    while(error > tol && iter < maxiter) 
    {
        double detJ = getJacobian(nodes, ndim, uvw, jac);
        jac.invert(invjac);        
        for (int i=0; i< ndim; i++)
        {
            for (int j=0; j< ndim; j++)
            {
                uvw[i] += invjac(i,j)*(xyz[j] - xyz_iter[j]);
            }
        }
        uvwToxyz(nodes,ndim,uvw,xyz_iter);
        error = 0;
        for (int i=0; i< ndim; i++)
        {
            error += (xyz[i] -xyz_iter[i])*(xyz[i]-xyz_iter[i]);
        }
        error = sqrt(error);        
        iter++;
    };
};

/*! get shape functions at a point given by isoparametric coordinates
  @param[in] nodes nodes list to access their coordinates
  @param[in] ndim Dimension of the domain
  @param[in] uvw Isoparametric coordinates
  @param[out] f Shape functions
*/
void classElements::getShapeFunctionsXYZ(const vector<classNodes *> &nodes, int ndim, const vector<double>& uvw, vector<double>& f) const
{
    // the same as the ones in the isoparametric space
    getShapeFunctionsUVW(uvw,f);
};

/*! get grad of shape function at a point given by isoparametric coordinates
  @param[in] nodes nodes list to access their coordinates
  @param[in] ndim Dimension of the domain
  @param[in] uvw Isoparametric coordinates
  @param[out] gradf Grad of shape functions
*/
void classElements::getGradShapeFunctionsXYZ(const vector<classNodes *> &nodes, int ndim, const vector<double>& uvw, vector<double>& gradf) const
{
    int numShapeFunctions = myNodes.size();
    vector<double> fgraduvw(numShapeFunctions*localNdim, 0); // grad uvw given in local dimension
    getGradShapeFunctionsUVW(uvw,fgraduvw);
    gradf.resize(numShapeFunctions*ndim);
    setAll(gradf,0);
    mMatrix jac(ndim,ndim), invJac(ndim,ndim);
    double detJ = getJacobian(nodes, ndim, uvw, jac); 
    jac.invert(invJac);
    for (int i=0; i< numShapeFunctions; i++)
    {
        const vector<double>&  xyz_node = nodes[myNodes[i]]->getXYZ();
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< localNdim; k++)
            {
                gradf[i*ndim+j] += fgraduvw[i*localNdim +k]*invJac(k,j);
            }
        }
    }
};

/*! get Jacobian at a point given by isoparametric coordinates
  @param[in] nodes nodes list to access their coordinates
  @param[in] ndim Dimension of the domain
  @param[in] uvw Isoparametric coordinates
  @param[out] jac Jacobian jac_kl = partial x_k/ partial xi_l
  @return Jacobian determinant 
*/
double classElements::getJacobian(const vector<classNodes *> &nodes, int ndim, const vector<double>& uvw, mMatrix& jac) const
{
    jac.resize(ndim, ndim, true);
    int numShapeFunctions = myNodes.size();
    vector<double> gradfuvw(numShapeFunctions*localNdim,0.);
    getGradShapeFunctionsUVW(uvw,gradfuvw); // grad in isoparametric space given in localNdim
    for (int i=0; i< numShapeFunctions; i++)
    {
        const vector<double>&  xyz_node = nodes[myNodes[i]]->getXYZ();
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< localNdim; k++)
            {
                jac(j,k) += gradfuvw[i*localNdim+k]*xyz_node[j];
            }
        }
    }
    double detJ = 0;
    auto prod_vec = [](const vector<double>& a, const vector<double>& b, vector<double>& c)
    {
      c.resize(3);
      c[2] = a[0] * b[1] - a[1] * b[0];
      c[1] = -a[0] * b[2] + a[2] * b[0];
      c[0] = a[1] * b[2] - a[2] * b[1];
    };
    
    if (localNdim < ndim)
    {
        if (localNdim == 0)
        {
          // a point
          detJ = 1;
          jac.setAll(0);
          for (int k=0; k< ndim; k++)
          {
            jac(k,k)=1.;
          }
          
        }
        else if (localNdim == 1)
        {
          detJ =0.;
          for (int k=0; k< ndim; k++)
          {
            detJ += jac(k,0)*jac(k,0);
          }
          detJ = sqrt(detJ);
          
          // regularize matrix
          vector<double> a(3), b(3);
          for (int k=0; k< ndim; k++)
          {
            a[k] = jac(k,0);
          }
          if((fabs(a[0]) >= fabs(a[1]) && fabs(a[0]) >= fabs(a[2])) || 
                (fabs(a[1]) >= fabs(a[0]) && fabs(a[1]) >= fabs(a[2]))) 
          {
            // if a0 or a1 larger than a2
            b[0] = a[1];
            b[1] = -a[0];
            b[2] = 0.;
          }
          else 
          {
            // a2 is the max value
            b[0] = 0.;
            b[1] = a[2];
            b[2] = -a[1];
          }
          double normb = norm2(b,3);
          scale(b,1./normb);
          for (int k=0; k< ndim; k++)
          {
              jac(k,1) = b[k];
          }
          
          if (ndim == 3)
          {
              vector<double> c(3);
              prod_vec(a, b, c);
              double normc = norm2(c,3);
              scale(c,1./normc);
              for (int k=0; k< ndim; k++)
              {
                  jac(k,2) = c[k];
              }
          }          
        }
        else if (localNdim == 2) 
        {
            auto square = [](double x){return x*x;};
            detJ = sqrt(square(jac(0,0) * jac(1,1) - jac(0,1)* jac(1,0)) +
              square(jac(2,0) * jac(0,1) - jac(0,0) * jac(2,1)) +
              square(jac(1,0) * jac(2,1) - jac(2,0) * jac(1,1)));
            // regularize matrix
            vector<double> a(3), b(3), c(3);
            a[0] = jac(0,0);
            a[1] = jac(1,0);
            a[2] = jac(2,0);
            b[0] = jac(0,1);
            b[1] = jac(1,1);
            b[2] = jac(2,1);
            prod_vec(a, b, c);
            double normc = norm2(c,3);
            scale(c,1./normc);
            jac(0,2) = c[0];
            jac(1,2) = c[1];
            jac(2,2) = c[2];
        }
    }
    else
    {
      detJ = jac.det();
    }
    return detJ;
};


/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM.
      @param[in] ind Label of th element
      @param[inout] GPs Array with all GPs in the domain
      @param[in] nodes Array with all nodes in the domain
      @param[in] ndim Dimension of the domain
    */
void classElements::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim)
{
    // clear existing GP lists
    myGPs.clear();
    int nGPs=0;
    vector<vector<double> > uvwAll; // isoparametric space
    vector<double> wAll; // weights
    int order = getDefaultIntegrationOrder();
    getIntegrationPoints(order, nGPs,uvwAll,wAll);
    mMatrix jac(ndim,ndim);
    vector<double> xyzIn(ndim,0.);
    for (int i=0; i< nGPs; i++)
    { 
        double J = getJacobian(nodes, ndim, uvwAll[i], jac);
        uvwToxyz(nodes,ndim,uvwAll[i],xyzIn);
        classGPs *tempGP = new classGPs(GPs.size(), xyzIn, wAll[i], J, ndim);
        tempGP->setNeighbours(myNodes, myNodes.size());
        
        vector<double> phi, dphi;
        getShapeFunctionsXYZ(nodes, ndim, uvwAll[i],phi);
        getGradShapeFunctionsXYZ(nodes, ndim, uvwAll[i], dphi);
        tempGP->setShapeFunctions(phi, dphi);
        tempGP->setCell(ind, this->getType());
        tempGP->setConsModel(constitutiveModel);
        myGPs.push_back(GPs.size());
        GPs.push_back(tempGP);
    };
};

/*! This function defines just the positions, weights and Jacobians of the GPs inside the element. 
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] orderd This value defines the number of GPs inside the element (dramatically important for MM simulations)
  @param[inout] mmGPs Array with all GPs in the MM domain
*/
void classElements::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int bulkElement)
{
    if (localNdim >= ndim)
    {
        ERROR("defineGPs with bulk element must be defined in sub-elements");
        exit(-1);
    }
    int currentIndex = GPs.size();
    defineGPs(ind, GPs, nodes, ndim);
    for (int j = currentIndex; j < GPs.size(); j++) 
    {
        GPs[j]->setBulkElement(bulkElement);
    }
};

/*! This function defines just the positions, weights and Jacobians of the GPs inside the element.
    @param[in] ind Label of th element
    @param[inout] GPs Array with all GPs in the domain
    @param[in] nodes Array with all nodes in the domain
    @param[in] ndim Dimension of the domain
    @param[in] orderd This value defines the number of GPs inside the element (dramatically important for MM simulations)
    @param[inout] mmGPs Array with all GPs in the MM domain
  */
void classElements::definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int orderd,
                        vector<classGPs *> &mmGPs)
{
    // clear existing GP lists
    myGPs.clear();
    int nGPs=0;
    vector<vector<double> > uvwAll; // isoparametric space
    vector<double> wAll; // weights
    if (orderd >0)
    {
        // constrained order
        getIntegrationPoints(orderd,nGPs,uvwAll,wAll);
    }
    else
    {
        // defaut order
        int order = getDefaultIntegrationOrderMM();
        getIntegrationPoints(order,nGPs,uvwAll,wAll);
    }
    
    mMatrix jac(ndim,ndim);
    vector<double> xyzIn(ndim,0.);
    for (int i = 0; i < nGPs; i++) 
    {
        double J = getJacobian(nodes, ndim, uvwAll[i], jac);
        uvwToxyz(nodes,ndim,uvwAll[i],xyzIn);
        classGPs *tempMMGPs = new classGPs(GPs.size(), xyzIn, wAll[i], J, ndim);        
        tempMMGPs->setCell(ind, this->getType());
        tempMMGPs->setConsModel(constitutiveModel);
        myGPs.push_back(GPs.size());
        mmGPs.push_back(tempMMGPs);
        GPs.push_back(tempMMGPs);
    };
};


/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM.
  @param[in] ind Label of th element
  @param[inout] tempGP Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] Point Point to be evaluated
*/
void classElements::defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim, const vector<double>& Point)
{
    vector<double> uvw;
    xyzTouvw(nodes,ndim,Point,uvw);
    tempGP = new classGPs(ind, Point, 1., 1., ndim);  // weight and jacobian will not be used
    tempGP->setNeighbours(myNodes, myNodes.size());
    tempGP->setCell(ind, this->getType());
    tempGP->setConsModel(constitutiveModel);
    vector<double> phi, dphi;
    getShapeFunctionsXYZ(nodes, ndim, uvw, phi);
    getGradShapeFunctionsXYZ(nodes, ndim, uvw, dphi);
    tempGP->setShapeFunctions(phi, dphi);
};
