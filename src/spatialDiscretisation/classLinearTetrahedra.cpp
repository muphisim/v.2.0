//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classLinearTetrahedra.cpp
  \brief This file contains the definition of the functions used of the classLinearTetrahedra
*/

#include "classLinearTetrahedra.h"


/*! This function gets the integration points over elements
  @param[in] order Integration order
  @param[out] npts number of integration points
  @param[out] uvw Isoparametric coordinates
  @param[out] weight all weights
*/
void classLinearTetrahedra::getIntegrationPoints(int order, int& npts, vector<vector<double> >& uvw, vector<double>& weight) const
{
    double xw[256][4];//1gp, 3dim+w
    GPsTetrahedron3DValues(order, npts, xw);//xi, eta, zeta, weights);
    uvw.clear();
    uvw.resize(npts,vector<double>(3,0));
    weight.resize(npts);
    for (int i=0; i< npts; i++)
    {
        uvw[i][0] = xw[i][0];
        uvw[i][1] = xw[i][1];
        uvw[i][2] = xw[i][2];
        weight[i] = xw[i][3];
    }
};

/*! This function gets the shape function in the isoparametric space
  @param[in] uvw Isoparametric coordinates
  @param[out] fuvw Shape functions
*/
void classLinearTetrahedra::getShapeFunctionsUVW(const vector<double>& uvw, vector<double>& fuvw) const
{
    fuvw.resize(4);
    fuvw[0] = -1.0 * uvw[0] - 1.0 * uvw[1] - 1.0 * uvw[2] + 1.0;
    fuvw[1] = 1.0 * uvw[0];
    fuvw[2] = 1.0 * uvw[1];
    fuvw[3] = 1.0 * uvw[2];
};

 /*! This function gets the shape function in the isoparametric space
  @param[in] uvw Isoparametric coordinates
  @param[out] gradfuvw Grad of shape functions
*/
void classLinearTetrahedra::getGradShapeFunctionsUVW(const vector<double>& uvw, vector<double>& gradfuvw) const
{
    gradfuvw.resize(12);
    //natural derivatives c++ style (i*ndim+j convention):
    //natural derivatives c++ style (i*ndim+j convention):
    gradfuvw[0] = -1.0 - 0.0 - 0.0 + 0.0;
    gradfuvw[3] = 1.0;
    gradfuvw[6] = 0.0;
    gradfuvw[9] = 0.0;
    gradfuvw[1] = -0.0 - 1.0 - 0.0 + 0.0;
    gradfuvw[4] = 0.0;
    gradfuvw[7] = 1.0;
    gradfuvw[10] = 0.0;
    gradfuvw[2] = -0.0 - 0.0 - 1.0 + 0.0;
    gradfuvw[5] = 0.0;
    gradfuvw[8] = 0.0;
    gradfuvw[11] = 1.0;
};


/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classLinearTetrahedra::getn0(vector<classNodes *>& nod, int ndim) {
    ERROR("Tetrahedra cannot be a subelement for the application of the pressure BCs");
    exit(-1);
}

/*! This function calculates a subelement from a given element. This function should not be called for segments
  @param[in] me Label of the element
  @param[in] S label of the edge (2D) or surface (3D) that is the subelement
  @param[out] subElement Subelement
  @param[in] nod Array with all nodes in the domain
*/
void classLinearTetrahedra::getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod,
                                          int bulkElement) {
    int numSs = myNodes.size();
    if (S > numSs) {
        ERROR("You are trying to access to a subElement that does not exist in your original element");
        exit(-1);
    }
    vector<int> mySubNodes;
    if (S == 1) 
    {
        mySubNodes.push_back(myNodes[2]);
        mySubNodes.push_back(myNodes[1]);
        mySubNodes.push_back(myNodes[0]);
    }
    else if (S==2)
    {
      mySubNodes.push_back(myNodes[3]);
      mySubNodes.push_back(myNodes[0]);
      mySubNodes.push_back(myNodes[1]);
    }
    else if (S==3)
    {
        mySubNodes.push_back(myNodes[2]);
        mySubNodes.push_back(myNodes[3]);
        mySubNodes.push_back(myNodes[1]);
    } 
    else 
    {
        mySubNodes.push_back(myNodes[2]);
        mySubNodes.push_back(myNodes[0]);
        mySubNodes.push_back(myNodes[3]);
    }

    //cout<< mySubNodes[0]<< " d "<< mySubNodes[1]<<endl;
    classElements *subElement2;
    subElement2 = new classLinearTriangles(me, "2dLinearTriangle", mySubNodes, 2, nod, bulkElement);
    subElements.push_back(subElement2);

}


/*! Constructor for classElements.
  @param[in] me Number of element
  @param[in] type Type of element
  @param[in] myNodes Nodes that build the element
  @param[in] ndim Dimension of the domain
  @param[in] nod Array with all nodes in the domain
  @param[in] bulkElement index of mother element
*/
classLinearTetrahedra::classLinearTetrahedra(int me, string type, const vector<int>& myNodes, int ndim,
                                             vector<classNodes *>& nod, int bulkElement) : classElements(me, type,
                                                                                                        myNodes, ndim,
                                                                                                        bulkElement) {
    int numNodes = this->myNodes.size();
    vector<vector<double> > nodesVec;
    for (int j = 0; j < numNodes; j++) {
        const vector<double>& currentNode = nod[myNodes[j]]->getXYZ();
        nodesVec.push_back(currentNode);
    }
    this->charLength = this->characteristicLength(ndim, nodesVec);
};

/*! This function calculates the characteristic lenght of the element. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ndim Dimension of the domain
  @param[in] nodesVec array containing postion vectors of nodes
*/
double classLinearTetrahedra::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    const vector<double>& p0 = nodesVec[0];
    const vector<double>& p1 = nodesVec[1];
    const vector<double>& p2 = nodesVec[2];
    const vector<double>& p3 = nodesVec[3];
    vector<double> tetraNodes;
    tetraNodes.resize(ndim * 4);
    for (int i = 0; i < ndim; i++) {
        tetraNodes[i * 4] = p0[i];
        tetraNodes[i * 4 + 1] = p1[i];
        tetraNodes[i * 4 + 2] = p2[i];
        tetraNodes[i * 4 + 3] = p3[i];

    }
    double radius = inscribedSphere(tetraNodes, ndim);
    return 2 * radius;
};