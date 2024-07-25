//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classLinearHexahedra.cpp
  \brief This file contains the definition of the functions used of the classLinearHexahedra
*/

#include "classLinearHexahedra.h"

/*! This function gets the integration points over elements
  @param[in] order Integration order
  @param[out] npts number of integration points
  @param[out] uvw Isoparametric coordinates
  @param[out] weight all weights
*/
void classLinearHexahedra::getIntegrationPoints(int order, int& npts, vector<vector<double> >& uvw, vector<double>& weight) const
{
    double xw[256][4];//1gp, 3dim+w
    GPsHexahedron3DValues(order, npts, xw);//xi, eta, zeta, weights);
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
void classLinearHexahedra::getShapeFunctionsUVW(const vector<double>& uvw, vector<double>& fuvw) const
{
    fuvw.resize(8);
    fuvw[0] = 0.125 * (1 - uvw[0]) * (1 - uvw[1]) * (1 - uvw[2]);
    fuvw[1] = 0.125 * (1 + uvw[0]) * (1 - uvw[1]) * (1 - uvw[2]);
    fuvw[2] = 0.125 * (1 + uvw[0]) * (1 + uvw[1]) * (1 - uvw[2]);
    fuvw[3] = 0.125 * (1 - uvw[0]) * (1 + uvw[1]) * (1 - uvw[2]);
    fuvw[4] = 0.125 * (1 - uvw[0]) * (1 - uvw[1]) * (1 + uvw[2]);
    fuvw[5] = 0.125 * (1 + uvw[0]) * (1 - uvw[1]) * (1 + uvw[2]);
    fuvw[6] = 0.125 * (1 + uvw[0]) * (1 + uvw[1]) * (1 + uvw[2]);
    fuvw[7] = 0.125 * (1 - uvw[0]) * (1 + uvw[1]) * (1 + uvw[2]);
};

 /*! This function gets the shape function in the isoparametric space
  @param[in] uvw Isoparametric coordinates
  @param[out] gradfuvw Grad of shape functions
*/
void classLinearHexahedra::getGradShapeFunctionsUVW(const vector<double>& uvw, vector<double>& gradfuvw) const
{
    gradfuvw.resize(24);
    //natural derivatives c++ style (i*ndim+j convention):
    gradfuvw[0] = -0.125 * (1 - uvw[1]) * (1 - uvw[2]);
    gradfuvw[3] = 0.125 * (1 - uvw[1]) * (1 - uvw[2]);
    gradfuvw[6] = 0.125 * (1 + uvw[1]) * (1 - uvw[2]);
    gradfuvw[9] = -0.125 * (1 + uvw[1]) * (1 - uvw[2]);
    gradfuvw[12] = -0.125 * (1 - uvw[1]) * (1 + uvw[2]);
    gradfuvw[15] = 0.125 * (1 - uvw[1]) * (1 + uvw[2]);
    gradfuvw[18] = 0.125 * (1 + uvw[1]) * (1 + uvw[2]);
    gradfuvw[21] = -0.125 * (1 + uvw[1]) * (1 + uvw[2]);
    gradfuvw[1] = -0.125 * (1 - uvw[0]) * (1 - uvw[2]);
    gradfuvw[4] = -0.125 * (1 + uvw[0]) * (1 - uvw[2]);
    gradfuvw[7] = 0.125 * (1 + uvw[0]) * (1 - uvw[2]);
    gradfuvw[10] = 0.125 * (1 - uvw[0]) * (1 - uvw[2]);
    gradfuvw[13] = -0.125 * (1 - uvw[0]) * (1 + uvw[2]);
    gradfuvw[16] = -0.125 * (1 + uvw[0]) * (1 + uvw[2]);
    gradfuvw[19] = 0.125 * (1 + uvw[0]) * (1 + uvw[2]);
    gradfuvw[22] = 0.125 * (1 - uvw[0]) * (1 + uvw[2]);
    gradfuvw[2] = -0.125 * (1 - uvw[0]) * (1 - uvw[1]);
    gradfuvw[5] = -0.125 * (1 + uvw[0]) * (1 - uvw[1]);
    gradfuvw[8] = -0.125 * (1 + uvw[0]) * (1 + uvw[1]);
    gradfuvw[11] = -0.125 * (1 - uvw[0]) * (1 + uvw[1]);
    gradfuvw[14] = 0.125 * (1 - uvw[0]) * (1 - uvw[1]);
    gradfuvw[17] = 0.125 * (1 + uvw[0]) * (1 - uvw[1]);
    gradfuvw[20] = 0.125 * (1 + uvw[0]) * (1 + uvw[1]);
    gradfuvw[23] = 0.125 * (1 - uvw[0]) * (1 + uvw[1]);
};


/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classLinearHexahedra::getn0(vector<classNodes *>& nod, int ndim){
    ERROR("Hexahedra cannot be a subelement for the application of the pressure BCs");
    exit(-1);
}

/*! This function calculates a subelement from a given element. This function should not be called for segments
  @param[in] me Label of the element
  @param[in] S label of the edge (2D) or surface (3D) that is the subelement
  @param[out] subElement Subelement
  @param[in] nod Array with all nodes in the domain
*/
void classLinearHexahedra::getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod, int bulkElement) {
    int numSs = myNodes.size();
    if (S > numSs) {
        ERROR("You are trying to access to a subElement that does not exist in your original element");
        exit(-1);
    }
    vector<int> mySubNodes;
    switch (S) {
        case 1: {
            mySubNodes.push_back(myNodes[0]);
            mySubNodes.push_back(myNodes[1]);
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[3]);
            break;
        }
        case 2: {
            mySubNodes.push_back(myNodes[4]);
            mySubNodes.push_back(myNodes[7]);
            mySubNodes.push_back(myNodes[6]);
            mySubNodes.push_back(myNodes[5]);
            break;
        }
        case 3: {
            mySubNodes.push_back(myNodes[0]);
            mySubNodes.push_back(myNodes[4]);
            mySubNodes.push_back(myNodes[5]);
            mySubNodes.push_back(myNodes[1]);
            break;
        }
        case 4: {
            mySubNodes.push_back(myNodes[1]);
            mySubNodes.push_back(myNodes[5]);
            mySubNodes.push_back(myNodes[6]);
            mySubNodes.push_back(myNodes[2]);
            break;
        }
        case 5: {
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[6]);
            mySubNodes.push_back(myNodes[7]);
            mySubNodes.push_back(myNodes[3]);
            break;
        }
        case 6: {
            mySubNodes.push_back(myNodes[3]);
            mySubNodes.push_back(myNodes[7]);
            mySubNodes.push_back(myNodes[4]);
            mySubNodes.push_back(myNodes[0]);
            break;
        }
    }

    //cout<< mySubNodes[0]<< " d "<< mySubNodes[1]<<endl;
    classElements *subElement2;
    subElement2 = new classLinearSquare(me, "2dLinearSquare", mySubNodes, 2, nod, bulkElement);
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
classLinearHexahedra::classLinearHexahedra(int me, string type,  const vector<int>& myNodes, int ndim,
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
double classLinearHexahedra::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    const vector<double>& p0 = nodesVec[0];
    const vector<double>& p1 = nodesVec[1];
    const vector<double>& p2 = nodesVec[2];
    const vector<double>& p3 = nodesVec[3];
    const vector<double>& p4 = nodesVec[4];
    const vector<double>& p5 = nodesVec[5];
    const vector<double>& p6 = nodesVec[6];
    const vector<double>& p7 = nodesVec[7];
    vector<double> hexaNodes;
    hexaNodes.resize(ndim * 8);
    for (int i = 0; i < ndim; i++) {
        hexaNodes[i * 8] = p0[i];
        hexaNodes[i * 8 + 1] = p1[i];
        hexaNodes[i * 8 + 2] = p2[i];
        hexaNodes[i * 8 + 3] = p3[i];
        hexaNodes[i * 8 + 4] = p4[i];
        hexaNodes[i * 8 + 5] = p5[i];
        hexaNodes[i * 8 + 6] = p6[i];
        hexaNodes[i * 8 + 7] = p7[i];
    }
    double radius = inscribedSphereHex(hexaNodes, ndim);
    return 2 * radius;
};
