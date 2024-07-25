//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classQuadTetrahedra.cpp
  \brief This file contains the definition of the functions used of the classQuadTetrahedra
*/

#include "classQuadTetrahedra.h"

/*! This function gets the integration points over elements
  @param[in] order Integration order
  @param[out] npts number of integration points
  @param[out] uvw Isoparametric coordinates
  @param[out] weight all weights
*/
void classQuadTetrahedra::getIntegrationPoints(int order, int& npts, vector<vector<double> >& uvw, vector<double>& weight) const
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
void classQuadTetrahedra::getShapeFunctionsUVW(const vector<double>& uvw, vector<double>& fuvw) const
{
    fuvw.resize(10);
    fuvw[0] = 2.0 * uvw[2] * uvw[2] - 1.0 * uvw[2];
    fuvw[1] = 2.0 * uvw[1] * uvw[1] - 1.0 * uvw[1];
    fuvw[2] = 2.0 * uvw[0] * uvw[0] - 1.0 * uvw[0];
    fuvw[3] = 2.0 * uvw[0] * uvw[0] + 4.0 * uvw[0] * uvw[1] + 4.0 * uvw[0] * uvw[2] - 3.0 * uvw[0] +
             2.0 * uvw[1] * uvw[1] + 4.0 * uvw[1] * uvw[2] - 3.0 * uvw[1] + 2.0 * uvw[2] * uvw[2] -
             3.0 * uvw[2] + 1.0;
    fuvw[4] = 4.0 * uvw[1] * uvw[2];
    fuvw[5] = 4.0 * uvw[0] * uvw[1];
    fuvw[6] = 4.0 * uvw[0] * uvw[2];
    fuvw[7] = -4.0 * uvw[0] * uvw[2] - 4.0 * uvw[1] * uvw[2] - 4.0 * uvw[2] * uvw[2] + 4.0 * uvw[2];
    fuvw[8] = -4.0 * uvw[0] * uvw[1] - 4.0 * uvw[1] * uvw[1] - 4.0 * uvw[1] * uvw[2] + 4.0 * uvw[1];
    fuvw[9] = -4.0 * uvw[0] * uvw[0] - 4.0 * uvw[0] * uvw[1] - 4.0 * uvw[0] * uvw[2] + 4.0 * uvw[0];
};

 /*! This function gets the shape function in the isoparametric space
  @param[in] uvw Isoparametric coordinates
  @param[out] gradfuvw Grad of shape functions
*/
void classQuadTetrahedra::getGradShapeFunctionsUVW(const vector<double>& uvw, vector<double>& gradfuvw) const
{
    gradfuvw.resize(30);
    //natural derivatives c++ style (i*ndim+j convention):
    gradfuvw[0] = 0.0 - 0.0;
    gradfuvw[3] = 0.0 - 0.0;
    gradfuvw[6] = 4.0 * uvw[0] - 1.0;
    gradfuvw[9] = 4.0 * uvw[0] + 4.0 * uvw[1] + 4.0 * uvw[2] - 3.0 + 0.0 + 0.0 - 0.0 + 0.0 - 0.0 + 0.0;
    gradfuvw[12] = 0.0;
    gradfuvw[15] = 4.0 * uvw[1];
    gradfuvw[18] = 4.0 * uvw[2];
    gradfuvw[21] = -4.0 * uvw[2] - 0.0 - 0.0 + 0.0;
    gradfuvw[24] = -4.0 * uvw[1] - 0.0 - 0.0 + 0.0;
    gradfuvw[27] = -8.0 * uvw[0] - 4.0 * uvw[1] - 4.0 * uvw[2] + 4.0;
    gradfuvw[1] = 0.0 - 0.0;
    gradfuvw[4] = 4.0 * uvw[1] - 1.0;
    gradfuvw[7] = 0.0 - 0.0;
    gradfuvw[10] = 0.0 + 4.0 * uvw[0] + 0.0 - 0.0 + 4.0 * uvw[1] + 4.0 * uvw[2] - 3.0 + 0.0 - 0.0 + 0.0;
    gradfuvw[13] = 4.0 * uvw[2];
    gradfuvw[16] = 4.0 * uvw[0];
    gradfuvw[19] = 0.0;
    gradfuvw[22] = -0.0 - 4.0 * uvw[2] - 0.0 + 0.0;
    gradfuvw[25] = -4.0 * uvw[0] - 8.0 * uvw[1] - 4.0 * uvw[2] + 4.0;
    gradfuvw[28] = -0.0 - 4.0 * uvw[0] - 0.0 + 0.0;
    gradfuvw[2] = 4.0 * uvw[2] - 1.0;
    gradfuvw[5] = 0.0 - 0.0;
    gradfuvw[8] = 0.0 - 0.0;
    gradfuvw[11] = 0.0 + 0.0 + 4.0 * uvw[0] - 0.0 + 0.0 + 4.0 * uvw[1] - 0.0 + 4.0 * uvw[2] - 3.0 + 0.0;
    gradfuvw[14] = 4.0 * uvw[1];
    gradfuvw[17] = 0.0;
    gradfuvw[20] = 4.0 * uvw[0];
    gradfuvw[23] = -4.0 * uvw[0] - 4.0 * uvw[1] - 8.0 * uvw[2] + 4.0;
    gradfuvw[26] = -0.0 - 0.0 - 4.0 * uvw[1] + 0.0;
    gradfuvw[29] = -0.0 - 0.0 - 4.0 * uvw[0] + 0.0;
};

/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classQuadTetrahedra::getn0(vector<classNodes *>& nod, int ndim) {
    ERROR("Tetrahedra cannot be a subelement for the application of the pressure BCs");
    exit(-1);
}

/*! This function calculates a subelement from a given element. This function should not be called for segments
  @param[in] me Label of the element
  @param[in] S label of the edge (2D) or surface (3D) that is the subelement
  @param[out] subElement Subelement
  @param[in] nod Array with all nodes in the domain
*/
void classQuadTetrahedra::getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod,
                                        int bulkElement) {
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
            mySubNodes.push_back(myNodes[4]);
            mySubNodes.push_back(myNodes[5]);
            mySubNodes.push_back(myNodes[6]);
            break;
        }
        case 2: {
            mySubNodes.push_back(myNodes[0]);
            mySubNodes.push_back(myNodes[3]);
            mySubNodes.push_back(myNodes[1]);
            mySubNodes.push_back(myNodes[7]);
            mySubNodes.push_back(myNodes[8]);
            mySubNodes.push_back(myNodes[4]);
            break;
        }
        case 3: {
            mySubNodes.push_back(myNodes[1]);
            mySubNodes.push_back(myNodes[3]);
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[8]);
            mySubNodes.push_back(myNodes[9]);
            mySubNodes.push_back(myNodes[5]);
            break;
        }
        case 4: {
            mySubNodes.push_back(myNodes[0]);
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[3]);
            mySubNodes.push_back(myNodes[6]);
            mySubNodes.push_back(myNodes[9]);
            mySubNodes.push_back(myNodes[7]);
            break;
        }
    }
    //for a given S_a:
    //quad tri cornernodes (counterclockwise) a0, a1, a2
    //quad tri midnodes (counterclockwise, first node right adjecent to a0) a3, a4, a5
    //so subelements S_a1 = {a0, a3, a5}, S_a2={a3, a1, a4}, S_a3={a5, a4, a2}, S_a4={a5, a3, a4}
    classElements *subElement2;
    subElement2 = new classQuadTriangles(me, "2dQuadTriangle", mySubNodes, 2, nod, bulkElement);
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
classQuadTetrahedra::classQuadTetrahedra(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
                                         int bulkElement) : classElements(me, type, myNodes, ndim, bulkElement) {
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
double classQuadTetrahedra::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    double minDist = 1E25;
    //get all nodes in an array

    //divide tetrahedron into subelements, calculate radius of the inscribed sphere for each subelement and pick the smallest
    vector<vector<int> > subelements = {{6, 5, 2, 9},
                                        {4, 8, 1, 5},
                                        {0, 4, 6, 7},
                                        {3, 5, 6, 9},
                                        {8, 3, 4, 5},
                                        {7, 3, 4, 6},
                                        {0, 4, 6, 7}};
    double radii[7];
    for (int j = 0; j < 7; j++)//tetrahedron subdivides into 7 subtetrahedras with cornernodes and midnodes
    {
        //each subelement can again be subdivided into 4 smaller subsubelements (all subelements are tetrahedras)
        //the subdivision of each tetrahedral element is such that each of the subsubelements has a base corresponding to one of
        //the triangular surfaces of the subelement, and the last node is center of the inscribed sphere. The sum of volumes of
        //subsubelements is equal to the volume of the subelement, and since v=1/3*a*h, where h of a subsubelement is the radius of the
        //sphere, we know that total volume of the subelement = 1/3*(s1+s2+s3+s4)*r
        const vector<double>& p0 = nodesVec[subelements[j][0]];
        const vector<double>& p1 = nodesVec[subelements[j][1]];
        const vector<double>& p2 = nodesVec[subelements[j][2]];
        const vector<double>& p3 = nodesVec[subelements[j][3]];//vertecies of current subelement in the loop, for convenience
        vector<double> tetraNodes;
        tetraNodes.resize(ndim * 4);
        for (int i = 0; i < ndim; i++) {
            tetraNodes[i * 4] = p0[i];
            tetraNodes[i * 4 + 1] = p1[i];
            tetraNodes[i * 4 + 2] = p2[i];
            tetraNodes[i * 4 + 3] = p3[i];

        }
        radii[j] = inscribedSphere(tetraNodes, ndim);
    }
    double rmin = radii[0];
    for (int i = 1; i < 7; i++) {
        if (rmin > radii[i])
            rmin = radii[i];
    }
    if (rmin < minDist) {
        minDist = 2 * rmin;//diameter!!! radius was way slow
    }
    return minDist;
};