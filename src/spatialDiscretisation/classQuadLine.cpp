//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classQuadLine.cpp
  \brief This file contains the definition of the functions for classQuadLines
*/

#include "classQuadLine.h"

/*! This function gets the integration points over elements
  @param[in] order Integration order
  @param[out] npts number of integration points
  @param[out] uvw Isoparametric coordinates
  @param[out] weight all weights
*/
void classQuadLine::getIntegrationPoints(int order, int& npts, vector<vector<double> >& uvw, vector<double>& weight) const
{
    double xw[256][2];//nr of gps, nr of abcissa+1weight/gp
    GPsLine1DValues(order, npts, xw); 
    uvw.clear();
    uvw.resize(npts,vector<double>(1,0));
    weight.resize(npts);
    for (int i=0; i< npts; i++)
    {
        uvw[i][0] = xw[i][0];
        weight[i] = xw[i][1];
    };
};

/*! This function gets the shape function in the isoparametric space
  @param[in] uvw Isoparametric coordinates
  @param[out] fuvw Shape functions
*/
void classQuadLine::getShapeFunctionsUVW(const vector<double>& uvw, vector<double>& fuvw) const
{
    fuvw.resize(2);
    //shape functions c++ style:
    fuvw[0] = 0.5 * uvw[0] * uvw[0] - 0.5 * uvw[0];
    fuvw[1] = 0.5 * uvw[0] * uvw[0] + 0.5 * uvw[0];
    fuvw[2] = -1.0 * uvw[0] * uvw[0] + 1.0;
};

 /*! This function gets the shape function in the isoparametric space
  @param[in] uvw Isoparametric coordinates
  @param[out] gradfuvw Grad of shape functions
*/
void classQuadLine::getGradShapeFunctionsUVW(const vector<double>& uvw, vector<double>& gradfuvw) const
{
    gradfuvw.resize(2);
    gradfuvw[0] = 1.0 * uvw[0] - 0.5;
    gradfuvw[1] = 1.0 * uvw[0] + 0.5;
    gradfuvw[2] = -2.0 * uvw[0] + 0.0;
};


/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classQuadLine::getn0(vector<classNodes *>& nod, int ndim) { //should work assuming midnode has index 2
    vector<double> n0(ndim, 0);
    const vector<double>& xyz1 = nod[myNodes[0]]->getXYZ();
    const vector<double>& xyz2 = nod[myNodes[1]]->getXYZ();

    n0[0] = (xyz2[1] - xyz1[1]);
    n0[1] = (xyz2[0] - xyz1[0]) * -1;

    double no = norm2(n0, ndim);
    for (int i = 0; i < ndim; i++) {
        n0[i] = n0[i] / no;
    }
    return n0;
}

/*! Constructor for classElements.
  @param[in] me Number of element
  @param[in] type Type of element
  @param[in] myNodes Nodes that build the element
  @param[in] ndim Dimension of the domain
  @param[in] nod Array with all nodes in the domain
  @param[in] bulkElement index of mother element
*/
classQuadLine::classQuadLine(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
                             int bulkElement) : classElements(me, type, myNodes, ndim, bulkElement) {
    int numNodes = this->myNodes.size();
    vector<vector<double>> nodesVec;
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
double classQuadLine::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    const vector<double>& xyz1 = nodesVec[0];
    const vector<double>& xyz2 = nodesVec[1];
    const vector<double>& xyz3 = nodesVec[2];
    double dist1 = calcDist(xyz1, xyz2, ndim);
    double dist2 = calcDist(xyz1, xyz3, ndim);
    double dist3 = calcDist(xyz2, xyz3, ndim);
    double dist = std::min(dist1, dist2);
    dist = std::min(dist, dist3);
    return dist;
};
