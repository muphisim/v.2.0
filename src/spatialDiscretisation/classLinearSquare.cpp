//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classLinearSquare.cpp
  \brief This file contains the definition of the functions for classLinearSquare
*/

#include "classLinearSquare.h"

/*! This function gets the integration points over elements
  @param[in] order Integration order
  @param[out] npts number of integration points
  @param[out] uvw Isoparametric coordinates
  @param[out] weight all weights
*/
void classLinearSquare::getIntegrationPoints(int order, int& npts, vector<vector<double> >& uvw, vector<double>& weight) const
{
    double xw[256][3];//1gp, 2dim+w
    GPsSquare2DValues(order, npts, xw);//xi, eta, zeta, weights);
    uvw.clear();
    uvw.resize(npts,vector<double>(2,0));
    weight.resize(npts);
    for (int i=0; i< npts; i++)
    {
        uvw[i][0] = xw[i][0];
        uvw[i][1] = xw[i][1];
        weight[i] = xw[i][2];
    }
};

/*! This function gets the shape function in the isoparametric space
  @param[in] uvw Isoparametric coordinates
  @param[out] fuvw Shape functions
*/
void classLinearSquare::getShapeFunctionsUVW(const vector<double>& uvw, vector<double>& fuvw) const
{
    fuvw.resize(4);
    fuvw[0] = 0.25 * (1 - uvw[0]) * (1 - uvw[1]);
    fuvw[1] = 0.25 * (1 + uvw[0]) * (1 - uvw[1]);
    fuvw[2] = 0.25 * (1 + uvw[0]) * (1 + uvw[1]);
    fuvw[3] = 0.25 * (1 - uvw[0]) * (1 + uvw[1]);

};

 /*! This function gets the shape function in the isoparametric space
  @param[in] uvw Isoparametric coordinates
  @param[out] gradfuvw Grad of shape functions
*/
void classLinearSquare::getGradShapeFunctionsUVW(const vector<double>& uvw, vector<double>& gradfuvw) const
{
    gradfuvw.resize(8);
    //natural derivatives c++ style (i*ndim+j convention):
    gradfuvw[0] = -0.25 * (1 - uvw[1]);
    gradfuvw[2] = 0.25 * (1 - uvw[1]);
    gradfuvw[4] = 0.25 * (1 + uvw[1]);
    gradfuvw[6] = -0.25 * (1 + uvw[1]);
    gradfuvw[1] = -0.25 * (1 - uvw[0]) ;
    gradfuvw[3] = -0.25 * (1 + uvw[0]) ;
    gradfuvw[5] = 0.25 * (1 + uvw[0]) ;
    gradfuvw[7] = 0.25 * (1 - uvw[0]) ;
};


/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classLinearSquare::getn0(vector<classNodes *>& nod, int ndim) {
    const vector<double>& xyz1 = nod[myNodes[0]]->getXYZ();
    const vector<double>& xyz2 = nod[myNodes[1]]->getXYZ();
    const vector<double>& xyz3 = nod[myNodes[2]]->getXYZ();

    vector<double> n0(ndim, 0);

    for (int i = 0; i < ndim; i++) {
        int coord1, coord2;
        if (i == 0) //nx
        {
            coord1 = 1;
            coord2 = 2;
        } else if (i == 1) //ny
        {
            coord1 = 2;
            coord2 = 0;
        } else {
            coord1 = 0;
            coord2 = 1;
        }
        n0[i] = (xyz2[coord1] - xyz1[coord1]) * (xyz3[coord2] - xyz1[coord2]) -
                (xyz3[coord1] - xyz1[coord1]) * (xyz2[coord2] - xyz1[coord2]);
    }
    double no = norm2(n0, ndim);
    for (int i = 0; i < ndim; i++) {
        n0[i] = n0[i] / no;
    }
    return n0;

}

/*! This function calculates a subelement from a given element. This function should not be called for segments
  @param[in] me Label of the element
  @param[in] S label of the edge (2D) or surface (3D) that is the subelement
  @param[out] subElement Subelement
  @param[in] nod Array with all nodes in the domain
*/
void classLinearSquare::getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod,
                                         int bulkElement) {
    // I know I am a Square
    int numSs = myNodes.size();
    if (S > numSs) {
        ERROR("You are trying to access to a subElement that does not exist in your original element");
        exit(-1);
    }
    vector<int> mySubNodes;
    if (S != 4) {
        mySubNodes.push_back(myNodes[S - 1]);
        mySubNodes.push_back(myNodes[S]);
    } else {
        mySubNodes.push_back(myNodes[S - 1]);
        mySubNodes.push_back(myNodes[0]);
    }
    classElements *subElement2;
    subElement2 = new classLinearLine(me, "XXX1dLinear", mySubNodes, 1, nod, bulkElement);//XXX call 1D element
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
classLinearSquare::classLinearSquare(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
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
double classLinearSquare::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    vector<double> squareNodes;
    squareNodes.resize(ndim * 4, 0);

    const vector<double>& p0 = nodesVec[0];
    const vector<double>& p1 = nodesVec[1];
    const vector<double>& p2 = nodesVec[2];
    const vector<double>& p3 = nodesVec[3];

    for (int i = 0; i < ndim; i++) {
        squareNodes[i * 4] = p0[i];
        squareNodes[i * 4 + 1] = p1[i];
        squareNodes[i * 4 + 2] = p2[i];
        squareNodes[i * 4 + 3] = p3[i];
    }
    double radius = inscribedCircleQuad(squareNodes, ndim);
    return 2 * radius;
};