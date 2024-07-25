//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file FSIAdapter.h
  \brief This file contains all functions related the FSI adapter using preCICE
*/


#ifndef _FSIAdapter_H_
#define _FSIAdapter_H_

#ifdef FSI
#include "SolverInterface.hpp"
#include <string>
#include <vector>
#include "classNodes.h"
#include "classElements.h"

class FSIAdapter
{
   protected:
    static FSIAdapter* _instance;
  public:
    static FSIAdapter* getInstance();
    static void clear();
    
  protected:
    precice::SolverInterface* _precice;
    // data
    std::string _meshName, _solverName, _preciceConfigFile;
    std::vector<std::string>_dataWriteName, _dataReadName, _forceInterfaceName;
    std::map<std::string, std::vector<double> > _forceInterfaceParameters;
    std::vector<int> _vertexIds; // list of boundary nodes
    int _vertexSize; // number of boundary nodes
    std::map<std::string, std::vector<double> > _readData;
    std::map<std::string, std::vector<double> > _readDataPrev;

    // data checkpoint of the solver
    double _dtIterFSI, _timeIterFSI;
    
    // interface data
    std::vector<int> _boundaryNodes;
    std::vector<classElements *> _boundaryElements;
    std::vector<std::vector <double> > _boundaryCoordinates;
    std::map<int, int> _nodeMap; // boundary node is map to _vertexIds
    
  public:
    FSIAdapter(const std::string config_file_name="../precice-config.xml",
               const std::string solver_name = "Solid", 
               const std::string mesh_name = "Solid-Mesh",
               const std::string data_write_name = "Displacement",
               const std::string data_read_name ="Force");
    ~FSIAdapter();
    void loadFSISettingFromFile(std::string fileName);
    
    void setBoundaryNodes(const std::vector<int>& boundary_nodes);
    void setBoundaryNodeCoordinates(const std::vector<vector <double> >& boundary_coordinates);
    void setBoundaryElements(const std::vector<classElements *>& boundaryElements);
    
    double initialize();
    void updateDisplacementInterface(double dt, const std::vector<double>& disp, int numNodes);
    void readBlockVectorData();
    void updateForceInterface(vector<classNodes *> &nod, vector<classElements *> &elem, int &ndim, double timeRun, double dt, 
                              const std::vector<double>& disp, std::vector<double>& force, int numNodes, double fact,
                              bool tangentEstimation, map<pair<int,int>,double>& stiffBC);
    double advance(double dt);
    bool isCouplingOngoing();
    void finalize();
    
    const std::string& actionReadIterationCheckpoint();
    const std::string& actionWriteIterationCheckpoint();
    
    bool isActionRequired(const std::string& action);
    void markActionFulfilled(const std::string& action);
    
                      
    void makeACheckPoint(double timeRun, double dt);
    void restoreACheckPoint(double& timeRun, double& dt);
};

#endif // FSI

#endif //_FSIAdapter_H_
