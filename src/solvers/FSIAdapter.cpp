//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file FSIAdapter.cpp
  \brief This file contains all functions related the FSI adapter using preCICE
*/


#include "FSIAdapter.h"
#include "classGPs.h"
#include "commandLine.h"
#include "boundaryConditions.h"

#ifdef FSI

/*! \brief constructor
  @param[in] config_file_name path to the configure file
  @param[in] solver_name name of the solver
  @param[in] mesh_name name of the mesh
  @param[in] data_write_name name of data for writing
  @param[in] data_read_name name of data for reading
*/
FSIAdapter::FSIAdapter(const std::string config_file_name,
               const std::string solver_name, 
               const std::string mesh_name,
               const std::string data_write_name,
               const std::string data_read_name):
               _preciceConfigFile(config_file_name),
               _solverName(solver_name),_meshName(mesh_name), 
               _dataWriteName(1,data_write_name), 
               _dataReadName(1,data_read_name), 
               _forceInterfaceName(1,"Force"),
               _vertexSize(0), _dtIterFSI(0.), _timeIterFSI(0.),
               _precice(NULL)
{
  
};

FSIAdapter* FSIAdapter::_instance = NULL;
FSIAdapter* FSIAdapter::getInstance()
{
  if (FSIAdapter::_instance == NULL)
  {
    FSIAdapter::_instance = new FSIAdapter();
  }
  return FSIAdapter::_instance;
}
void FSIAdapter::clear()
{
  if (FSIAdapter::_instance != NULL)
  {
    delete FSIAdapter::_instance;
    FSIAdapter::_instance = NULL;
  }
}

FSIAdapter::~FSIAdapter()
{
  if (_precice !=NULL)
  {
    delete _precice;
    _precice = NULL;
  }
}

/*! \brief Load FSI setting from file
  @param[in] fileName filename in current directory
*/

void FSIAdapter::loadFSISettingFromFile(std::string fileName)
{

  printf("read FSI setting file: %s \n",fileName.c_str());
  FILE *pFile = fopen(fileName.c_str(), "r");
  if(pFile == NULL) 
  {
    printf("File %s could not be opened\n", fileName.c_str());
    return;
  }
  char what[256];
  while(!feof(pFile)) 
  {
    if(fscanf(pFile, "%s", what) != 1) 
    {
      fclose(pFile);
      return;
    }
    if(what[0] == '#') 
    {
      char buffer[1024];
      if(fgets(buffer, sizeof(buffer), pFile) == nullptr)
      {
        printf("intructions cannot be read!!!\n");
      }
    }
    else if (!strcmp(what, "*PRECICE_CONFIG_FILE")) 
    {
      // load geometry
      char line[256];
      if(fscanf(pFile, "%s", line) != 1) 
      {
        fclose(pFile);
        return;
      }
      printf("precice-config-file: %s\n",line);
      _preciceConfigFile  = line;
    }
    else if(!strcmp(what, "*SOLVER_NAME")) 
    {
      // load geometry
      char line[256];
      if(fscanf(pFile, "%s", line) != 1) 
      {
        fclose(pFile);
        return;
      }
      printf("solver name: %s\n",line);
      _solverName  = line;
      
    }
    else if(!strcmp(what, "*MESH_NAME")) 
    {
      char line[256];
      if(fscanf(pFile, "%s", line) != 1) 
      {
        fclose(pFile);
        return;
      }
      printf("mesh name: %s\n",line);
      _meshName = line;
    }
    else if(!strcmp(what, "*DATA_WRITE")) 
    {
      int numData;
      if(fscanf(pFile, "%d", &numData) != 1) 
      {
        fclose(pFile);
        return;
      }
      _dataWriteName.clear();
      for (int i=0; i< numData; i++)
      {
        char line[256];
        if(fscanf(pFile, "%s", line) != 1) 
        {
          fclose(pFile);
          return;
        }
        _dataWriteName.push_back(line);
      }
      for (int i=0; i< _dataWriteName.size(); i++)
      {
        printf("data write to other solver %d %s \n",i,_dataWriteName[i].c_str());
      }
    }
    else if(!strcmp(what, "*DATA_READ")) 
    {
      int numData;
      if(fscanf(pFile, "%d", &numData) != 1) 
      {
        fclose(pFile);
        return;
      }
      _dataReadName.clear();
      for (int i=0; i< numData; i++)
      {
        char line[256];
        if(fscanf(pFile, "%s", line) != 1) 
        {
          fclose(pFile);
          return;
        }
        _dataReadName.push_back(line);
      }
      for (int i=0; i< _dataReadName.size(); i++)
      {
        printf("data read from other solver: %d %s \n",i,_dataReadName[i].c_str());
      }
    }
    else if(!strcmp(what, "*FORCE_INTERFACE")) 
    {
      int numData;
      if(fscanf(pFile, "%d", &numData) != 1) 
      {
        fclose(pFile);
        return;
      }
      _forceInterfaceName.clear();
      for (int i=0; i< numData; i++)
      {
        char line[256];
        int numVars;
        if(fscanf(pFile, "%s %d", &line,&numVars) != 2) 
        {
          fclose(pFile);
          return;
        }
        std::string rhsName= line;
        _forceInterfaceName.push_back(rhsName);
        if (numVars > 0)
        {
          _forceInterfaceParameters[rhsName].resize(numVars);
          for (int j=0; j< numVars; j++)
          {
            if(fscanf(pFile, "%lf", &_forceInterfaceParameters[rhsName][j]) != 1) 
            {
              fclose(pFile);
              return;
            }
          }
        }
      }
      for (int i=0; i< _forceInterfaceName.size(); i++)
      {
        printf("force interface: %d %s \n",i,_forceInterfaceName[i].c_str());
        if (_forceInterfaceParameters.find(_forceInterfaceName[i]) != _forceInterfaceParameters.end())
        {
          printVector(_forceInterfaceParameters[_forceInterfaceName[i]],"parameters for "+_forceInterfaceName[i]+":");
        }
        else
        {
          printf("no parameter\n");
        }
      }
    }
    else
    {
      printf("Input %s is not valid\n", what);
    }
  }
  fclose(pFile);
};

void FSIAdapter::setBoundaryNodes(const vector<int>& boundary_nodes)
{
  _boundaryNodes = boundary_nodes;
};
void FSIAdapter::setBoundaryNodeCoordinates(const vector<vector <double> >& boundary_coordinates)
{
  _boundaryCoordinates = boundary_coordinates;
};

void FSIAdapter::setBoundaryElements(const std::vector<classElements *>& boundaryElements)
{
  _boundaryElements = boundaryElements;
};

/*! \brief initialize function
*/
double FSIAdapter::initialize()
{
  
  if (_precice!=NULL)
    delete _precice;
  _precice = new precice::SolverInterface(_solverName,_preciceConfigFile,GeneralOptions::commRank,GeneralOptions::commSize);
  
  _vertexSize = _boundaryNodes.size();
  int dim = _precice->getDimensions();
  std::vector<double> boundary_coords(_vertexSize*dim); // coords of coupling vertices
  for (int j=0; j <_vertexSize; j++)
  {
    std::vector<double>& xyz = _boundaryCoordinates[j];
    if (dim == 2) 
    { //2D
      boundary_coords[dim*j] = xyz[0];
      boundary_coords[dim*j + 1] = xyz[1];
    } else if (dim == 3) 
    { //3D
      boundary_coords[dim*j] = xyz[0];
      boundary_coords[dim*j + 1] = xyz[1];
      boundary_coords[dim*j + 2] = xyz[2];
    }
  };
  
  int mesh_id  = _precice->getMeshID(_meshName);
  _vertexIds.resize(_vertexSize);
  _precice->setMeshVertices(mesh_id,_vertexSize, boundary_coords.data(), _vertexIds.data());
  
  _nodeMap.clear();
  for (int in=0; in < _vertexSize; in++)
  {
    _nodeMap[_boundaryNodes[in]] = in;
  }
  /*
  printVector(_vertexIds,"_vertexIds");
  printVector(_boundaryNodes,"_boundaryNodes");
  for (int ie=0; ie< _boundaryElements.size(); ie++)
  {
    std::string text = "element:" + std::to_string(_boundaryElements[ie]->getMe());
    printVector(_boundaryElements[ie]->getMyNodes(), text);
  } 
   */
   
  
  for (int i=0; i< _dataReadName.size(); i++)
  {
    _readData[_dataReadName[i]].resize(_vertexSize*dim,0.);
    _readDataPrev[_dataReadName[i]].resize(_vertexSize*dim,0.);
  }
  return _precice->initialize();
};

/*! \brief update the displacement at interface and send to other solver
  @param[in] dt time step
  @param[in] disp vector of all degrees of freedom of unknown field
  @param[in] numNodes number of nodes of the domain
*/
void FSIAdapter::updateDisplacementInterface(double dt, const std::vector<double>& disp, int numNodes)
{
  if (_precice->isWriteDataRequired(dt))
  {
    int dim = _precice->getDimensions();
    for (int iw =0; iw < _dataWriteName.size(); iw++)
    {
      if (_dataWriteName[iw] == "Displacement")
      {
        std::vector<double> displacements(_vertexSize*dim,0.);
        // displacements : vector containing the value of the displacement on the boundary nodes
        for (int j=0; j <_vertexSize; j++)
        {
          int i = _boundaryNodes[j];
          if (dim == 2) 
          { //2D
            displacements[dim*j] = disp[i];
            displacements[dim*j + 1] = disp[i + numNodes];
          } 
          else if (dim == 3) 
          { //3D
            displacements[dim*j] = disp[i];
            displacements[dim*j + 1] = disp[i + numNodes];
            displacements[dim*j + 2] = disp[i + 2*numNodes];
          }
        }
        int mesh_id  = _precice->getMeshID(_meshName);
        int write_id = _precice->getDataID(_dataWriteName[iw], mesh_id);
        _precice->writeBlockVectorData(write_id, _vertexSize, _vertexIds.data(), displacements.data());
      }
      else
      {
        ERROR("write data %s has not been implemented",_dataWriteName[iw].c_str());
      }
    }
  };
}

/*! \brief read data from the other solver
*/
void FSIAdapter::readBlockVectorData()
{
  if (_precice->isReadDataAvailable()) 
  {
    int mesh_id  = _precice->getMeshID(_meshName);
    for (int i=0; i< _dataReadName.size(); i++)
    {
      int read_id  = _precice->getDataID(_dataReadName[i], mesh_id);
      _readDataPrev[_dataReadName[i]] = _readData[_dataReadName[i]];
      _precice->readBlockVectorData(read_id, _vertexSize, _vertexIds.data(), _readData[_dataReadName[i]].data());
    }
  };
}

/*! \brief update the external force vector with the data from other solver
  @param[in] force external force vector
  @param[in] numNodes number of nodes of the domain
*/
void FSIAdapter::updateForceInterface(vector<classNodes *> &nod, vector<classElements *> &elem, int &ndim, double timeRun, double dt, 
                                    const std::vector<double>& disp, std::vector<double>& force, int numNodes, double fact,
                                    bool tangentEstimation, map<pair<int,int>,double>& stiffBC)
{  
  int dim = _precice->getDimensions();
  for (int dd=0; dd< _forceInterfaceName.size(); dd++)
  { 
    if (_forceInterfaceName[dd]=="Force")
    {
      if (_readData.find("Force") == _readData.end())
      {
        ERROR("Force is not present in read data");
        exit(-1);
      }
      std::vector<double>& _forces = _readData["Force"];
      std::vector<double>& _forcesPrev = _readDataPrev["Force"];
      for (int j=0; j <_vertexSize; j++)
      {
        int i = _boundaryNodes[j];
        if (dim == 2) 
        { //2D
          force[i] += (1.-fact)*_forcesPrev[dim * j] + fact*_forces[dim * j];
          force[i + numNodes] += (1.-fact)*_forcesPrev[dim * j + 1] + fact*_forces[dim * j + 1];
        }
        else if (dim == 3) 
        { //3D
          force[i] += (1.-fact)*_forcesPrev[dim * j] + fact*_forces[dim * j];
          force[i + numNodes] += (1.-fact)*_forcesPrev[dim * j + 1] + fact*_forces[dim * j + 1];
          force[i + 2 * numNodes] += (1.-fact)*_forcesPrev[dim * j + 2] + fact*_forces[dim * j + 2];
        }
      };
    }
    else if (_forceInterfaceName[dd]=="Stress")
    {
      if (_readData.find("Stress") == _readData.end())
      {
        ERROR("Presure is not present in read data");
        exit(-1);
      }
      const std::vector<double>& Pcur = _readData["Stress"];
      const std::vector<double>& Pprev = _readDataPrev["Stress"];
      for (int i=0; i< _boundaryElements.size(); i++)
      {
        std::vector<classElements* > eleList(1,_boundaryElements[i]);
        eleList[0]->setMe(0);
        std::vector<double> Pf(3,0);
        std::vector<double> Pi(3,0);
        const std::vector<int>& eleNodes = _boundaryElements[i]->getMyNodes();
        int numNodes = eleNodes.size();
        for (int k=0; k< eleNodes.size(); k++)
        {
          int pos = _nodeMap[eleNodes[k]];
          Pf[0] += Pcur[ndim*pos]/numNodes;
          Pf[1] += Pcur[ndim*pos+1]/numNodes;
          Pi[0] += Pprev[ndim*pos]/numNodes;
          Pi[1] += Pprev[ndim*pos+1]/numNodes;
          if (ndim == 3)
          {
            Pf[2] += Pcur[ndim*pos+2]/numNodes;
            Pi[2] += Pprev[ndim*pos+2]/numNodes;
          }
        }
        double Tx = (1-fact)*Pi[0] + fact*Pf[0];
        double Ty = (1-fact)*Pi[1] + fact*Pf[1];
        double Tz = (1-fact)*Pi[2] + fact*Pf[2];
        classTractionRamp* bcTraction = new classTractionRamp("TRACTION RAMP",eleList,Tx, Ty, Tz,timeRun-dt,timeRun,Pi[0],Pi[1],Pi[2]);
        BCsTractionRamp(disp, nod, elem, bcTraction, ndim, force, timeRun, dt);
        delete bcTraction;
      };
    }
    else if (_forceInterfaceName[dd]=="Pressure")
    {
      if (_readData.find("Pressure") == _readData.end())
      {
        ERROR("Presure is not present in read data");
        exit(-1);
      }
      const std::vector<double>& Pcur = _readData["Pressure"];
      const std::vector<double>& Pprev = _readDataPrev["Pressure"];
      
      for (int i=0; i< _boundaryElements.size(); i++)
      {
        std::vector<classElements* > eleList(1,_boundaryElements[i]);
        eleList[0]->setMe(0);
        double Pf(0.), Pi(0);
        const std::vector<int>& eleNodes = _boundaryElements[i]->getMyNodes();
        int numNodes = eleNodes.size();
        for (int k=0; k< eleNodes.size(); k++)
        {
          int pos = _nodeMap[eleNodes[k]];
          Pf += Pcur[pos*ndim]/numNodes;
          Pi += Pprev[pos*ndim]/numNodes;
        }
        double P = (1-fact)*Pi + fact*Pf;
        classPressureRamp* bcTraction = new classPressureRamp("PRESSURE RAMP",eleList,P,timeRun-dt,timeRun,Pi);
        BCsPressureRamp(disp, nod, elem, bcTraction, ndim, force, timeRun, dt, tangentEstimation, stiffBC);
        delete bcTraction;
      };
    }
    else if (_forceInterfaceName[dd]=="StarlingFlux_0")
    {
      if (_readData.find("Pressure") == _readData.end())
      {
        ERROR("Presure is not present in read data");
        exit(-1);
      }
      
      if (_forceInterfaceParameters.find("StarlingFlux_0")==_forceInterfaceParameters.end())
      {
        ERROR("parameters of StarlingFlux cannot be found");
        exit(-1);
      };
      double Lp = _forceInterfaceParameters["StarlingFlux_0"][0];
      INFO("StarlingFlux_0 Lp = %e",Lp);
      
      const std::vector<double>& Pcur = _readData["Pressure"];
      const std::vector<double>& Pprev = _readDataPrev["Pressure"];
      
      for (int i=0; i< _boundaryElements.size(); i++)
      {
        std::vector<classElements* > eleList(1,_boundaryElements[i]);
        eleList[0]->setMe(0);
        double Pf(0.), Pi(0);
        const std::vector<int>& eleNodes = _boundaryElements[i]->getMyNodes();
        int numNodes = eleNodes.size();
        for (int k=0; k< eleNodes.size(); k++)
        {
          int pos = _nodeMap[eleNodes[k]];
          Pf += Pcur[pos*ndim]/numNodes;
          Pi += Pprev[pos*ndim]/numNodes;
        }
        double P = (1-fact)*Pi + fact*Pf;
        classConvection* bcConvection = new classConvection("CONVECTION",eleList,Lp, P, timeRun-dt,timeRun,0);
        BCsConvection(disp, force, nod, elem, bcConvection, ndim, timeRun, dt, tangentEstimation, stiffBC);
        delete bcConvection;
      };
    }
    else
    {
      ERROR("%s has not been implemented",_forceInterfaceName[dd].c_str());
      exit(-1);
    }
  }
};

/*! \brief advance function in preCICE
  @param[in] dt time step
*/
double FSIAdapter::advance(double dt)
{
  return _precice->advance(dt);
};

bool FSIAdapter::isCouplingOngoing()
{
  return _precice->isCouplingOngoing();
}

/*! \brief finalize when finishing coupling
*/
void FSIAdapter::finalize()
{
    _precice->finalize();
};

const std::string& FSIAdapter::actionReadIterationCheckpoint()
{
  return precice::constants::actionReadIterationCheckpoint(); 
}
const std::string& FSIAdapter::actionWriteIterationCheckpoint()
{
  return precice::constants::actionWriteIterationCheckpoint();
}

bool FSIAdapter::isActionRequired(const std::string& action)
{
  return _precice->isActionRequired(action);
}
void FSIAdapter::markActionFulfilled(const std::string& action)
{
  _precice->markActionFulfilled(action);
}

/*! \brief make a checkpoint of the solver
  @param[in] timeRun time
  @param[in] dt time step
*/
void FSIAdapter::makeACheckPoint(double timeRun, double dt)
{
  _timeIterFSI = timeRun;
  _dtIterFSI = dt;
};

/*! \brief load a checkpoint of the solver
  @param[out] timeRun time
  @param[out] dt time step
*/
void FSIAdapter::restoreACheckPoint(double& timeRun, double& dt)
{
  timeRun = _timeIterFSI;
  dt = _dtIterFSI;
};

#endif //FSI