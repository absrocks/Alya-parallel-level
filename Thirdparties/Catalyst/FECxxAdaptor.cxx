
// Adaptor for getting Fortran simulation code into ParaView CoProcessor.

// CoProcessor specific headers
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkImageData.h"

// Fortran specific header
#include "vtkCPPythonAdaptorAPI.h"



#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPAdaptorAPI.h"
#include "vtkFloatArray.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"

//#include <cstdint>
                                   


namespace
{
  vtkSmartPointer<vtkUnstructuredGrid> VTKGrid;
  //-------------------------------------------------------------------------  
  void BuildVTKGrid(unsigned int numberOfPoints, 
		    double* pointsData,
		    unsigned int numberOfCells,
		    unsigned int* cellsTypes,
		    unsigned int icount,
		    uint64_t* cellsData,
		    uint64_t* offset)
    
  {   
    // create the points information
    vtkNew<vtkDoubleArray> pointArray;
    pointArray->SetNumberOfComponents(3);
    pointArray->SetArray(pointsData, static_cast<vtkIdType>(numberOfPoints*3), 1);
    
    vtkNew<vtkPoints> points;
    points->SetData(pointArray.GetPointer());
    
    VTKGrid=vtkSmartPointer<vtkUnstructuredGrid>::New();
    VTKGrid->SetPoints(points.GetPointer());
    /////////////////////////////////////////////////////

    std::vector<vtkIdType> Cells(cellsData, cellsData + icount);
    VTKGrid->Allocate( static_cast<vtkIdType>( Cells.size() ) );

    //std::cout<< icount <<" = icount \n";
    //std::cout<< Cells.size() <<" = Cells.size() \n";
    //std::cout<<"\n";
    //std::cout<<"VTK_TETRA \n";
    //std::cout<< " i = "<< i <<"  \n";
    //std::cout<<" Cells.data()+n_cell*i+1=" << *Cells.data()+n_cell*i+1 <<" \n";
    //std::cout<<"\n";
    //std::cout<<" offset[i]= "<< offset[i] <<" \n";
    //std::cout<<" count= "<< count <<" \n";
    //std::cout<<" offset[i]+count=" << offset[i]+count <<" \n";º   

    
    for(int i=0 ; i<numberOfCells; i++ )
      {	
	if(cellsTypes[i]==10)
	  {
	    int n_cell = 4;
	    VTKGrid->InsertNextCell(VTK_TETRA, n_cell,  Cells.data()+ offset[i]);
	  }
	if(cellsTypes[i]==12)
	  {
	    int n_cell = 8;
	    VTKGrid->InsertNextCell(VTK_HEXAHEDRON, n_cell,  Cells.data()+ offset[i]);
	  }
	if(cellsTypes[i]==13)
	  {
	    int n_cell = 6;
	    VTKGrid->InsertNextCell(VTK_WEDGE, n_cell,  Cells.data()+ offset[i]);
	  }
	if(cellsTypes[i]==14)
	  {
	    int n_cell = 5;
	    VTKGrid->InsertNextCell(VTK_PYRAMID, n_cell,  Cells.data()+ offset[i]);
	  }
      }
   
	///////////////////////////////////////////////////
  }
    //------------------------------------------------------------------------------
    void UpdateVTKAttributes(unsigned int numberOfPoints,
			     double* velocityData,
			     double* pressureData)
    {
      if(VTKGrid->GetPointData()->GetNumberOfArrays() == 0)
	{
	  // velocity array
	  vtkNew<vtkDoubleArray> velocity;
	  velocity->SetName("velocity");
	  velocity->SetNumberOfComponents(3);
	  velocity->SetNumberOfTuples(static_cast<vtkIdType>(numberOfPoints));
	  VTKGrid->GetPointData()->AddArray(velocity.GetPointer());
	  //	}
      //if(VTKGrid->GetCellData()->GetNumberOfArrays() == 0)
      //	{
	  // pressure array
	  vtkNew<vtkDoubleArray> pressure;
	  pressure->SetName("pressure");
	  pressure->SetNumberOfComponents(1);
	  //VTKGrid->GetCellData()->AddArray(pressure.GetPointer());// ASSIGNED BY CELL
	  VTKGrid->GetPointData()->AddArray(pressure.GetPointer());// ASSIGNED BY POINT
	}
      vtkDoubleArray* velocity = vtkDoubleArray::SafeDownCast(
							      VTKGrid->GetPointData()->GetArray("velocity"));
      // The velocity array is ordered as vx0,vx1,vx2,..,vy0,vy1,vy2,..,vz0,vz1,vz2,..
      // so we need to create a full copy of it with VTK's ordering of
      // vx0,vy0,vz0,vx1,vy1,vz1,..
      for(unsigned int i=0;i<numberOfPoints;i++)
	{
	  double values[3] = {velocityData[i], velocityData[i+numberOfPoints],
			      velocityData[i+2*numberOfPoints]};
	  velocity->SetTupleValue(i, values);
	}

      vtkDoubleArray* pressure = vtkDoubleArray::SafeDownCast(
							    VTKGrid->GetPointData()->GetArray("pressure"));
      // The pressure array is a scalar array so we can reuse
      // memory as long as we ordered the points properly.
      pressure->SetArray(pressureData, static_cast<vtkIdType>(numberOfPoints), 1);
    }
    //------------------------------------------------------------------------------------------------
    void BuildVTKDataStructures(unsigned int numberOfPoints, 
				double* pointsData,
				unsigned int numberOfCells,
				unsigned int* cellsTypes,
				unsigned int icount,
				uint64_t* cellsData,
				uint64_t* offset,
				double* velocity,
				double* pressure)

    {
      if(VTKGrid == NULL)
	{
	  // The grid structure isn't changing so we only build it
	  // the first time it's needed. If we needed the memory
	  // we could delete it and rebuild as necessary.
	  VTKGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	  BuildVTKGrid(numberOfPoints, pointsData, numberOfCells,
		       cellsTypes, icount, cellsData, offset);
	  
	}
      //----------------------------------------------------------------------------------------------
      UpdateVTKAttributes(numberOfPoints, velocity, pressure);
      
    }
  }
  
extern "C" void catalystcoprocess_(unsigned int* numberOfPoints, 
				   double* pointsData,
				   unsigned int* numberOfCells, 
				   unsigned int* cellsTypes,
				   unsigned int* icount,
				   uint64_t* cellsData,
				   uint64_t* offset,
				   double* time,
				   unsigned int* timeStep,
				   double* pressure,
				   double* velocity)
{
  /*   
  cout << "numberOfPoints = " <<  *numberOfPoints << endl; 
  for(unsigned int i=0;i<*numberOfPoints;i++)
    {
      cout << "pointsData[I] = " << i << "||||" <<  pointsData[i] << endl;
    }

  cout << "numberOfCells = " <<  *numberOfCells << endl; 
  for(unsigned int i=0;i<*numberOfCells;i++)
    {
      cout << " cellsTypes[i]= " << i << "||||" <<  cellsTypes[i] << endl;
      cout << " offset[i]= " << i << "||||" <<  offset[i] << endl;
    }
  
  cout << "icount = " <<  *icount << endl; 
  for(unsigned int i=0;i<*icount;i++)
    {
      cout << " cellsData[i]= " << i << "||||" <<  cellsData[i] << endl;
    }
  */

 

  vtkCPProcessor* processor = vtkCPAdaptorAPI::GetCoProcessor();
  vtkCPDataDescription* dataDescription = vtkCPAdaptorAPI::GetCoProcessorData();
  

  if(processor == NULL || dataDescription == NULL)
    {
    cerr << "ERROR: Catalyst not properly initialized.\n";
    return;
    }

  dataDescription->AddInput("input");
  dataDescription->SetTimeData(*time, *timeStep);

  int lastTimeStep = 0;
  
  if(lastTimeStep == true)
  {
    //assume that we want to all the pipelines to execute if it
    // is the last time step.
    dataDescription->ForceOutputOn();
  }
  if(processor->RequestDataDescription(dataDescription) != 0)
    {
      BuildVTKDataStructures(*numberOfPoints, pointsData, *numberOfCells,
			     cellsTypes, *icount, cellsData, offset, velocity, pressure);
   
    dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
    processor->CoProcess(dataDescription);
    }
 
}
