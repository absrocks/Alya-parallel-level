/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfInterface.cxx,v 1.1 2009-12-17 18:17:29 kwleiter Exp $  */
/*  Date : $Date: 2009-12-17 18:17:29 $ */
/*  Version : $Revision: 1.1 $ */
/*                                                                 */
/*  Author:                                                        */
/*     Kenneth Leiter                                              */
/*     kenneth.leiter@arl.army.mil                                 */
/*     US Army Research Laboratory                                 */
/*     Aberdeen Proving Ground, MD                                 */
/*                                                                 */
/*     Copyright @ 2009 US Army Research Laboratory                */
/*     All Rights Reserved                                         */
/*     See Copyright.txt or http://www.arl.hpc.mil/ice for details */
/*                                                                 */
/*  Modified by:                                                   */
/*     Raul de la Cruz                                             */
/*     delacruz@bsc.es                                             */
/*     Barcelona Supercomputing Center - CNS                       */
/*     Barcelona, Spain                                            */
/*     February, March, April, May 2013                            */
/*                                                                 */
/*     This software is distributed WITHOUT ANY WARRANTY; without  */
/*     even the implied warranty of MERCHANTABILITY or FITNESS     */
/*     FOR A PARTICULAR PURPOSE.  See the above copyright notice   */
/*     for more information.                                       */
/*                                                                 */
/*******************************************************************/

#include <Xdmf.h>
#include <XdmfSet.h>

#include <stdio.h>
#include <sstream>
#include <map>
#include <stack>
#include <vector>

#include "XdmfInterface.h"


/**
 *
 * Initialize a new Xdmf file.
 *
 */
XdmfInterface::XdmfInterface(char * outputName, int externalHDF)
{
  myDOM = new XdmfDOM();
  myRoot = new XdmfRoot();
  myDomain = new XdmfDomain();
  myTopology = NULL;
  myGeometry = NULL;
  currentTime = -1;
  Debug = XDMF_FALSE;

  myRoot->SetDOM(myDOM);
  myRoot->Build();
  myRoot->Insert(myDomain);
  myName = outputName;
  this->externalHDF = externalHDF;
}

XdmfInterface::~XdmfInterface()
{
  this->Destroy();
}

/**
 *
 * Destroy objects from XdmfInterface.
 *
 */
void
XdmfInterface::Destroy()
{
  currentTime = -1;
  externalHDF = -1;

  /* myGeometry and myTopology are deleted
   * by each Uniform Grid destructor class */
  myGeometry = NULL;
  myTopology = NULL;

  while (!myDataItems.empty())
  {
    delete myDataItems.top();
    myDataItems.pop();
  }

  while (!myAttributes.empty())
  {
    delete myAttributes.back();
    myAttributes.pop_back();
  }

  while (!myInformations.empty())
  {
    delete myInformations.back();
    myInformations.pop_back();
  }

  while (!myCollections.empty())
  {
    delete myCollections.top();
    myCollections.pop();
  }

  delete myDOM;
  delete myRoot;
  delete myDomain;
}

/**
 *
 * Set a time to be assigned to the next grid.
 *
 */
void
XdmfInterface::SetTime(XdmfFloat64 t)
{
  currentTime = t;
}

/**
 *
 * Return the currentTime
 *
 */
XdmfFloat64
XdmfInterface::GetTime()
{
  return currentTime;
}

/**
 *
 * Set external HDF flag.
 *
 */
void
XdmfInterface::SetExternalHDF(XdmfInt32 externalHDF)
{
  this->externalHDF = externalHDF;
}

/**
 *
 * Return the externalHDF flag.
 *
 */
XdmfInt32
XdmfInterface::GetExternalHDF()
{
  return externalHDF;
}

/**
 *
 * Set the Debug status
 *
 */
void
XdmfInterface::SetDebug(XdmfInt32 value)
{
  Debug = value;
  SetGlobalDebug( value );
}

/**
 *
 * Return the Debug status
 *
 */
XdmfInt32
XdmfInterface::GetDebug()
{
  return Debug;
}

/**
 *
 * Add a regular Grid or Collection to the XdmfDOM.
 * Collections can be 'Spatial' or 'Temporal' type.
 * Nested collections are supported.
 *
 * GridType can be:
 *    XDMF_GRID_UNIFORM
 *    XDMF_GRID_COLLECTION
 *    XDMF_GRID_TREE
 *    XDMF_GRID_SUBSET
 *    XDMF_GRID_UNSET
 *
 * CollectionType can be:
 *    XDMF_GRID_COLLECTION_TEMPORAL
 *    XDMF_GRID_COLLECTION_SPATIAL
 *    XMDF_GRID_COLLECTION_UNSET
 *
 */
XdmfInt32
XdmfInterface::AddGrid(char * gridName, XdmfInt32 gridType, XdmfInt32 collectionType)
{
  /* Sanity check */
  /* Makes no sense to add any kind of Grid inside an Uniform Grid or other
   * collection element. Uniform Grid must be the leaf of each Grid branch */
  if (!myCollections.empty()) {
    if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Grid") == 0) {
      XdmfErrorMessage("Cannot add a Grid inside a non Grid element (" <<
                       myCollections.top()->GetElementType() << ")");
      return(XDMF_FAIL);
    }
    if (((XdmfGrid*)myCollections.top())->GetGridType() == XDMF_GRID_UNIFORM) {
      XdmfErrorMessage("Cannot add a Grid inside an Uniform Grid");
      return(XDMF_FAIL);
    }
  }

  /* Set Grid to be deleted by parents by default, however
   * this grid will be deleted by CloseGrid if called before */
  XdmfGrid * currentGrid = new XdmfGrid();
  currentGrid->SetGridType(gridType);
  currentGrid->SetDeleteOnGridDelete(true);
  XdmfDebug("Adding Grid " << gridName << " as type " << currentGrid->GetGridTypeAsString());

  if (gridType == XDMF_GRID_COLLECTION) {
    if (currentGrid->SetCollectionType(collectionType) == XDMF_FAIL) {
      XdmfErrorMessage("Can not set collection type to " << collectionType);
      return(XDMF_FAIL);
    }
  }


  /* If we try to write over the same grid name, modify the current grid name */
  std::stringstream totalGridName;
  std::string stdGridName = gridName;
  if (myGridNames.find(stdGridName) == myGridNames.end())
  {
    myGridNames[stdGridName] = 1;
    totalGridName << gridName;
  }
  else
  {
    myGridNames[stdGridName]++;
    totalGridName << gridName << "_" << myGridNames[stdGridName];
  }
  currentGrid->SetName(totalGridName.str().c_str());


  /* Add new Grid to the collection stack and add it to the DOM structure */
  if (myCollections.empty())
    myDomain->Insert(currentGrid);
  else
    ((XdmfGrid*)myCollections.top())->Insert(currentGrid);

  currentGrid->Build();
  myCollections.push(currentGrid);

  /* Set as the last Grid added to the collection */
  lastGrid = currentGrid;


  /* DOM and XdmfXmlNode Element are transfered to/created in children after calling
   * Parent->Insert( children ) and XdmfElement::Insert( ) function members
   *
   * Insert call: Create XdmfXmlNode objects and set element to a specific type,
   *              Transfers also DOM to descendents
   * Build call: Traverses DOM tree (XdmfXmlNodes) and fills their attributes,
   *             and CDATA for LightData and H5 files for HeavyData */

  /* Be careful: Topology and Geometry elements are added automatically by Xdmf
   * internals when Grid is inserted into a Domain or another Grid */

  return(XDMF_SUCCESS);
}

/**
 *
 * Close the current open collection.  If within a nested collection, close
 * the most deeply nested collection.
 *
 * Although this call is not mandatory because resources are freed when
 * ~XdmfInterface is called, it is recommendable to call
 * CloseGrid after WriteGrid and finishing to build the grid.
 *
 */
int
XdmfInterface::CloseGrid()
{
  /* Sanity check - top element is a Grid? */
  /* Check if there is at least one Grid/Attribute in the collection */
  if (myCollections.empty()) {
    XdmfErrorMessage("No Grid added, please add them to the collection before continue");
    return(XDMF_FAIL);
  }

  /* Close previous elements till reach the first Grid in stack */
  if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Grid") == 0) {
    XdmfErrorMessage("Top element in Collection stack is not a Grid ("
                     << myCollections.top()->GetElementType() << ") : " <<
                     "check your call tree of Add/Close* methods to XdmfInterface");
    XdmfErrorMessage("Recovering status by closing any other opened elements");

    /* Don't delete elements, only pop from stack (they are deleted by parent Grid) */
    while ( !myCollections.empty() && (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Grid") == 0 )) {
      myCollections.pop();
    }
  }

  /* Destroy current top Grid */
  /* It is better to destroy Grid here in order to free memory as soon as
   * Grid DOM/Xml data has been built (after WriteGrid) instead of flagging
   * SetGridOnDelete(1) to be destroyed after parent Grid is deleted */
  if (!myCollections.empty()) {
    delete ((XdmfGrid*)myCollections.top());
    myCollections.pop();

    /* Remove the last Grid in parent Grid element (Children list) */
    if (!myCollections.empty()) {
      ((XdmfGrid*)myCollections.top())->SetNumberOfChildren(((XdmfGrid*)myCollections.top())->GetNumberOfChildren() - 1);

      /* Set the parent of the closed Grid as the lastGrid opened */
      lastGrid = (XdmfGrid*)myCollections.top();
    }
  }
  else {
    XdmfErrorMessage("Trying to close/delete a non-existent Grid");
    return(XDMF_FAIL);
  }
  
  return(XDMF_SUCCESS);
}


/**
 *
 * Set the topology type to be assigned to the next grid.
 * XDMF_INT_32 or XDMF_INT_64 type are supported for Topology --> INTEGER*4 | INTEGER*8
 *
 */
XdmfInt32
XdmfInterface::SetGridTopology(void ** reference, XdmfInt32 topologyType, XdmfInt32 numberType, XdmfInt32 noeRank,
                               XdmfInt64 * numberOfElements, XdmfInt32 dimRank, XdmfInt64 * dimensions, XdmfPointer conns)
{
  std::stringstream shapeElems, shapeDims;

  if (numberOfElements) {
    for (int i= 0; i< noeRank; i++)
      shapeElems << numberOfElements[i] << " ";
    XdmfDebug("Shape for NumberOfElements is: " << shapeElems.str());
  }
  else shapeElems << "";

  if (dimensions) {
    for (int i= 0; i< dimRank; i++)
      shapeDims << dimensions[i] << " ";
    XdmfDebug("Shape for Dimensions is: " << shapeDims.str());
  }
  else shapeDims << "";

  return this->SetGridTopologyFromShape(reference, topologyType, numberType, (char*) shapeElems.str().c_str(),
                                        (char*) shapeDims.str().c_str(), conns);
}


/**
 *
 * Set the topology type to be assigned to the next grid.
 * XDMF_INT_32 or XDMF_INT_64 type are supported for Topology --> INTEGER*4 | INTEGER*8
 *        XDMF_NOTOPOLOGY     0x0
 *        XDMF_POLYVERTEX     0x1
 *        XDMF_POLYLINE       0x2
 *        XDMF_POLYGON        0x3
 *        XDMF_TRI            0x4
 *        XDMF_QUAD           0x5
 *        XDMF_TET            0x6
 *        XDMF_PYRAMID        0x7
 *        XDMF_WEDGE          0x8
 *        XDMF_HEX            0x9
 *        XDMF_EDGE_3         0x0022
 *        XDMF_TRI_6          0x0024
 *        XDMF_QUAD_8         0x0025
 *        XDMF_TET_10         0x0026
 *        XDMF_PYRAMID_13     0x0027
 *        XDMF_WEDGE_15       0x0028
 *        XDMF_WEDGE_18       0x0029
 *        XDMF_HEX_20         0x0030
 *        XDMF_HEX_24         0x0031
 *        XDMF_HEX_27         0x0032
 *        XDMF_MIXED          0x0070
 *        XDMF_2DSMESH        0x0100
 *        XDMF_2DRECTMESH     0x0101
 *        XDMF_2DCORECTMESH   0x0102
 *        XDMF_3DSMESH        0x1100
 *        XDMF_3DRECTMESH     0x1101
 *        XDMF_3DCORECTMESH   0x1102
 *
 *
 * Structured:
 *  * 2DSMesh, 2DRectMesh, 2DCORECTMesh, 3DSMesh, 3DRectMesh and 3DCORECTMesh
 *      - NumberOfElements: NumberOfElements (topology is implicit)
 *
 * Unstructured:
 *  * POLY* topologies
 *      - NodesPerElement:  NodesPerElement
 *      - NumberOfElements: NumberOfElements
 *      - Dimensions:       Dimensions (NodesPerElement * NumberOfElements)
 *
 *  * Mixed topologies
 *      - NumberOfElements: NumberOfElements
 *      - Dimensions:       Total number of nodes in mesh including elements descriptors
 *
 *  * Triangle, Quadrilateral, Tetrahedron, Pyramid, Wedge, Hexahedron,
 *    Edge_3, Triangle_6, Quadrilateral_8, Tetrahedron_10, Pyramid_13,
 *    Wedge_15, Wedge_18, Hexahedron_20, Hexahedron_24, Hexahedron_27
 *      - NumberOfElements: NumberOfElements (Dimensions/NodesPerElement)
 *      - Dimensions:       Total number of nodes in mesh including elements descriptors
 *
 */
XdmfInt32
XdmfInterface::SetGridTopologyFromShape(void ** reference, XdmfInt32 topologyType, XdmfInt32 numberType,
                                        char * numberOfElements, char * dimensions, XdmfPointer conns)
{
  /* Sanity check */
  /* Check if there is at least one Grid in the collection */
  if (myCollections.empty()) {
    XdmfErrorMessage("Uniform Grid not added, please add one to the collection before continue");
    return(XDMF_FAIL);
  }

  /* Check if last element added to collection stack is a Grid */
  if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Grid") == 0) {
    XdmfErrorMessage("Cannot set Topology inside a non Grid element (" <<
                     myCollections.top()->GetElementType() << ")");
    return(XDMF_FAIL);
  }

  /* Check if last Grid added is not a collection type and is uniform type */
  if (((XdmfGrid*)myCollections.top())->GetGridType() != XDMF_GRID_UNIFORM) {
    XdmfErrorMessage("Cannot set Topology to a Collection type Grid, Grid must be Uniform type");
    return(XDMF_FAIL);
  }

  /* Check if topology was already set, never delete topology object, it belongs to Grid object */
  if (!myTopology) myTopology = ((XdmfGrid*)myCollections.top())->GetTopology();

  /* Topology */
  myTopology->SetTopologyType(topologyType);

  /* NumberOfElements */
  myTopology->GetShapeDesc()->SetShapeFromString(numberOfElements);
  XdmfDebug("NumberOfElements fixed to " << myTopology->GetNumberOfElements() << " for " <<
             myTopology->GetTopologyTypeAsString() << " type");


  /* Sanity check - NodesPerElement | NumberOfElements | Dimensions */
  /* Compare topology element nodes and number of elements in data to write */
  XdmfDataDesc shape;
  switch( myTopology->GetTopologyType() ) {
    case XDMF_2DSMESH:
    case XDMF_2DRECTMESH:
    case XDMF_2DCORECTMESH:
      /* Sanity check for 2D structured meshes */
      shape.SetShapeFromString(numberOfElements);
      if (shape.GetRank() != 2) {
        XdmfErrorMessage( "Topology type " << myTopology->GetTopologyTypeAsString() << " require a Dimension " <<
                          "Rank of 2 for NumberOfElements. Current NumberOfElement dimensions are: " << numberOfElements );
        return(XDMF_FAIL);
      }
      break;
    case XDMF_3DSMESH:
    case XDMF_3DRECTMESH:
    case XDMF_3DCORECTMESH:
      /* Sanity check for 3D structured meshes */
      shape.SetShapeFromString(numberOfElements);
      if (shape.GetRank() != 3) {
        XdmfErrorMessage( "Topology type " << myTopology->GetTopologyTypeAsString() << " require a Dimension " <<
                          "Rank of 3 for NumberOfElements. Current NumberOfElement dimensions are: " << numberOfElements );
        return(XDMF_FAIL);
      }
      break;
    case XDMF_MIXED:
      /* Requires communication among processors, maybe is better not to do it because performance issues */
      // GetCellOffsets()
      break;
    case XDMF_POLYVERTEX:
    case XDMF_POLYLINE:
    case XDMF_POLYGON: // POLY* elements
      /* NodesPerElement (Only for POLY* topology) - Try to guess */
      shape.SetShapeFromString(dimensions);
      myTopology->SetNodesPerElement( shape.GetNumberOfElements() / myTopology->GetNumberOfElements() );
      XdmfDebug("NodesPerElement fixed to " << myTopology->GetNodesPerElement() << " for " <<
                myTopology->GetTopologyTypeAsString() << " type");
    default: // The rest of topology types
      // Dimensions = (NodesPerElement * NumberOfElements)
      shape.SetShapeFromString(dimensions);
      if (shape.GetNumberOfElements() != myTopology->GetNodesPerElement() * myTopology->GetNumberOfElements()) {
        XdmfErrorMessage("Check grid Topology parameters: Dimensions (" << shape.GetShapeAsString() <<
                         ") != NodesPerElement (" << myTopology->GetNodesPerElement() << ") * NumberOfElements (" <<
                         myTopology->GetShapeDesc()->GetShapeAsString() << ") for " << myTopology->GetTopologyTypeAsString() << " type");
        return(XDMF_FAIL);
      }
      break;
  }


  /* Use reference XdmfXmlNode object as reference for new object */
  if ((reference) && (*reference != NULL)) {
    XdmfXmlNode xmlNode = (XdmfXmlNode) *reference;

    if (myTopology->CheckForReference(myTopology->GetElement(), xmlNode) != xmlNode) {
      XdmfErrorMessage("Error checking object Reference for Topology " << myTopology << ". Falling in default mode");

      /* Fallback to non-reference mode */
      *reference = NULL;
    }
    else XdmfDebug("Setting Topology " << myTopology << " Reference as " << myDOM->GetPath(xmlNode) << " of XdmfXmlNode " << xmlNode);
  }

  /* Otherwise build a new object if topology is Unstructured */
  if (((!reference) || (reference && !*reference)) && (myTopology->GetClass() != XDMF_STRUCTURED)) {
    if (externalHDF) {
      // Insert DataItem element from an existing HDF5 dataset
      XdmfInjectXml( myTopology, conns );
    }
    else {
      XdmfArray * myConnections = myTopology->GetConnectivity();
      myConnections->SetShapeFromString( dimensions );

      /* Set Array data */
      if (numberType == XDMF_INT32_TYPE) {
        myConnections->SetNumberType(XDMF_INT32_TYPE);
        myConnections->SetValues(0, (XdmfInt32*)conns, myConnections->GetNumberOfElements());
      }
      else if (numberType == XDMF_INT64_TYPE) {
        myConnections->SetNumberType(XDMF_INT64_TYPE);
        myConnections->SetValues(0, (XdmfInt64*)conns, myConnections->GetNumberOfElements());
      }
      else {
        XdmfErrorMessage("Number Type: " << numberType << " is not valid");
        return(XDMF_FAIL);
      }
    }
  }

  /* Topology built-in Grid Class is used, so there is no need to add this object to Grid */
  /* Add Topology to the outter most Uniform Grid */
  // myCollections.top()->SetTopology(myTopology); /* Add Topology into Xdmf Grid class */
  // myCollections.top()->Insert(myTopology);      /* Add Topology into DOM tree */
  // myCollections.top()->Build();                 /* Add XML/DOM attributes to the XdmfXmlNode object */

  /* Return XdmfXmlNode reference to this object */
  if (reference) {
    *reference = (void*) myTopology->GetElement();
    XdmfDebug("Returning Topology reference " << *reference);
  }

  return(XDMF_SUCCESS);
}


/**
 *
 * Set the geometry type to be assigned to the next grid.
 *
 */
XdmfInt32
XdmfInterface::SetGridGeometry(void ** reference, XdmfInt32 geometryType, XdmfInt32 numberType,
                               XdmfInt32 rank, XdmfInt64 * dimensions, XdmfPointer points)
{
  std::stringstream shape;

  if (dimensions) {
    for (int i= 0; i< rank; i++)
      shape << dimensions[i] << " ";
    XdmfDebug("Shape for Dimensions is: " << shape.str());
  }
  else shape << "";

  return this->SetGridGeometryFromShape(reference, geometryType, numberType,
                                        (char*) shape.str().c_str(), points);
}


/**
 *
 * Set the geometry type to be assigned to the next grid.
 *        XDMF_GEOMETRY_NONE          0
 *        XDMF_GEOMETRY_XYZ           1
 *        XDMF_GEOMETRY_XY            2
 *        XDMF_GEOMETRY_X_Y_Z         3
 *        XDMF_GEOMETRY_X_Y           4
 *        XDMF_GEOMETRY_VXVYVZ        5
 *        XDMF_GEOMETRY_ORIGIN_DXDYDZ 6
 *        XDMF_GEOMETRY_VXVY          7
 *        XDMF_GEOMETRY_ORIGIN_DXDY   8
 *
 * HDF5URLs: HDF5URL;HDF5URL;HDF5URL
 *
 */
XdmfInt32
XdmfInterface::SetGridGeometryFromShape(void ** reference, XdmfInt32 geometryType, XdmfInt32 numberType,
                                        char * dimensions, XdmfPointer points)
{
  /* Sanity check */
  /* Check if there is at least one Grid in the collection */
  if (myCollections.empty()) {
    XdmfErrorMessage("Uniform Grid not added, please add one to the collection before continue");
    return(XDMF_FAIL);
  }

  /* Check if last element added to collection stack is a Grid */
  if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Grid") == 0) {
    XdmfErrorMessage("Cannot set Geometry inside a non Grid element (" <<
                     myCollections.top()->GetElementType() << ")");
    return(XDMF_FAIL);
  }

  /* Check if last Grid added is not a collection type and is uniform type */
  if (((XdmfGrid*)myCollections.top())->GetGridType() != XDMF_GRID_UNIFORM) {
    XdmfErrorMessage("Cannot set Geometry to a Collection type Grid, Grid must be Uniform type");
    return(XDMF_FAIL);
  }

  /* Check if geometry was already set, never delete topology object, it belongs to Grid object */
  if (!myGeometry) myGeometry = ((XdmfGrid*)myCollections.top())->GetGeometry();

  /* Sanity check between Geometry and Topology data */
  if (!myTopology) {
    XdmfErrorMessage("Topology must be set before than Geometry for the Uniform Grid");
    return(XDMF_FAIL);
  }
  XdmfInt32 topologyType = myTopology->GetTopologyType();
  switch(geometryType) {
    case XDMF_GEOMETRY_XY:
    case XDMF_GEOMETRY_X_Y:
      if ((topologyType != XDMF_POLYVERTEX) && (topologyType != XDMF_POLYLINE) && (topologyType != XDMF_POLYGON) &&
          (topologyType != XDMF_TRI) && (topologyType != XDMF_QUAD) && (topologyType != XDMF_2DSMESH)) {
        XdmfErrorMessage("Topology must be [XDMF_POLYVERTEX | XDMF_POLYLINE | XDMF_POLYGON | XDMF_TRI | XDMF_QUAD | XDMF_2DSMESH]" <<
                         "in order to support XDMF_GEOMETRY_XY or XDMF_GEOMETRY_X_Y geometries (Current Topology is " <<
                         myTopology->GetTopologyTypeAsString() << ")" );
        return(XDMF_FAIL);
      }
      break;
    case XDMF_GEOMETRY_XYZ:
    case XDMF_GEOMETRY_X_Y_Z:
      if ((myTopology->GetClass() == XDMF_STRUCTURED) && (topologyType != XDMF_3DSMESH)) {
        XdmfErrorMessage("Topology must be [XDMF_UNSTRUCTURED | XDMF_3DSMESH]" <<
                         "in order to support XDMF_GEOMETRY_XYZ or XDMF_GEOMETRY_X_Y_Z geometries (Current Topology is " <<
                         myTopology->GetTopologyTypeAsString() << ")" );
        return(XDMF_FAIL);
      }
      break;
    case XDMF_GEOMETRY_VXVY:
      if (topologyType != XDMF_2DRECTMESH) {
        XdmfErrorMessage("Topology must be [XDMF_2DRECTMESH]" <<
                         "in order to support XDMF_GEOMETRY_VXVY geometry (Current Topology is " <<
                         myTopology->GetTopologyTypeAsString() << ")" );
        return(XDMF_FAIL);
      }
      break;
    case XDMF_GEOMETRY_VXVYVZ:
      if (topologyType != XDMF_3DRECTMESH) {
        XdmfErrorMessage("Topology must be [XDMF_3DRECTMESH]" <<
                         "in order to support XDMF_GEOMETRY_VXVYVZ geometry (Current Topology is " <<
                         myTopology->GetTopologyTypeAsString() << ")" );
        return(XDMF_FAIL);
      }
      break;
    case XDMF_GEOMETRY_ORIGIN_DXDY:
      if (topologyType != XDMF_2DCORECTMESH) {
        XdmfErrorMessage("Topology must be [XDMF_2DCORECTMESH]" <<
                         "in order to support XDMF_GEOMETRY_ORIGIN_DXDY geometry (Current Topology is " <<
                         myTopology->GetTopologyTypeAsString() << ")" );
        return(XDMF_FAIL);
      }
      break;
    case XDMF_GEOMETRY_ORIGIN_DXDYDZ:
      if (topologyType != XDMF_3DCORECTMESH) {
        XdmfErrorMessage("Topology must be [XDMF_3DCORECTMESH]" <<
                         "in order to support XDMF_GEOMETRY_ORIGIN_DXDYDZ geometry (Current Topology is " <<
                         myTopology->GetTopologyTypeAsString() << ")" );
        return(XDMF_FAIL);
      }
      break;
  }
  myGeometry->SetGeometryType(geometryType);

  /* Set Geometry Shape */
  XdmfDataDesc shape;
  shape.SetShapeFromString(dimensions);
  myGeometry->SetNumberOfPoints(shape.GetNumberOfElements());

  /* Use reference XdmfXmlNode object as reference for new object */
  if ((reference) && (*reference != NULL)) {
    XdmfXmlNode xmlNode = (XdmfXmlNode) *reference;

    if (myGeometry->CheckForReference(myGeometry->GetElement(), xmlNode) != xmlNode) {
      XdmfErrorMessage("Error checking object Reference for Geometry " << myGeometry << ". Falling in default mode");

      /* Fallback to non-reference mode */
      *reference = NULL;
    }
    else XdmfDebug("Setting Geometry " << myGeometry << " Reference as " << myDOM->GetPath(xmlNode) << " of XdmfXmlNode " << xmlNode);
  }

  /* Otherwise build a new object */
  if ((!reference) || (reference && !*reference)) {
    if (myGeometry->GetGeometryType() == XDMF_GEOMETRY_ORIGIN_DXDYDZ) {

      if (myGeometry->GetNumberOfPoints() != 6) {
        XdmfErrorMessage("Geometry Type XDMF_GEOMETRY_ORIGIN_DXDYDZ must have 6 points (Ox,Oy,Oz) (Dx,Dy,Dz)");
        return(XDMF_FAIL);
      }

      switch(numberType) {
        case XDMF_FLOAT32_TYPE:
          myGeometry->SetOrigin( ((XdmfFloat32*)points)[0], ((XdmfFloat32*)points)[1], ((XdmfFloat32*)points)[2] );
          myGeometry->SetDxDyDz( ((XdmfFloat32*)points)[3], ((XdmfFloat32*)points)[4], ((XdmfFloat32*)points)[5] );
          break;
        case XDMF_FLOAT64_TYPE:
          myGeometry->SetOrigin( ((XdmfFloat64*)points)[0], ((XdmfFloat64*)points)[1], ((XdmfFloat64*)points)[2] );
          myGeometry->SetDxDyDz( ((XdmfFloat64*)points)[3], ((XdmfFloat64*)points)[4], ((XdmfFloat64*)points)[5] );
          break;
        default:
          XdmfErrorMessage("Geometry numberType: %d not supported, must be XDMF_FLOAT32_TYPE or XDMF_FLOAT64_TYPE");
          return(XDMF_FAIL);
      }
    }
    else if (myGeometry->GetGeometryType() == XDMF_GEOMETRY_ORIGIN_DXDY) {

      if (myGeometry->GetNumberOfPoints() != 3) {
        XdmfErrorMessage("Geometry Type XDMF_GEOMETRY_ORIGIN_DXDY must have 3 points (Ox,Oy) (Dx,Dy)");
        return XDMF_FAIL;
      }

      switch(numberType) {
        case XDMF_FLOAT32_TYPE:
          myGeometry->SetOrigin( ((XdmfFloat32*)points)[0], ((XdmfFloat32*)points)[1], 0 );
          myGeometry->SetDxDyDz( ((XdmfFloat32*)points)[2], ((XdmfFloat32*)points)[3], 0 );
          break;
        case XDMF_FLOAT64_TYPE:
          myGeometry->SetOrigin( ((XdmfFloat64*)points)[0], ((XdmfFloat64*)points)[1], 0 );
          myGeometry->SetDxDyDz( ((XdmfFloat64*)points)[2], ((XdmfFloat64*)points)[3], 0 );
          break;
        default:
          XdmfErrorMessage("Geometry numberType: %d not supported, must be XDMF_FLOAT32_TYPE or XDMF_FLOAT64_TYPE");
          return(XDMF_FAIL);
      }
    }
    else if (externalHDF) { // External HDF only makes sense if Geometry is not ORIGIN_DXDYDZ
      // Insert DataItem element from an existing HDF5 dataset

      /* Parse different URL for XDMF_GEOMETRY_X_Y | XDMF_GEOMETRY_X_Y_Z (2 HDF5 URLs)
       * and XDMF_GEOMETRY_VXVY | XDMF_GEOMETRY_VXVYVZ (3 HDF5 URLs) */
      XdmfInjectXml( myGeometry, points );
    }
    else {
      XdmfArray * myPoints = myGeometry->GetPoints();
      myPoints->SetNumberType(numberType);

      /* Dimensions refers to the total number of elements in points array */
      myPoints->SetShapeFromString(dimensions);

      WriteToXdmfArray(myPoints, points);
    }
  }

  /* Geometry built-in Grid Class is used, so there is no need to add this object to Grid */
  /* Add Geometry to the outter most Uniform Grid */
  // myCollections.top()->SetGeometry(myGeometry); /* Add Geometry into Xdmf Grid class */
  // myCollections.top()->Insert(myGeometry);      /* Add Geometry into DOM tree */
  // myCollections.top()->Build();                 /* Add XML/DOM attributes to the XdmfXmlNode object */

  /* Return XdmfXmlNode reference to this object */
  if (reference) {
    *reference = (void*) myGeometry->GetElement();
    XdmfDebug("Returning Geometry reference " << *reference);
  }

  return(XDMF_SUCCESS);
}


/**
 *
 * Add an attribute to be written to the next grid.  Multiple attributes can
 * be added and written to a single grid.
 * // Number Types
 *         XDMF_INT8_TYPE      1
 *         XDMF_INT16_TYPE     6
 *         XDMF_INT32_TYPE     2
 *         XDMF_INT64_TYPE     3
 *         XDMF_FLOAT32_TYPE   4
 *         XDMF_FLOAT64_TYPE   5
 *         XDMF_UINT8_TYPE     7
 *         XDMF_UINT16_TYPE    8
 *         XDMF_UINT32_TYPE    9
 *
 * // Value Types
 *         XDMF_ATTRIBUTE_TYPE_NONE       0x00
 *         XDMF_ATTRIBUTE_TYPE_SCALAR     0x01
 *         XDMF_ATTRIBUTE_TYPE_VECTOR     0x02
 *         XDMF_ATTRIBUTE_TYPE_TENSOR     0x03
 *         XDMF_ATTRIBUTE_TYPE_MATRIX     0x04
 *         XDMF_ATTRIBUTE_TYPE_TENSOR6    0x05
 *         XDMF_ATTRIBUTE_TYPE_GLOBALID   0x06
 *         XDMF_ATTRIBUTE_TYPE_UNIFORM    0x00 // By default Attributes are uniform
 *         XDMF_ATTRIBUTE_TYPE_COLLECTION 0x10 // Attribute composed of several DataItems
 *         XDMF_ATTRIBUTE_TYPE_MASK       0x0F // Evaluates type of Single DataItem
 *
 *    i.e.: XDMF_ATTRIBUTE_TYPE_SCALAR | XDMF_ATTRIBUTE_TYPE_COLLECTION
 *
 * // Where Values are Assigned
 *         XDMF_ATTRIBUTE_CENTER_GRID  0
 *         XDMF_ATTRIBUTE_CENTER_CELL  1
 *         XDMF_ATTRIBUTE_CENTER_FACE  2
 *         XDMF_ATTRIBUTE_CENTER_EDGE  3
 *         XDMF_ATTRIBUTE_CENTER_NODE  4
 *
 */
XdmfInt32
XdmfInterface::AddGridAttribute(void ** reference, char * attributeName, XdmfInt32 numberType, XdmfInt32 attributeCenter,
                                XdmfInt32 attributeType, XdmfInt32 rank, XdmfInt64 *dimensions, XdmfPointer data)
{
  std::stringstream shape;
  for (int i= 0; i< rank; i++)
    shape << dimensions[i] << " ";
  XdmfDebug("Shape is: " << shape.str());
  return this->AddGridAttributeFromShape(reference, attributeName, numberType, attributeCenter,
                                         attributeType, (char*) shape.str().c_str(), NULL, data);
}


/**
 *
 * Add an attribute with shape to be written to the next grid.  Multiple attributes can
 * be added and written to a single grid. Dimensions are specified using shape string rather than int.
 * Extra argument supports setting units string of attribute.
 *
 */
XdmfInt32
XdmfInterface::AddGridAttributeFromShape(void ** reference, char * attributeName, XdmfInt32 numberType,
                                         XdmfInt32 attributeCenter, XdmfInt32 attributeType, char * shape,
                                         char * units, XdmfPointer data)
{
  /* Sanity check */
  /* Check if there is at least one Grid in the collection */
  if (myCollections.empty()) {
    XdmfErrorMessage("Grid not added, please add one to the collection before continue");
    return(XDMF_FAIL);
  }

  /* Check if last element added to collection stack is a Grid */
  if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Grid") == 0) {
    XdmfErrorMessage("Cannot add an Attribute inside a non Grid element (" <<
                     myCollections.top()->GetElementType() << ")");
    return(XDMF_FAIL);
  }

  XdmfAttribute * currAttribute = new XdmfAttribute();
  currAttribute->SetName(attributeName);
  if (units) currAttribute->SetUnits(units);
  currAttribute->SetAttributeCenter(attributeCenter);
  currAttribute->SetAttributeType(attributeType);
  currAttribute->SetDeleteOnGridDelete(true);

  /* Attribute must be inserted before checking reference
   * because XdmfXmlNode is created once is inserted */

  /* Add Attribute to the outter most Grid */
  ((XdmfGrid*)myCollections.top())->Insert(currAttribute);      /* Add Attribute into Xdmf Grid and DOM tree */

  /* Use reference XdmfXmlNode object as reference for new object */
  if ((reference) && (*reference != NULL)) {
    XdmfXmlNode xmlNode = (XdmfXmlNode) *reference;

    if (currAttribute->CheckForReference(currAttribute->GetElement(), xmlNode) != xmlNode) {
      XdmfErrorMessage("Error checking object Reference for Attribute " << currAttribute << ". Falling in default mode");

      /* Fallback to non-reference mode */
      *reference = NULL;
    }
    else XdmfDebug("Setting Attribute " << currAttribute << " Reference as " << myDOM->GetPath(xmlNode) << " of XdmfXmlNode " << xmlNode);
  }

  /* Otherwise build a new object if not a Collection Attribute */
  if (((!reference) || (reference && !*reference)) && (!currAttribute->GetIsMultiple())) {
    if (externalHDF) {
      // Insert DataItem element from an existing HDF5 dataset
      XdmfDebug("Injecting XML dataitem from " << (char*)data);
      XdmfInjectXml( currAttribute, data );
      XdmfDebug("Despues Injecting XML dataitem from " << (char*)data);
    }
    else {
      XdmfArray * array = currAttribute->GetValues();
      array->SetNumberType(numberType);
      if (array->SetShapeFromString(shape) == XDMF_FAIL) { // dims via shape string
        XdmfErrorMessage("Can not set array shape to: " << shape);
        return(XDMF_FAIL);
      }
      WriteToXdmfArray(array, data);
    }
  }

  /* Only build object if it is not going to be a Uniform Grid (leaf of the tree).
   * Otherwise, this Attribute will be built in WriteGrid() call by the Build()
   * method of the Uniform Grid class.
   */
  if (((XdmfGrid*)myCollections.top())->GetGridType() != XDMF_GRID_UNIFORM)
    currAttribute->Build();                    /* Add XML/DOM attributes to the XdmfXmlNode object */
  else myAttributes.push_back(currAttribute);

  /* Add Attribute to the collection stack if Collection element */
  if (currAttribute->GetIsMultiple()) myCollections.push(currAttribute);

  /* Return XdmfXmlNode reference to this object */
  if (reference) {
    *reference = (void*) currAttribute->GetElement();
    XdmfDebug("Returning Attribute reference " << *reference);
  }

  return(XDMF_SUCCESS);
}


/**
 *
 * Close the current open Attribute collection. If within a nested collection, close
 * the most deeply nested collection.
 *
 * Although this call is not mandatory because resources are freed when
 * ~XdmfInterface is called, it is recommendable to call
 * CloseAttribute before closing Grid and WriteGrid.
 *
 */
int
XdmfInterface::CloseAttribute()
{
  /* Sanity check */
  /* Check if there is at least one Grid/Attribute in the collection */
  if (myCollections.empty()) {
    XdmfErrorMessage("Neither Grid nor Attribute added, please add them to the collection before continue");
    return(XDMF_FAIL);
  }

  /* Check if last element added to collection stack is an Attribute */
  if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Attribute") == 0) {
    XdmfErrorMessage("Top element in Collection stack is not an Attribute ("
                     << myCollections.top()->GetElementType() << ") : " <<
                     "check your call tree of Add/Close* methods to XdmfInterface");
    XdmfErrorMessage("Recovering status by closing any other opened elements");

    /* Don't delete elements, only pop from stack (they are deleted by parent Grid) */
    while ( !myCollections.empty() && (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Attribute") == 0 )) {
      myCollections.pop();
    }
  }

  if (!myCollections.empty()) {
    /* Remove current top Attribute */
    /* Attributes are removed by their Grid parents in Destroy function */
    myCollections.pop();
  }
  else {
    XdmfErrorMessage("Trying to close/delete a non-existent Attribute");
    return(XDMF_FAIL);
  }

  return(XDMF_SUCCESS);
}


/**
 *
 * Add a DataItem to be written to the previous Attribute element. Multiple dataitems can
 * be added and written into a single Attribute (Collection Attribute).
 * // Number Types
 *         XDMF_INT8_TYPE      1
 *         XDMF_INT16_TYPE     6
 *         XDMF_INT32_TYPE     2
 *         XDMF_INT64_TYPE     3
 *         XDMF_FLOAT32_TYPE   4
 *         XDMF_FLOAT64_TYPE   5
 *         XDMF_UINT8_TYPE     7
 *         XDMF_UINT16_TYPE    8
 *         XDMF_UINT32_TYPE    9
 *
 * // Item Types
 *         XDMF_ITEM_UNIFORM        0x00
 *         XDMF_ITEM_HYPERSLAB      0x01
 *         XDMF_ITEM_COORDINATES    0x02
 *         XDMF_ITEM_FUNCTION       0x03
 *         XDMF_ITEM_COLLECTION     0x14
 *         XDMF_ITEM_TREE           0x15
 *
 *         XDMF_ITEM_MASK           0xF0    // Evaluates to a Single Array ?
 *
 * // Format available - Only makes sense in XDMF_ITEM_UNIFORM
 *         XDMF_FORMAT_XML      0
 *         XDMF_FORMAT_HDF      1
 *         XDMF_FORMAT_MYSQL    2
 *         XDMF_FORMAT_BINARY   3
 *
 * Reference is not returned back for a DataItem when externalHDF is enabled
 * and it is a XDMF_ITEM_UNIFORM (Single Array) because there is not a real
 * XdmfDataItem added to the XdmfInterface Tree actually.
 * Future improvement? Add a real XdmfDataItem in this case?
 */
XdmfInt32
XdmfInterface::AddDataItem(void ** reference, char * itemName, XdmfInt32 numberType, XdmfInt32 itemType,
                           XdmfInt32 format, char * function, XdmfInt32 rank, XdmfInt64 *dimensions, XdmfPointer data)
{
  std::stringstream shape;
  for (int i= 0; i< rank; i++)
    shape << dimensions[i] << " ";
  XdmfDebug("Shape is: " << shape.str());
  return this->AddDataItemFromShape(reference, itemName, numberType, itemType, format,
                                    function, (char*) shape.str().c_str(), data);
}


XdmfInt32
XdmfInterface::AddDataItemFromShape(void ** reference, char * itemName, XdmfInt32 numberType, XdmfInt32 itemType,
                                    XdmfInt32 format, char * function, char * shape, XdmfPointer data)
{
  /* Sanity check */
  /* Check if there is at least one Grid in the collection */
  if (myCollections.empty()) {
    XdmfErrorMessage("Attribute not added, please add one to the collection before continue");
    return(XDMF_FAIL);
  }

  /* Check if last element added to collection stack is an Attribute or another DataItem */
  if ((XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Attribute") == 0) &&
      (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "DataItem") == 0)) {
    XdmfErrorMessage("Cannot add a DataItem inside a non Attribute or DataItem element (" <<
                     myCollections.top()->GetElementType() << ")");
    return(XDMF_FAIL);
  }

  XdmfDataItem * currDataItem = new XdmfDataItem();
  if (itemName) currDataItem->SetName(itemName);
  currDataItem->SetItemType(itemType);
  /* Only set format (XML, HDF) if DataItem is Uniform */
  if (!currDataItem->GetIsMultiple())
    currDataItem->SetFormat(format);
  if (itemType == XDMF_ITEM_FUNCTION) {
    if (function) currDataItem->SetFunction(function);
    else {
      XdmfErrorMessage("Please set a function operation for a XDMF_ITEM_FUNCTION DataItem");
      return(XDMF_FAIL);
    }
  }
  currDataItem->SetDeleteOnGridDelete(true);
  currDataItem->GetDataDesc()->SetShapeFromString(shape);

  /* DataItem must be inserted before checking reference because XdmfXmlNode is
   * created once is inserted.  For DataItem we can not use "Insert" method until
   * we know if is going to be a real reference or an externalHDF reference.
   * Otherwise the DataItem would be inserted twice under its parent.  For this
   * purpose we use the trick of creating an isolated XdmfXmlNode, use it for
   * checking reference and destroy it afterwards. For not externalHDF cases we
   * finally insert a new DataItem element using Insert method. */
  XdmfXmlNode xmlNode = myCollections.top()->GetDOM()->InsertNew(myCollections.top()->GetElement(),
                                                                 currDataItem->GetElementName());

  /* Use reference XdmfXmlNode object as reference for new object */
  if ((reference) && (*reference != NULL)) {
    XdmfXmlNode xmlNode = (XdmfXmlNode) *reference;

    if (currDataItem->CheckForReference(currDataItem->GetElement(), xmlNode) != xmlNode) {
      XdmfErrorMessage("Error checking object Reference for DataItem " << currDataItem << ". Falling in default mode");

      /* Fallback to non-reference mode */
      *reference = NULL;
    }
    else XdmfDebug("Setting DataItem " << currDataItem << " Reference as " << myDOM->GetPath(xmlNode) << " of XdmfXmlNode " << xmlNode);
  }

  /* Delete temporal XdmfXmlNode element added for reference checking */
  myCollections.top()->GetDOM()->DeleteNode( xmlNode );


  /* Otherwise build a new object if is a Uniform DataItem (not Collection) */
  if (((!reference) || (reference && !*reference)) && (!currDataItem->GetIsMultiple())) {
    if (externalHDF) {
      // Inject data to parent (XdmfAttribute or XdmfDataItem)
      // Insert DataItem element from an existing HDF5 dataset
      XdmfDebug("Injecting XML dataitem from " << (char*)data);
      XdmfInjectXml( myCollections.top(), data );
    }
    else {
      XdmfArray * array = currDataItem->GetArray();
      array->SetNumberType(numberType);
      if (array->SetShapeFromString(shape) == XDMF_FAIL) { // dims via shape string
        XdmfErrorMessage("Can not set array shape to: " << shape);
        return(XDMF_FAIL);
      }
      WriteToXdmfArray(array, data);
    }
  }

  /* Only if not running in externalHDF mode or not a leaf DataItem, add XdmfDataItem to the Xdmf class structure */
  if ((!externalHDF) || (currDataItem->GetIsMultiple())) {

    /* Add DataItem to the inner most Attribute/DataItem parent */
    /* Add DataItem into parent Xdmf Element and DOM tree */
    if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Attribute") == 0)
      ((XdmfAttribute*)myCollections.top())->Insert(currDataItem);    /* Add DataItem into parent Xdmf Element and DOM tree */
    else ((XdmfDataItem*)myCollections.top())->Insert(currDataItem);

    /* Only build object if it is not going to be a Uniform Grid (leaf of the tree) */
    if (lastGrid->GetGridType() != XDMF_GRID_UNIFORM) {
      currDataItem->Build();                    /* Add XML/DOM attributes to the XdmfXmlNode object */
    }
    /* Add DataItem to the stack of DataItems */
    else myDataItems.push(currDataItem);

    /* Add Attribute to the collection stack if Collection element */
    if (currDataItem->GetIsMultiple()) {
      myCollections.push(currDataItem);

      /* Also delete Array object in DataItem because this object does not
       * contain any array actually. Solves Dimension attribute */
      currDataItem->SetArray(NULL);
    }

    /* Return XdmfXmlNode reference to this object */
    if (reference) {
      *reference = (void*) currDataItem->GetElement();
      XdmfDebug("Returning DataItem reference " << *reference);
    }
  }
  /* Otherwise delete XdmfDataItem element because it has been injected through raw Xml */
  else {
    delete currDataItem;
  }

  return(XDMF_SUCCESS);
}


/**
 *
 * Close the current open DataItem collection. If within a nested collection, close
 * the most deeply nested collection.
 *
 * Although this call is not mandatory because resources are freed when
 * ~XdmfInterface is called, it is recommendable to call
 * CloseDataItem before closing Grid and WriteGrid.
 *
 */
int
XdmfInterface::CloseDataItem()
{
  /* Sanity check */
  /* Check if there is at least one Grid/Attribute/DataItem in the collection */
  if (myCollections.empty()) {
    XdmfErrorMessage("Neither Grid, Attribute nor DataItem added, please add them to the collection before continue");
    return(XDMF_FAIL);
  }

  /* Check if last element added to collection stack is a DataItem */
  if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "DataItem") == 0) {
    XdmfErrorMessage("Top element in Collection stack is not a DataItem ("
                     << myCollections.top()->GetElementType() << ") : " <<
                     "check your call tree of Add/Close* methods to XdmfInterface");
    return(XDMF_FAIL);
  }

  /* Remove current top DataItem */
  /* DataItems are removed by their Grid parents in Destroy function */
  myCollections.pop();

  return(XDMF_SUCCESS);
}


/**
 *
 * Add an Information element to be written to the last Grid/Attribute/DataItem element
 * of the current collection.
 * If we are not within a collection add to the top level domain
 *
 */
XdmfInt32
XdmfInterface::AddInformation(void ** reference, char * informationName, char * value)
{
  /* Create XdmfInformation object */
  XdmfInformation * currInfo = new XdmfInformation();
  currInfo->SetName(informationName);
  currInfo->SetDeleteOnGridDelete(true);

  /* Information must be inserted before checking reference
   * because XdmfXmlNode is created once is inserted */

  /* Add Information to the outter most Grid */
  if (!myCollections.empty()) myCollections.top()->Insert(currInfo);    /* Add Information into Xdmf Grid and DOM tree */
  else                        myDomain->Insert(currInfo);

  /* Use reference XdmfXmlNode object as reference for new object */
  if ((reference) && (*reference != NULL)) {
    XdmfXmlNode xmlNode = (XdmfXmlNode) *reference;

    if (currInfo->CheckForReference(currInfo->GetElement(), xmlNode) != xmlNode) {
      XdmfErrorMessage("Error checking object Reference for Information " << currInfo << ". Falling in default mode");

      /* Fallback to non-reference mode */
      *reference = NULL;
    }
    else XdmfDebug("Setting Information " << currInfo << " Reference as " << myDOM->GetPath(xmlNode) << " of XdmfXmlNode " << xmlNode);
  }

  /* Otherwise build a new object */
  if ((!reference) || (reference && !*reference)) currInfo->SetValue(value);

  /* Only build object if it is not going to be a Uniform Grid (leaf of the tree) */
  if (!myCollections.empty()) {
    if (lastGrid->GetGridType() != XDMF_GRID_UNIFORM) {
      currInfo->Build();                    /* Add XML/DOM attributes to the XdmfXmlNode object */
    }
    else myInformations.push_back(currInfo);
  }
  else currInfo->Build(); //myDomain->Build();

  /* Return XdmfXmlNode reference to this object */
  if (reference) {
    *reference = (void*) currInfo->GetElement();
    XdmfDebug("Returning Information reference " << *reference);
  }

  return(XDMF_SUCCESS);
}


/**
 *
 * Write out "generic" data to XDMF.  This writes out data to the end of the
 * the current collection.  It is independent of any grids.
 *
 * Where Ids are Assigned
 *   XDMF_SET_TYPE_UNSET -1
 *   XDMF_SET_TYPE_NODE   1
 *   XDMF_SET_TYPE_CELL   2
 *      - Only one DataItem is written
 *
 *   XDMF_SET_TYPE_FACE   3
 *      - Two DataItems are written:
 *          1st: Define CellIds
 *          2nd: Define FaceIds in given CellIds
 *
 *   XDMF_SET_TYPE_EDGE   4
 *      - Three DataItems are written:
 *          1st: Define CellIds
 *          2nd: Define FaceIds in given CellIds
 *          3rd: Define Edge in given FaceIds
 *
 *  When externalHDF is enabled several HDF5 URIs are separated through semicolons:
 *    i.e.: test.h5:/CellIds;test.h5:/FaceIds;test.h5:/EdgeIds;test.h5:/Data
 *
 */
XdmfInt32
XdmfInterface::AddGridSet(void ** reference, char * setName, XdmfInt32 numberType, XdmfInt32 setType,
                          XdmfInt32 attributeType, XdmfInt32 rank, XdmfInt64 *dimensions,
                          XdmfPointer cellIds, XdmfPointer faceIds, XdmfPointer ids, XdmfPointer data)
{
  std::stringstream shape;
  for (int i= 0; i< rank; i++)
    shape << dimensions[i] << " ";
  XdmfDebug("Shape is: " << shape.str());
  return this->AddGridSetFromShape(reference, setName, numberType, setType, attributeType,
                                   (char*) shape.str().c_str(), cellIds, faceIds, ids, data);
}


/**
 *
 * Write out "generic" data to XDMF.  This writes out data to the end of the
 * the current collection.  It is independent of any grids.
 *
 */
XdmfInt32
XdmfInterface::AddGridSetFromShape(void ** reference, char * setName, XdmfInt32 numberType, XdmfInt32 setType,
                                   XdmfInt32 attributeType, char * shape, XdmfPointer cellIds, XdmfPointer faceIds,
                                   XdmfPointer ids, XdmfPointer data)
{
  /* Sanity check */
  /* Check if there is at least one Grid in the collection */
  if (myCollections.empty()) {
    XdmfErrorMessage("Grid not added, please add one to the collection before continue adding Set");
    return(XDMF_FAIL);
  }

  /* Check if last element added to collection stack is a Grid */
  if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Grid") == 0) {
    XdmfErrorMessage("Cannot add a Set inside a non Grid element (" <<
                     myCollections.top()->GetElementType() << ")");
    return(XDMF_FAIL);
  }

  XdmfSet * currSet = new XdmfSet();
  currSet->SetName(setName);
  currSet->SetSetType(setType);

  /* Check whether SetType and number of Ids matches */
  if (!ids) {
    XdmfErrorMessage("Ids vector is mandatory for " << currSet->GetSetTypeAsString() << " set.");
    return(XDMF_FAIL);
  }
  switch( currSet->GetSetType() ) {
    case XDMF_SET_TYPE_NODE:
    case XDMF_SET_TYPE_CELL:
      break;
    case XDMF_SET_TYPE_EDGE:
      if (!faceIds) {
        XdmfErrorMessage("FaceIds vector is mandatory for " << currSet->GetSetTypeAsString() << " Sets.");
        return(XDMF_FAIL);
      }
    case XDMF_SET_TYPE_FACE:
      if (!cellIds) {
        XdmfErrorMessage("CellIds vector is mandatory for " << currSet->GetSetTypeAsString() << " Sets.");
        return(XDMF_FAIL);
      }
      break;
    default:
      XdmfErrorMessage("Unknown Set type, please check Xdmf manual.");
      return(XDMF_FAIL);
  }

  currSet->GetShapeDesc()->SetShapeFromString(shape);
  currSet->SetDeleteOnGridDelete(true);

  /* Set must be inserted before checking reference
   * because XdmfXmlNode is created once is inserted */

  /* Add Set to the outter most Grid */
  myCollections.top()->Insert(currSet);      /* Add Set into Xdmf Grid and DOM tree */

  /* Use reference XdmfXmlNode object as reference for new object */
  if ((reference) && (*reference != NULL)) {
    XdmfXmlNode xmlNode = (XdmfXmlNode) *reference;

    if (currSet->CheckForReference(currSet->GetElement(), xmlNode) != xmlNode) {
      XdmfErrorMessage("Error checking object Reference for Set " << currSet << ". Falling in default mode");

      /* Fallback to non-reference mode */
      *reference = NULL;
    }
    else XdmfDebug("Setting Set " << currSet << " Reference as " << myDOM->GetPath(xmlNode) << " of XdmfXmlNode " << xmlNode);
  }

  /* Otherwise build a new object */
  if ((!reference) || (reference && !*reference)) {
    if (externalHDF) {
      // Insert DataItem/s and Attribute element from existing HDF5 datasets
      std::stringstream dataXml;

      switch( currSet->GetSetType() ) {
        case XDMF_SET_TYPE_NODE:
        case XDMF_SET_TYPE_CELL:
          break;
        case XDMF_SET_TYPE_FACE:
        case XDMF_SET_TYPE_EDGE:
          XdmfDebug("Injecting XML dataitem from " << (char*)cellIds);
          XdmfGetXml( cellIds, dataXml );

          /* Only Edge type must write Faces ids */
          if (currSet->GetSetType() == XDMF_SET_TYPE_FACE) break;

          XdmfDebug("Injecting XML dataitem from " << (char*)faceIds);
          XdmfGetXml( faceIds, dataXml );

          break;
        default:
          break;
      }

      XdmfDebug("Injecting XML dataitem from " << (char*)ids);
      XdmfGetXml( ids, dataXml );

      /* If there are also Attributes for current Set add them */
      if (data) {
        XdmfDebug("Injecting XML dataitem from " << (char*)data);
        XdmfGetXml( data, dataXml );
      }

      /* Inject raw Xml to the current Set */
      currSet->SetDataXml( (char*)dataXml.str().c_str() );

    }
    else {
      // Copy Elements from Set to XdmfArray
      XdmfArray *idsArray;
      std::stringstream heavyDataName;

      switch( currSet->GetSetType() ) {
        case XDMF_SET_TYPE_NODE:
        case XDMF_SET_TYPE_CELL:
          break;
        case XDMF_SET_TYPE_EDGE:

          idsArray = currSet->GetFaceIds();
          idsArray->SetNumberType(numberType);
          if (idsArray->SetShapeFromString(shape) == XDMF_FAIL) { // dims via shape string
            XdmfErrorMessage("Can not set faceIdsArray shape to: " << shape);
            return(XDMF_FAIL);
          }
          heavyDataName << myName << ".h5:/" << setName << "FaceIds";

          idsArray->SetHeavyDataSetName(heavyDataName.str().c_str());
          WriteToXdmfArray(idsArray, faceIds);

        case XDMF_SET_TYPE_FACE:

          idsArray = currSet->GetCellIds();
          idsArray->SetNumberType(numberType);
          if (idsArray->SetShapeFromString(shape) == XDMF_FAIL) { // dims via shape string
            XdmfErrorMessage("Can not set cellIdsArray shape to: " << shape);
            return(XDMF_FAIL);
          }
          heavyDataName.str("");
          heavyDataName.clear();
          heavyDataName << myName << ".h5:/" << setName << "CellIds";

          idsArray->SetHeavyDataSetName(heavyDataName.str().c_str());
          WriteToXdmfArray(idsArray, faceIds);

          break;
        default:
          break;
      }

      idsArray = currSet->GetIds();
      idsArray->SetNumberType(numberType);
      if (idsArray->SetShapeFromString(shape) == XDMF_FAIL) { // dims via shape string
        XdmfErrorMessage("Can not set idsArray shape to: " << shape);
        return(XDMF_FAIL);
      }
      heavyDataName.str("");
      heavyDataName.clear();
      heavyDataName << myName << ".h5:/" << setName << "Ids";

      idsArray->SetHeavyDataSetName(heavyDataName.str().c_str());
      WriteToXdmfArray(idsArray, ids);

      /* If there are also Attributes for current Set add them */
      if (data) {
        XdmfAttribute *dataAttribute = new XdmfAttribute();
        std::stringstream dataName;
        dataName << setName << "Data";
        dataAttribute->SetName(dataName.str().c_str());

        switch( currSet->GetSetType() ) {
          case XDMF_SET_TYPE_NODE: dataAttribute->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);break;
          case XDMF_SET_TYPE_CELL: dataAttribute->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_CELL);break;
          case XDMF_SET_TYPE_FACE: dataAttribute->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_FACE);break;
          case XDMF_SET_TYPE_EDGE: dataAttribute->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_EDGE);break;
          default: break;
        }

        dataAttribute->SetAttributeType(attributeType);
        dataAttribute->SetDeleteOnGridDelete(true);

        XdmfArray *array = dataAttribute->GetValues();
        array->SetNumberType(numberType);
        if (array->SetShapeFromString(shape) == XDMF_FAIL) { // dims via shape string
          XdmfErrorMessage("Can not set dataAttribute shape to: " << shape);
          return(XDMF_FAIL);
        }

        heavyDataName << myName << ".h5:/" << setName << "Data";
        array->SetHeavyDataSetName(heavyDataName.str().c_str());

        WriteToXdmfArray(array, data);

        /* Add Attribute into Set */
        currSet->Insert(dataAttribute);
      }
    }
  }

  /* Only build object if it is not going to be a Uniform Grid (leaf of the tree) */
  if (((XdmfGrid*)myCollections.top())->GetGridType() != XDMF_GRID_UNIFORM) {
    currSet->Build();                    /* Add XML/DOM attributes to the XdmfXmlNode object */
  }

  /* Return XdmfXmlNode reference to this object */
  if (reference) {
    *reference = (void*) currSet->GetElement();
    XdmfDebug("Returning Set reference " << *reference);
  }

  return(XDMF_SUCCESS);
}

/**
 *
 * Read a Xdmf file into the current XdmfDOM.  Must call XdmfReadGrid() to read in associated geometry
 * topology, and attributes.
 *
 */
void
XdmfInterface::ReadFile(char * filePath)
{
  // Clear state and start over before reading file
  this->Destroy();

  myDOM = new XdmfDOM();
  myRoot = new XdmfRoot();
  myDomain = new XdmfDomain();

  myDOM->Parse(filePath);
  myDomain->SetElement(myDOM->FindElement("Domain"));
  myRoot->SetElement(myDOM->GetRoot());

  // Perhaps we should support collections more on this part?
  while (!myCollections.empty())
  {
    delete myCollections.top();
    myCollections.pop();
  }
  myGridPaths.clear();
  myGridNames.clear();
  this->ReadFilePriv(myDomain->GetElement());
}

void
XdmfInterface::ReadFilePriv(XdmfXmlNode currElement)
{
  XdmfGrid currGrid = XdmfGrid();
  for (int i = 0; i < myDOM->FindNumberOfElements("Grid", currElement); i++)
  {
    currGrid.SetDOM(myDOM);
    currGrid.SetElement(myDOM->FindElement("Grid", i, currElement));
    currGrid.Update();
    if (currGrid.GetGridType() != XDMF_GRID_COLLECTION)
    {
      myGridPaths.push_back(myDOM->GetPath(currGrid.GetElement()));
      std::string gridName = currGrid.GetName();
      if (gridName.find_last_of("_") != std::string::npos)
      {
        try
        {
          atoi(gridName.substr(gridName.find_last_of("_") + 1, gridName.size()
              - gridName.find_last_of("_")).c_str());
          gridName = gridName.substr(0, gridName.find_last_of("_"));
        }
        catch (int)
        {
        }
      }
      std::string stdGridName = gridName;
      if (myGridNames.find(stdGridName) == myGridNames.end())
      {
        myGridNames[stdGridName] = 1;
      }
      else
      {
        myGridNames[stdGridName]++;
      }
    }
    this->ReadFilePriv(currGrid.GetElement());
  }
}

/**
 *
 * Read a grid in the current XdmfDOM into XdmfGeometry, XdmfTopology, and XdmfAttribute elements.
 * An XdmfReadGrid() followed by a XdmfWriteGrid() will make a copy of the grid.
 *
 */
void
XdmfInterface::ReadGrid(char * gridName)
{
  ReadGridPriv(gridName, myDomain->GetElement());
}

/**
 *
 * Read a grid in the current XdmfDOM into XdmfGeometry, XdmfTopology, and XdmfAttribute elements.
 * An XdmfReadGrid() followed by a XdmfWriteGrid() will make a copy of the grid.
 *
 */
void
XdmfInterface::ReadGridAtIndex(int gridIndex)
{
  ReadGridPriv(myGridPaths[gridIndex].c_str());
}

/**
 *
 * Helper function for XdmfReadGrid.  Ensures that all grids are traversed and that the method works
 * even within collections.
 *
 */
void
XdmfInterface::ReadGridPriv(char * gridName, XdmfXmlNode currElement)
{
  XdmfGrid currGrid = XdmfGrid();
  for (int i = 0; i < myDOM->FindNumberOfElements("Grid", currElement); i++)
  {
    currGrid.SetDOM(myDOM);
    currGrid.SetElement(myDOM->FindElement("Grid", i, currElement));
    currGrid.Update();
    if (currGrid.GetGridType() != XDMF_GRID_COLLECTION)
    {
      if (strcmp(gridName, currGrid.GetName()) == 0)
      {
        return this->ReadGridPriv(myDOM->GetPath(currGrid.GetElement()));
      }
    }
    this->ReadGridPriv(gridName, currGrid.GetElement());
  }
}

/**
 *
 * Helper function for XdmfReadGrid.  Ensures that all grids are traversed and that the method works
 * even within collections.
 *
 */
void
XdmfInterface::ReadGridPriv(XdmfConstString gridPath)
{
  XdmfGrid currGrid = XdmfGrid();
  currGrid.SetDOM(myDOM);
  currGrid.SetElement(myDOM->FindElementByPath(gridPath));
  currGrid.Update();

  delete myGeometry;
  delete myTopology;

  myGeometry = new XdmfGeometry();
  myGeometry->SetGeometryType(currGrid.GetGeometry()->GetGeometryType());
  myGeometry->SetNumberOfPoints(currGrid.GetGeometry()->GetNumberOfPoints());
  myGeometry->SetPoints(currGrid.GetGeometry()->GetPoints()->Clone());

  myTopology = new XdmfTopology();
  myTopology->SetTopologyType(currGrid.GetTopology()->GetTopologyType());
  myTopology->SetNumberOfElements(currGrid.GetTopology()->GetNumberOfElements());
  myTopology->SetConnectivity(currGrid.GetTopology()->GetConnectivity()->Clone());

  while (!myAttributes.empty())
  {
    delete myAttributes.back();
    myAttributes.pop_back();
  }

  for (int j = 0; j < currGrid.GetNumberOfAttributes(); j++)
  {
    currGrid.GetAttribute(j)->Update();
    XdmfAttribute * currAttribute = new XdmfAttribute();
    currAttribute->SetName(currGrid.GetAttribute(j)->GetName());
    currAttribute->SetAttributeCenter(currGrid.GetAttribute(j)->GetAttributeCenter());
    currAttribute->SetAttributeType(currGrid.GetAttribute(j)->GetAttributeType());
    currAttribute->SetDeleteOnGridDelete(true);
    if (!externalHDF) currAttribute->SetValues(currGrid.GetAttribute(j)->GetValues()->Clone());
    myAttributes.push_back(currAttribute);
  }

  for (int j = 0; j < currGrid.GetNumberOfInformations(); j++)
  {
    currGrid.GetInformation(j)->UpdateInformation();
    XdmfInformation * currInformation = new XdmfInformation();
    currInformation->SetName(currGrid.GetInformation(j)->GetName());
    currInformation->SetValue(currGrid.GetInformation(j)->GetValue());
    currInformation->SetDeleteOnGridDelete(true);
    myInformations.push_back(currInformation);
  }
}

/**
 *
 * Returns the number of grids in the current open file.  This ignores collections.
 *
 */
XdmfInt32
XdmfInterface::GetNumberOfGrids()
{
  return myGridPaths.size();
}

/**
 *
 * Returns the number of elements in the current open grid (the current active XdmfTopology Element). This is either
 * from a current read-in file or from a created but unwritten grid.  If no topology element is present, return -1.
 *
 */
XdmfInt32
XdmfInterface::GetNumberOfElements(XdmfInt64 * dimensions)
{
  if (myTopology != NULL)
    return (myTopology->GetShapeDesc()->GetShape( dimensions ));
  else return -1;
}

/**
 *
 * Returns the number of points in the current open grid (the current active XdmfGeometry Element).  This is either
 * from a current read-in file or from a created but unwritten grid.  If no geometry element is present, return -1.
 *
 */
XdmfInt64
XdmfInterface::GetNumberOfPoints()
{
  if (myGeometry != NULL)
    return (myGeometry->GetNumberOfPoints());
  else return -1;
}

/**
 *
 * Reads the point values from the current geometry into the passed array
 * pointer.  If the geometry has not been created no values are read.
 *
 */
void
XdmfInterface::ReadPointValues(XdmfInt32 numberType, XdmfInt32 startIndex, XdmfPointer arrayToFill,
                               XdmfInt32 numberOfValues, XdmfInt32 arrayStride, XdmfInt32 valuesStride)
{
  if (myGeometry != NULL)
  {
    this->ReadFromXdmfArray(myGeometry->GetPoints(), numberType, startIndex, arrayToFill,
        numberOfValues, arrayStride, valuesStride);
  }
}

/**
 *
 * Returns the number of values in the specified attribute.  Iterates over all
 * current open attributes to find the specified attribute name and returns the
 * number of values it contains.  If no attribute is found, return -1.
 *
 */
XdmfInt32
XdmfInterface::GetNumberOfAttributeValues(char * attributeName)
{
  /* Figure out how to obtain Shape using externalHDF */
  if (externalHDF) return -1;

  for (unsigned int i = 0; i < myAttributes.size(); i++)
    if (strcmp(myAttributes[i]->GetName(), attributeName) == 0)
      return ((XdmfInt32) myAttributes[i]->GetValues()->GetNumberOfElements());

  return -1;
}

/**
 *
 * Reads the values from the specified attribute into the passed array pointer.
 * If the attribute cannot be found, no values are read.
 *
 */
void
XdmfInterface::ReadAttributeValues(char * attributeName, XdmfInt32 numberType, XdmfInt32 startIndex,
    XdmfPointer arrayToFill, XdmfInt32 numberOfValues, XdmfInt32 arrayStride, XdmfInt32 valuesStride)
{
  if (externalHDF) return;

  for (unsigned int i = 0; i < myAttributes.size(); i++)
  {
    if (strcmp(myAttributes[i]->GetName(), attributeName) == 0)
    {
      this->ReadFromXdmfArray(myAttributes[i]->GetValues(), numberType, startIndex, arrayToFill,
          numberOfValues, arrayStride, valuesStride);
    }
  }
}

/**
 *
 * Return the values from the specified information element.  If the
 * information element cannot be found, no values are passed.  Information
 * elements at the top level domain are searched first, followed by the
 * currently loaded grid.
 *
 */
const char *
XdmfInterface::ReadInformationValue(char * informationName)
{
  // TODO: Make this work better for collections as well!
  for (unsigned int i = 0; i < myInformations.size(); i++)
  {
    if (strcmp(informationName, myInformations[i]->GetName()) == 0)
      return myInformations[i]->GetValue();
  }

  for (int i = 0; i < myDOM->FindNumberOfElements("Information", myDomain->GetElement()); i++)
  {
    XdmfInformation currInfo = XdmfInformation();
    currInfo.SetDOM(myDOM);
    currInfo.SetElement(myDOM->FindElement("Information", i, myDomain->GetElement(), 0));
    currInfo.UpdateInformation();
    if (strcmp(informationName, currInfo.GetName()) == 0)
      return currInfo.GetValue();
  }

  return NULL;
}

/**
 *
 * Add a grid to the XdmfDOM.  Assign the current topology, geometry, and grid
 * attributes to grid.  If within a collection, add grid to the collection,
 * otherwise add to the top level domain.  Assign time value if value is
 * nonnegative.
 *
 */
XdmfInt32
XdmfInterface::WriteGrid()
{
  XdmfGrid * grid;

  /* Check if there is at least one Grid in the collection */
  if (myCollections.empty()) {
    XdmfErrorMessage("Uniform Grid not added, please add one to the collection before continue");
    return(XDMF_FAIL);
  }

  /* Check if last element added to collection stack is a Grid */
  /* Close previous elements till reach the first Grid in stack */
  if (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Grid") == 0) {
    XdmfErrorMessage("Top element in Collection stack is not a Grid ("
                     << myCollections.top()->GetElementType() << ") : " <<
                     "check your call tree of Add/Close* methods to XdmfInterface");
    XdmfErrorMessage("Recovering status by closing any other opened elements");

    /* Don't delete elements, only pop from stack (they are deleted by parent Grid) */
    while ( !myCollections.empty() && (XDMF_WORD_CMP(myCollections.top()->GetElementType(), "Grid") == 0 )) {
      myCollections.pop();
    }

    /* Check if there is at least one Grid in the collection */
    if (myCollections.empty()) {
      XdmfErrorMessage("Uniform Grid not added, please add one to the collection before continue");
      return(XDMF_FAIL);
    }
  }

  /* Warn if top Grid in collection is not Uniform type */
  grid = (XdmfGrid*)myCollections.top();
  if (grid->GetGridType() != XDMF_GRID_UNIFORM) {
    XdmfErrorMessage("Cannot write Grid that is not Uniform type");
    return(XDMF_FAIL);
  }

  if (myTopology == NULL)
  {
    XdmfErrorMessage("Must set a topology before the grid can be written");
    return(XDMF_FAIL);
  }

  if (myGeometry == NULL)
  {
    XdmfErrorMessage("Must set a geometry before the grid can be written");
    return(XDMF_FAIL);
  }


  // Set Topology Name
  // Only change HDF5 name if not externalHDF is supplied (otherwise creates connectivity array)
  if ((myTopology->GetClass() != XDMF_STRUCTURED) &&
      (myTopology->GetTopologyType() != XDMF_POLYVERTEX) && (!externalHDF))
  {
    // Modify HDF5 names so we aren't writing over top of our data!
    std::stringstream topologyDataName;
    topologyDataName << myName << ".h5:/" << grid->GetName() << "/CONN";
    myTopology->GetConnectivity()->SetHeavyDataSetName(topologyDataName.str().c_str());
  }

  // Set Geometry Name
  // Only change HDF5 name if not externalHDF is supplied (otherwise creates points array)
  if (!externalHDF) {
    std::stringstream geometryDataName;
    geometryDataName << myName << ".h5:/" << grid->GetName() << "/XYZ";
    myGeometry->GetPoints()->SetHeavyDataSetName(geometryDataName.str().c_str());
  }


  /* Add Time information */
  /* TODO: Future improvements to support more complex Time objects */
  XdmfTime * t = grid->GetTime(); //new XdmfTime();
  if (currentTime >= 0)
  {
    t->SetTimeType(XDMF_TIME_SINGLE);
    t->SetValue(currentTime);
    grid->SetBuildTime(true);
    currentTime = -1;
  }

  /* XdmfInformation object is deleted in Grid destructor */
  while (myInformations.size() > 0)
  {
    myInformations.pop_back();
  }

  /* XdmfDataItem object is deleted in Grid destructor */
  while (myDataItems.size() > 0)
  {
    myDataItems.pop();
  }

  /* XdmfAttribute object is deleted in Grid destructor */
  while (myAttributes.size() > 0)
  {
    XdmfAttribute * currAttribute = myAttributes.back();
    if (!externalHDF) {
      std::stringstream attributeDataName;
      attributeDataName << myName << ".h5:/" << grid->GetName() << "/"
          << currAttribute->GetName();
      currAttribute->GetValues()->SetHeavyDataSetName(attributeDataName.str().c_str());
    }

    myAttributes.pop_back();
  }

  /* Build DOM structure / LightData / HeavyData for Grid */
  grid->Build();

  /* Add XPath to the GridPath vector */
  myGridPaths.push_back(myDOM->GetPath(grid->GetElement()));

  /* XdmfTopology and XdmfGeometry will be deleted in Grid destructor */
  myTopology = NULL;
  myGeometry = NULL;

  return(XDMF_SUCCESS);
}

/**
 *
 * Write constructed Xdmf file to disk with filename created upon initialization
 *
 */
void
XdmfInterface::WriteToFile()
{
  std::stringstream dataName;
  dataName << myName << ".xmf";
  myDOM->Write(dataName.str().c_str());
}

/**
 *
 * Print current XdmfDOM to console
 *
 */
void
XdmfInterface::Serialize()
{
  cout << myDOM->Serialize() << endl;
}

/**
 *
 * Copy current XdmfDOM to memory pointed by charPointer
 *
 */
void
XdmfInterface::GetDOM(char * charPointer)
{
  strcpy(charPointer, myDOM->Serialize());
}

/**
 *
 * Helper function to write different datatypes to an XdmfArray.
 *
 */
void
XdmfInterface::WriteToXdmfArray(XdmfArray * array, XdmfPointer data)
{
  /* Do nothing if HDF is read/written externally.
   * XdmfArray.SetValues copies data into a temporary
   * buffer which it is not needed in this case.
   * Safety return
   */
  if (externalHDF) return;

  switch (array->GetNumberType())
    {
  case XDMF_INT8_TYPE:
    array->SetValues(0, (XdmfInt8*) data, array->GetNumberOfElements());
    return;
  case XDMF_INT16_TYPE:
    array->SetValues(0, (XdmfInt16*) data, array->GetNumberOfElements());
    return;
  case XDMF_INT32_TYPE:
    array->SetValues(0, (XdmfInt32*) data, array->GetNumberOfElements());
    return;
  case XDMF_INT64_TYPE:
    array->SetValues(0, (XdmfInt64*) data, array->GetNumberOfElements());
    return;
  case XDMF_FLOAT32_TYPE:
    array->SetValues(0, (XdmfFloat32*) data, array->GetNumberOfElements());
    return;
  case XDMF_FLOAT64_TYPE:
    array->SetValues(0, (XdmfFloat64*) data, array->GetNumberOfElements());
    return;
  case XDMF_UINT8_TYPE:
    array->SetValues(0, (XdmfUInt8*) data, array->GetNumberOfElements());
    return;
  case XDMF_UINT16_TYPE:
    array->SetValues(0, (XdmfUInt16*) data, array->GetNumberOfElements());
    return;
  case XDMF_UINT32_TYPE:
    array->SetValues(0, (XdmfUInt32*) data, array->GetNumberOfElements());
    return;
  default:
    array->SetValues(0, (XdmfFloat64*) data, array->GetNumberOfElements());
    return;
    }
}

/**
 *
 * Helper function to read different datatypes from an XdmfArray.
 *
 */
void
XdmfInterface::ReadFromXdmfArray(XdmfArray * array, char * numberType, XdmfInt32 * startIndex,
    XdmfPointer * arrayToFill, XdmfInt32 * numberOfValues, XdmfInt32 * arrayStride,
    XdmfInt32 * valuesStride)
{
  /* Do nothing if HDF is read/written externally.
   * XdmfArray.SetValues copies data into a temporary
   * buffer which it is not needed in this case.
   * Safety return
   */
  if (externalHDF) return;

  XdmfArray a = XdmfArray();
  a.SetNumberTypeFromString(numberType);
  switch (a.GetNumberType())
    {
  case XDMF_INT8_TYPE:
    array->GetValues(*startIndex, (XdmfInt8*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
  case XDMF_INT16_TYPE:
    array->GetValues(*startIndex, (XdmfInt16*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
  case XDMF_INT32_TYPE:
    array->GetValues(*startIndex, (XdmfInt32*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
  case XDMF_INT64_TYPE:
    array->GetValues(*startIndex, (XdmfInt64*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
  case XDMF_FLOAT32_TYPE:
    array->GetValues(*startIndex, (XdmfFloat32*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
  case XDMF_FLOAT64_TYPE:
    array->GetValues(*startIndex, (XdmfFloat64*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
  case XDMF_UINT8_TYPE:
    array->GetValues(*startIndex, (XdmfUInt8*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
  case XDMF_UINT16_TYPE:
    array->GetValues(*startIndex, (XdmfUInt16*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
  case XDMF_UINT32_TYPE:
    array->GetValues(*startIndex, (XdmfUInt32*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
  default:
    array->GetValues(*startIndex, (XdmfFloat64*) arrayToFill, *numberOfValues, *arrayStride,
        *valuesStride);
    return;
    }
}

/**
 *
 * Helper function to read different datatypes from an XdmfArray.
 *
 */
void
XdmfInterface::ReadFromXdmfArray(XdmfArray * array, XdmfInt32 numberType, XdmfInt32 startIndex,
    XdmfPointer arrayToFill, XdmfInt32 numberOfValues, XdmfInt32 arrayStride, XdmfInt32 valuesStride)
{
  /* Do nothing if HDF is read/written externally.
   * XdmfArray.SetValues copies data into a temporary
   * buffer which it is not needed in this case.
   * Safety return
   */
  if (externalHDF) return;

  XdmfArray a = XdmfArray();
  a.SetNumberType(numberType);
  switch (a.GetNumberType())
    {
  case XDMF_INT8_TYPE:
    array->GetValues(startIndex, (XdmfInt8*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
  case XDMF_INT16_TYPE:
    array->GetValues(startIndex, (XdmfInt16*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
  case XDMF_INT32_TYPE:
    array->GetValues(startIndex, (XdmfInt32*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
  case XDMF_INT64_TYPE:
    array->GetValues(startIndex, (XdmfInt64*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
  case XDMF_FLOAT32_TYPE:
    array->GetValues(startIndex, (XdmfFloat32*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
  case XDMF_FLOAT64_TYPE:
    array->GetValues(startIndex, (XdmfFloat64*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
  case XDMF_UINT8_TYPE:
    array->GetValues(startIndex, (XdmfUInt8*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
  case XDMF_UINT16_TYPE:
    array->GetValues(startIndex, (XdmfUInt16*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
  case XDMF_UINT32_TYPE:
    array->GetValues(startIndex, (XdmfUInt32*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
  default:
    array->GetValues(startIndex, (XdmfFloat64*) arrayToFill, numberOfValues, arrayStride,
        valuesStride);
    return;
    }
}


//
// C++ will mangle the name based on the argument list. This tells the
// compiler not to mangle the name so we can call it from 'C' (but
// really Fortran in this case)
//
extern "C"
{

  /**
   *
   * Initialize a new Xdmf file.
   *
   */
  void XDMF_UTILS_DLL
  XdmfInit(long * pointer, char * outputName, int externalHDF)
  {
    XdmfInterface * myPointer = new XdmfInterface(outputName, externalHDF);
    *pointer = (long) myPointer;

    XdmfUtilsDebug("Inside XdmfInit - outputName is: " << outputName);
  }

  void XDMF_UTILS_DLL
  XdmfSetTime(long pointer, XdmfFloat64 t)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->SetTime(t);
  }

  XdmfFloat64 XDMF_UTILS_DLL
  XdmfGetTime(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->GetTime();
  }

  void XDMF_UTILS_DLL
  XdmfSetExternalHDF(long pointer, XdmfInt32 externalHDF)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->SetExternalHDF(externalHDF);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfGetExternalHDF(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->GetExternalHDF();
  }

  void XDMF_UTILS_DLL
  XdmfSetDebug(long pointer, XdmfInt32 value)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->SetDebug(value);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfGetDebug(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->GetDebug();
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfAddGrid(long pointer, char * gridName, XdmfInt32 gridType, XdmfInt32 collectionType)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->AddGrid(gridName, gridType, collectionType);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfCloseGrid(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->CloseGrid();
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfSetGridTopology(long pointer, void ** reference, XdmfInt32 topologyType, XdmfInt32 numberType, XdmfInt32 noeRank,
                      XdmfInt64 * numberOfElements, XdmfInt32 dimRank, XdmfInt64 * dimensions, XdmfPointer conns)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->SetGridTopology(reference, topologyType, numberType, noeRank,
                                      numberOfElements, dimRank, dimensions, conns);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfSetGridTopologyFromShape(long pointer, void ** reference, XdmfInt32 topologyType, XdmfInt32 numberType,
                               char * numberOfElements, char * dimensions, XdmfPointer conns)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->SetGridTopologyFromShape(reference, topologyType, numberType,
                                               numberOfElements, dimensions, conns);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfSetGridGeometry(long pointer, void ** reference, XdmfInt32 geometryType, XdmfInt32 numberType,
                      XdmfInt32 rank, XdmfInt64 * dimensions, XdmfPointer points)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->SetGridGeometry(reference, geometryType, numberType, rank, dimensions, points);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfSetGridGeometryFromShape(long pointer, void ** reference, XdmfInt32 geometryType, XdmfInt32 numberType,
                               char * dimensions, XdmfPointer points)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->SetGridGeometryFromShape(reference, geometryType, numberType, dimensions, points);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfAddGridAttribute(long pointer, void ** reference, char * attributeName, XdmfInt32 numberType,
                       XdmfInt32 attributeCenter, XdmfInt32 attributeType, XdmfInt32 rank,
                       XdmfInt64 *dimensions, XdmfPointer data)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->AddGridAttribute(reference, attributeName, numberType, attributeCenter,
                                       attributeType, rank, dimensions, data);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfAddGridAttributeFromShape(long pointer, void ** reference, char * attributeName, XdmfInt32 numberType,
                                XdmfInt32 attributeCenter, XdmfInt32 attributeType, char * shape,
                                char * units, XdmfPointer data)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->AddGridAttributeFromShape(reference, attributeName, numberType, attributeCenter,
                                                attributeType, shape, units, data);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfCloseAttribute(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->CloseAttribute();
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfAddDataItem(long pointer, void ** reference, char * itemName, XdmfInt32 numberType, XdmfInt32 itemType,
                  XdmfInt32 format, char * function, XdmfInt32 rank, XdmfInt64 *dimensions, XdmfPointer data)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->AddDataItem(reference, itemName, numberType, itemType,
                                  format, function, rank, dimensions, data);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfAddDataItemFromShape(long pointer, void ** reference, char * itemName, XdmfInt32 numberType,
                           XdmfInt32 itemType, XdmfInt32 format, char * function, char * shape, XdmfPointer data)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->AddDataItemFromShape(reference, itemName, numberType, itemType,
                                           format, function, shape, data);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfCloseDataItem(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->CloseDataItem();
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfAddInformation(long pointer, void ** reference, char * informationName, char * value)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->AddInformation(reference, informationName, value);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfAddGridSet(long pointer, void ** reference, char * setName, XdmfInt32 numberType, XdmfInt32 setType,
                 XdmfInt32 attributeType, XdmfInt32 rank, XdmfInt64 *dimensions, XdmfPointer cellIds,
                 XdmfPointer faceIds, XdmfPointer ids, XdmfPointer data)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->AddGridSet(reference, setName, numberType, setType, attributeType,
                                 rank, dimensions, cellIds, faceIds, ids, data);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfAddGridSetFromShape(long pointer, void ** reference, char * setName, XdmfInt32 numberType, XdmfInt32 setType,
                          XdmfInt32 attributeType, char * shape, XdmfPointer cellIds, XdmfPointer faceIds,
                          XdmfPointer ids, XdmfPointer data)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->AddGridSetFromShape(reference, setName, numberType, setType, attributeType,
                                          shape, cellIds, faceIds, ids, data);
  }

  void XDMF_UTILS_DLL
  XdmfReadFile(long pointer, char * filePath)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->ReadFile(filePath);
  }

  void XDMF_UTILS_DLL
  XdmfReadGrid(long pointer, char * gridName)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->ReadGrid(gridName);
  }

  void XDMF_UTILS_DLL
  XdmfReadGridAtIndex(long pointer, int gridIndex)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->ReadGridAtIndex(gridIndex);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfGetNumberOfGrids(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->GetNumberOfGrids();
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfGetNumberOfElements(long pointer, XdmfInt64 * dimensions)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->GetNumberOfElements(dimensions);
  }

  XdmfInt64 XDMF_UTILS_DLL
  XdmfGetNumberOfPoints(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->GetNumberOfPoints();
  }

  void XDMF_UTILS_DLL
  XdmfReadPointValues(long pointer, XdmfInt32 numberType, XdmfInt32 startIndex,
                      XdmfPointer arrayToFill, XdmfInt32 numberOfValues,
                      XdmfInt32 arrayStride, XdmfInt32 valuesStride)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->ReadPointValues(numberType, startIndex, arrayToFill, numberOfValues, arrayStride,
                               valuesStride);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfGetNumberOfAttributeValues(long pointer, char * attributeName)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->GetNumberOfAttributeValues(attributeName);
  }

  void XDMF_UTILS_DLL
  XdmfReadAttributeValues(long pointer, char * attributeName, XdmfInt32 numberType,
                          XdmfInt32 startIndex, XdmfPointer arrayToFill, XdmfInt32 numberOfValues,
                          XdmfInt32 arrayStride, XdmfInt32 valuesStride)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->ReadAttributeValues(attributeName, numberType, startIndex, arrayToFill,
                                   numberOfValues, arrayStride, valuesStride);
  }

  const char * XDMF_UTILS_DLL
  XdmfReadInformationValue(long pointer, char * informationName)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->ReadInformationValue(informationName);
  }

  XdmfInt32 XDMF_UTILS_DLL
  XdmfWriteGrid(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    return myPointer->WriteGrid();
  }

  void XDMF_UTILS_DLL
  XdmfWriteToFile(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->WriteToFile();
  }

  void XDMF_UTILS_DLL
  XdmfSerialize(long pointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->Serialize();
  }

  void XDMF_UTILS_DLL
  XdmfGetDOM(long pointer, char * charPointer)
  {
    XdmfInterface * myPointer = (XdmfInterface *) pointer;
    myPointer->GetDOM(charPointer);
  }

  /**
   *
   * Close XdmfInterface interface and clean memory
   *
   */
  void XDMF_UTILS_DLL
  XdmfClose(long * pointer)
  {
    XdmfInterface ** myPointer = (XdmfInterface **) pointer;
    delete *myPointer;
    *myPointer = NULL;
  }
}

