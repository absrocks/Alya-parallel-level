/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfInterface.h,v 1.2 2009-12-17 18:17:30 kwleiter Exp $  */
/*  Date : $Date: 2009-12-17 18:17:30 $ */
/*  Version : $Revision: 1.2 $ */
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
/*     February 2013                                               */
/*                                                                 */
/*     This software is distributed WITHOUT ANY WARRANTY; without  */
/*     even the implied warranty of MERCHANTABILITY or FITNESS     */
/*     FOR A PARTICULAR PURPOSE.  See the above copyright notice   */
/*     for more information.                                       */
/*                                                                 */
/*******************************************************************/

#include <string>
#include <iostream>

#ifndef XDMFINTERFACE_H_
#define XDMFINTERFACE_H_

using std::cerr;
using std::cout;
using std::cin;
using std::endl;

#if defined(WIN32) && !defined(XDMFSTATIC)

// Windows and DLL configuration
#if defined(XdmfUtils_EXPORTS)
  #define XDMF_UTILS_DLL __declspec(dllexport)
#else
  #define XDMF_UTILS_DLL __declspec(dllimport)
#endif

#else

// Linux or static configuration
#define XDMF_UTILS_DLL

#endif

#ifdef XdmfDebug
  #undef XdmfDebug
#endif
#define XdmfDebug(x) \
  { if ( this->Debug || this->myDOM->GetGlobalDebug() ) { \
      cerr << "XDMF Debug : " << __FILE__ << " line " << __LINE__ << " ("<< x << ")" << "\n"; \
    } \
  }

#define XdmfUtilsDebug(x) \
  { if ( ((XdmfInterface*)*pointer)->GetDebug() ) { \
      cerr << "XDMF Debug : " << __FILE__ << " line " << __LINE__ << " ("<< x << ")" << "\n"; \
    } \
  }


/* Injects raw Xml into an object for the DataItem Elements (can be several).
 * data is a char* to HDF5 dataset (i.e.: test.h5:/CONN;test.h5:/GEOM)
 * Used when externalHDF = TRUE
 */
/*
   XdmfString new_dataset = malloc(strlen("SERIAL:") + strlen(dataset) + 1);
   sprintf( new_dataset, "SERIAL:%s", dataset );
   XdmfString new_dataset = malloc(strlen("SERIAL:")valuesHDF.SetUseSerialFile(true);
   XdmfString new_dataset = malloc(strlen("SERIAL:")std::string newDataset = "SERIAL:";
   XdmfString new_dataset = malloc(strlen("SERIAL:")newDataset.append( (const char*)dataset );
*/
#define XdmfInjectXml( object, dataset ) \
  { \
    XdmfValuesHDF valuesHDF; \
    char *str, *token, *saveptr, *dataXml, *mydataset; \
    if ( (mydataset = (char*) malloc( strlen( (char*)dataset ) + 1 )) == NULL) { \
      XdmfErrorMessage("Error allocating memory when injecting Xml for local dataset " << dataset); \
      return(XDMF_FAIL); \
    } \
    strcpy( mydataset, (char*)dataset ); \
    /* Parse different HDF5 URIs (if any) */ \
    for (str = mydataset, token = str; token != NULL; str = NULL) { \
      if ((token = strtok_r( str, ";", &saveptr )) != NULL) { \
        if ((dataXml = valuesHDF.DataItemFromHDF( (const char*)token, true )) == NULL ) { \
          XdmfErrorMessage("Can not generate DataItem element for " << token << " dataset"); \
          return(XDMF_FAIL); \
        } \
        else { \
          /* Inject raw Xml and force to Build Data XML */ \
          XdmfDebug("Injecting into " << object->GetElementType() << " raw Xml: " << dataXml); \
          object->SetDataXml( dataXml ); \
          object->BuildFromDataXml(0, 1); \
        } \
      } \
    } \
    free(mydataset); \
  }


/* Returns the raw Xml string DataItem Elements in dataxml (can be several).
 * dataset is a char* to HDF5 dataset (i.e.: test.h5:/CONN;test.h5:/GEOM)
 * dataxml must be of type std::stringstream data type to support << append.
 * Used when externalHDF = TRUE
 */
#define XdmfGetXml( dataset, dataxml ) \
  { \
    char *str, *token, *saveptr, *xml; \
    XdmfValuesHDF valuesHDF; \
    for (str = (char*)dataset; token != NULL; str = NULL) { \
      if ((token = strtok_r( str, ";", &saveptr )) != NULL) { \
        if ((xml = valuesHDF.DataItemFromHDF( (const char*)token, true )) == NULL ) { \
          XdmfErrorMessage("Can not generate DataItem element for " << token << " dataset"); \
          return(XDMF_FAIL); \
        } \
        else dataxml << xml;\
      } \
    } \
  }


/* Macro that enables function symbols to be aliased to
 * multiple names (name mangling)
 *
 * Use as follows just before the header of the functions:
 *  F_SYMS(foo__,foo_,FOO,foo,(void *foo_param1, int *foo_param2, double *foo_param3, ... ))
 */
#define F_SYMS(r1,r2,r3,orig,params) \
    void r1 params __attribute__ ((alias (#orig))); \
    void r2 params __attribute__ ((alias (#orig))); \
    void r3 params __attribute__ ((alias (#orig)));


class XDMF_UTILS_DLL XdmfInterface
{
public:
  XdmfInterface(char * outputName, XdmfInt32 externalHDF);
  ~XdmfInterface();
  void        SetTime(XdmfFloat64 t);
  XdmfFloat64 GetTime();
  void        SetExternalHDF(XdmfInt32 externalHDF);
  XdmfInt32   GetExternalHDF();
  void        SetDebug(XdmfInt32 value);
  XdmfInt32   GetDebug();
  XdmfInt32   AddGrid(char * gridName, XdmfInt32 gridType, XdmfInt32 collectionType);
  XdmfInt32   CloseGrid();
  XdmfInt32   SetGridTopology(void ** reference, XdmfInt32 topologyType, XdmfInt32 numberType, XdmfInt32 noeRank,
                              XdmfInt64 * numberOfElements, XdmfInt32 dimRank, XdmfInt64 * dimensions, XdmfPointer conns);
  XdmfInt32   SetGridTopologyFromShape(void ** reference, XdmfInt32 topologyType, XdmfInt32 numberType,
                                       char * numberOfElements, char * dimensions, XdmfPointer conns);
  XdmfInt32   SetGridGeometry(void ** reference, XdmfInt32 geometryType, XdmfInt32 numberType,
                              XdmfInt32 rank, XdmfInt64 * dimensions, XdmfPointer points);
  XdmfInt32   SetGridGeometryFromShape(void ** reference, XdmfInt32 geometryType, XdmfInt32 numberType,
                                       char * dimensions, XdmfPointer points);
  XdmfInt32   AddGridAttribute(void ** reference, char * attributeName, XdmfInt32 numberType, XdmfInt32 attributeCenter,
                               XdmfInt32 attributeType, XdmfInt32 rank, XdmfInt64 *dimensions, XdmfPointer data);
  XdmfInt32   AddGridAttributeFromShape(void ** reference, char * attributeName, XdmfInt32 numberType,
                                        XdmfInt32 attributeCenter, XdmfInt32 attributeType, char * shape,
                                        char * units, XdmfPointer data);
  XdmfInt32   CloseAttribute();
  XdmfInt32   AddDataItem(void ** reference, char * itemName, XdmfInt32 numberType, XdmfInt32 itemType,
                          XdmfInt32 format, char * function, XdmfInt32 rank, XdmfInt64 *dimensions, XdmfPointer data);
  XdmfInt32   AddDataItemFromShape(void ** reference, char * itemName, XdmfInt32 numberType, XdmfInt32 itemType,
                                   XdmfInt32 format, char * function, char * shape, XdmfPointer data);
  XdmfInt32   CloseDataItem();
  XdmfInt32   AddInformation(void ** reference, char * informationName, char * value);
  XdmfInt32   AddGridSet(void ** reference, char * setName, XdmfInt32 numberType, XdmfInt32 setType,
                         XdmfInt32 attributeType, XdmfInt32 rank, XdmfInt64 *dimensions,
                         XdmfPointer cellIds, XdmfPointer faceIds, XdmfPointer ids, XdmfPointer data);
  XdmfInt32   AddGridSetFromShape(void ** reference, char * setName, XdmfInt32 numberType, XdmfInt32 setType,
                                  XdmfInt32 attributeType, char * shape, XdmfPointer cellIds, XdmfPointer faceIds,
                                  XdmfPointer ids, XdmfPointer data);

  void        ReadFile(char * filePath);
  void        ReadGrid(char * gridName);
  void        ReadGridAtIndex(int gridIndex);
  XdmfInt32   GetNumberOfGrids();
  XdmfInt32   GetNumberOfElements(XdmfInt64 * dimensions);
  XdmfInt64   GetNumberOfPoints();
  void        ReadPointValues(XdmfInt32 numberType, XdmfInt32 startIndex, XdmfPointer arrayToFill,
                              XdmfInt32 numberOfValues, XdmfInt32 arrayStride, XdmfInt32 valuesStride);
  XdmfInt32   GetNumberOfAttributeValues(char * attributeName);
  void        ReadAttributeValues(char * attributeName, XdmfInt32 numberType, XdmfInt32 startIndex,
                                  XdmfPointer arrayToFill, XdmfInt32 numberOfValues,
                                  XdmfInt32 arrayStride, XdmfInt32 valuesStride);
  const char* ReadInformationValue(char * informationName);
  XdmfInt32   WriteGrid();
  void        WriteToFile();
  void        Serialize();
  void        GetDOM(char * charPointer);

private:
  void        Destroy();
  void        ReadFilePriv(XdmfXmlNode currElement);
  void        ReadGridPriv(char * gridName, XdmfXmlNode currElement);
  void        ReadGridPriv(XdmfConstString gridPath);
  void        WriteToXdmfArray(XdmfArray * array, XdmfPointer data);
  void        ReadFromXdmfArray(XdmfArray * array, char * numberType, XdmfInt32 * startIndex,
                                XdmfPointer * arrayToFill, XdmfInt32 * numberOfValues, XdmfInt32 * arrayStride,
                                XdmfInt32 * valuesStride);
  void        ReadFromXdmfArray(XdmfArray * array, XdmfInt32 numberType, XdmfInt32 startIndex,
                                XdmfPointer arrayToFill, XdmfInt32 numberOfValues,
                                XdmfInt32 arrayStride, XdmfInt32 valuesStride);

  XdmfInt32                       Debug;
  XdmfDOM                       * myDOM;
  XdmfRoot                      * myRoot;
  XdmfDomain                    * myDomain;
  XdmfTopology                  * myTopology;
  XdmfGeometry                  * myGeometry;
  XdmfGrid                      * lastGrid;
  /* Stack of collection elements added to the Xdmf
   * tree. They are added when used any Add* method and
   * removed from the stack when Close* method is called.
   * Collection elements can be Grids/Attributes/DataItems.
   */
  std::stack<XdmfElement*>        myCollections;
  std::vector<XdmfAttribute*>     myAttributes;
  /* Current Tree of DataItems being built. They are
   * deleted by their owners (parent Grids or DataItems)
   * through SetDeleteOnGridDelete() method.
   * This stack is actually used to know where next
   * DataItem element should hang from
   */
  std::stack<XdmfDataItem*>       myDataItems; // Might be removed?
  std::vector<XdmfInformation*>   myInformations;
  /* Map changed from const char* to std::string type
   * in order to support strings and not char addresses */
  std::map<std::string, int>      myGridNames;
  std::vector<std::string>        myGridPaths;
  std::string                     myName;
  double                          currentTime;
  /* Specifies whether HDF has to be read/written
   * externally or Xdmf library does it
   */
  XdmfInt32                       externalHDF;
};

#endif /* XDMFINTERFACE_H_ */

