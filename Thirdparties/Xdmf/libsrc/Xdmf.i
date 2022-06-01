/*
  To create C++ files:
ICE_INCLUDES="-I${Xdmf_SOURCE_DIR} -I${Xdmf_SOURCE_DIR}/libsrc -I${Xdmf_BINARY_DIR}/libsrc -I${Xdmf_BINARY_DIR}/Ice/libsrc"

XdmfTcl.cxx:
swig -v -c++ -make_default -includeall -tcl -prefix Xdmf  -namespace ${ICE_INCLUDES} -o XdmfTcl.cxx Xdmf.i

XdmfPython.cxx:
swig -v -c++ -make_default -includeall -python -shadow $(ICE_INCLUDES) -o XdmfPython.cxx ${srcdir}/Xdmf.i

XdmfJava.cxx:
swig -v -c++ -make_default -includeall -shadow -java $(ICE_INCLUDES) -o XdmfJava.cxx ${srcdir}/Xdmf.i;  

XdmfC.cxx:
swig -v -c++ -includeall -c $(ICE_INCLUDES) -o XdmfC.cxx Xdmf.i
*/

%rename("GetUniqueNameDef") XdmfObject::GetUniqueName();
%rename("GetUniqueDef") GetUnique();

%rename("BuildFromDataXmlDef") XdmfElement::BuildFromDataXml();
%rename("InsertDef") XdmfElement::Insert();
%rename("SetElementDef") XdmfElement::SetElement(XdmfXmlNode Element);
%rename("GetValuesDef") XdmfElement::GetValues();
%rename("SetNumberTypeDef") XdmfDataDesc::SetNumberType(int);
%rename("SetNumberTypeFromStringDef") XdmfDataDesc::SetNumberTypeFromString(char const *);
%rename("GetCoordinatesDef") XdmfDataDesc::GetCoordinates();
%rename("GetCoordinatesDef2") XdmfDataDesc::GetCoordinates(long long);
%rename("GetCoordinatesAsStringDef") XdmfDataDesc::GetCoordinatesAsString();
%rename("GetCoordinatesAsStringDef2") XdmfDataDesc::GetCoordinatesAsString(long long);
%rename("SetCompressionDef") XdmfDataDesc::SetCompression();
%rename("AddCompoundMemberFromStringDef") XdmfDataDesc::AddCompoundMemberFromString(char const*, char const*, char const*);

/* Missing more constructors */
%rename("GetDataPointerDef") XdmfArray::GetDataPointer();
%rename("ResetDef") XdmfArray::Reset();

%rename("SetValuesDef") XdmfArray::SetValues(long long, char const*);
%rename("SetValuesDef2") XdmfArray::SetValues(long long, char const*, long long);
%rename("SetValuesDef3") XdmfArray::SetValues(long long, char const*, long long, long long);

%rename("SetValuesFromArrayDef") XdmfArray::SetValues(long long, XdmfArray*, long long);
%rename("SetValuesFromArrayDef2") XdmfArray::SetValues(long long, XdmfArray*, long long, long long);
%rename("SetValuesFromArrayDef3") XdmfArray::SetValues(long long, XdmfArray*, long long, long long, long long);

%rename("GetValuesAsInt8Def") XdmfArray::GetValuesAsInt8(long long, char*, long long);
%rename("GetValuesAsInt8Def2") XdmfArray::GetValuesAsInt8(long long, char*, long long, long long);
%rename("SetValuesFromInt8Def") XdmfArray::SetValuesFromInt8(long long, char*, long long);
%rename("SetValuesFromInt8Def2") XdmfArray::SetValuesFromInt8(long long, char*, long long, long long);

%rename("GetValuesAsInt32Def") XdmfArray::GetValuesAsInt32(long long, int*, long long);
%rename("GetValuesAsInt32Def2") XdmfArray::GetValuesAsInt32(long long, int*, long long, long long);
%rename("SetValuesFromInt32Def") XdmfArray::SetValuesFromInt32(long long, int*, long long);
%rename("SetValuesFromInt32Def2") XdmfArray::SetValuesFromInt32(long long, int*, long long, long long);

%rename("GetValuesAsInt64Def") XdmfArray::GetValuesAsInt64(long long, long long*, long long);
%rename("GetValuesAsInt64Def2") XdmfArray::GetValuesAsInt64(long long, long long*, long long, long long);
%rename("SetValuesFromInt64Def") XdmfArray::SetValuesFromInt64(long long, long long*, long long);
%rename("SetValuesFromInt64Def2") XdmfArray::SetValuesFromInt64(long long, long long*, long long, long long);

%rename("GetValuesAsFloat32Def") XdmfArray::GetValuesAsFloat32(long long, float*, long long);
%rename("GetValuesAsFloat32Def2") XdmfArray::GetValuesAsFloat32(long long, float*, long long, long long);
%rename("SetValuesFromFloat32Def") XdmfArray::SetValuesFromFloat32(long long, float*, long long);
%rename("SetValuesFromFloat32Def2") XdmfArray::SetValuesFromFloat32(long long, float*, long long, long long);

%rename("GetValuesAsFloat64Def") XdmfArray::GetValuesAsFloat64(long long, double*, long long);
%rename("GetValuesAsFloat64Def2") XdmfArray::GetValuesAsFloat64(long long, double*, long long, long long);
%rename("SetValuesFromFloat64Def") XdmfArray::SetValuesFromFloat64(long long, double*, long long);
%rename("SetValuesFromFloat64Def2") XdmfArray::SetValuesFromFloat64(long long, double*, long long, long long);

%rename("GetValuesDef") XdmfArray::GetValues();
%rename("GetValuesDef2") XdmfArray::GetValues(long long);
%rename("GetValuesDef3") XdmfArray::GetValues(long long, long long);

%rename("GenerateDef") XdmfArray::Generate(double, double);
%rename("GenerateDef2") XdmfArray::Generate(double, double, long long);

%rename("CloneDef") XdmfArray::Clone(long long);
%rename("CloneDef2") XdmfArray::Clone(long long, long long);

%rename("ReferenceDef") XdmfArray::Reference();
%rename("ReferenceDef2") XdmfArray::Reference(long long);
%rename("GetNextOlderArrayDef") GetNextOlderArray(long long);

%rename("__ParseDef") XdmfDOM::__Parse(char const*);
%rename("ParseDef") XdmfDOM::Parse(char const*);
%rename("GetNumberOfChildrenDef") XdmfDOM::GetNumberOfChildren(XdmfXmlNode);
%rename("IsChildDef") XdmfDOM::IsChild(XdmfXmlNode);
%rename("SerializeDef") XdmfDOM::Serialize();
%rename("WriteDef") XdmfDOM::Write();
%rename("CreateDef") XdmfDOM::Create();
%rename("CreateDef2") XdmfDOM::Create(char const*);
%rename("FindElementDef") XdmfDOM::FindElement(char const*);
%rename("FindElementDef2") XdmfDOM::FindElement(char const*, int);
%rename("FindElementDef3") XdmfDOM::FindElement(char const*, int, XdmfXmlNode);
%rename("FindNextElementDef") XdmfDOM::FindNextElement(char const*, XdmfXmlNode);
%rename("FindDataElementDef") XdmfDOM::FindDataElement();
%rename("FindDataElementDef2") XdmfDOM::FindDataElement(int);
%rename("FindDataElementDef3") XdmfDOM::FindDataElement(int, XdmfXmlNode);
%rename("FindElementByAttributeDef") XdmfDOM::FindElementByAttribute(char const*, char const*);
%rename("FindElementByAttributeDef2") XdmfDOM::FindElementByAttribute(char const*, char const*, int);
%rename("FindNumberOfElementsDef") XdmfDOM::FindNumberOfElements(char const*);
%rename("FindNumberOfElementsByAttributeDef") XdmfDOM::FindNumberOfElementsByAttribute(char const*, char const*);

%rename("GetArrayDef") XdmfDataItem::GetArray();
%rename("GetDataValuesDef") XdmfDataItem::GetDataValues();
%rename("GetDataValuesDef1") XdmfDataItem::GetDataValues(long long);
%rename("GetDataValuesDef2") XdmfDataItem::GetDataValues(long long, long long);
%rename("SetDataValuesDef") XdmfDataItem::SetDataValues(long long, char const*);
%rename("SetDataValuesDef2") XdmfDataItem::SetDataValues(long long, char const*, long long);

%rename("ReadDef") XdmfValues::Read();
%rename("WriteDef") XdmfValues::Write(XdmfArray*);

%rename("OpenDef") XdmfHeavyData::Open();
%rename("OpenDef2") XdmfHeavyData::Open(char const*);
%rename("ReadDef") XdmfHeavyData::Read();

%rename("CdDef") XdmfHDF::Cd();
%rename("CreateDatasetDef") XdmfHDF::CreateDataset();
%rename("CopyArrayDef") CopyArray(XdmfArray*);

%rename("GetConnectivityDef") XdmfTopology::GetConnectivity();
%rename("GetConnectivityDef2") XdmfTopology::GetConnectivity(XdmfArray*);
%rename("GetCellOffsetsDef") XdmfTopology::GetCellOffsets();

%rename("GetPointsDef") XdmfGeometry::GetPoints();
%rename("SetVectorXDef") XdmfGeometry::SetVectorX(XdmfArray*);
%rename("SetVectorYDef") XdmfGeometry::SetVectorY(XdmfArray*);
%rename("SetVectorZDef") XdmfGeometry::SetVectorZ(XdmfArray*);

%rename("FindGridsAtTimeDef") XdmfGrid::FindGridsAtTime(XdmfTime*, XdmfArray*);
%rename("FindGridsAtTimeDef2") XdmfGrid::FindGridsAtTime(XdmfTime*, XdmfArray*,double);

%rename("EvaluateDef") XdmfTime::Evaluate(XdmfGrid*);
%rename("EvaluateDef2") XdmfTime::Evaluate(XdmfGrid*, XdmfArray*);
%rename("EvaluateDef3") XdmfTime::Evaluate(XdmfGrid*, XdmfArray*, int);

%rename("IsValidRange") XdmfTime::IsValid(double, double);

%rename("GetIdsRef") XdmfSet::GetIds();
%rename("GetCellIdsRef") XdmfSet::GetCellIds();
%rename("GetFaceIdsRef") XdmfSet::GetFaceIds();

%rename("GetIdsRef") XdmfMap::GetIds();
%rename("GetMapIndexRef") XdmfMap::GetMapIndex();
%rename("GetMapDataRef") XdmfMap::GetMapData();

%rename("ConfigureUniformDef") XdmfDsm::ConfigureUniform(XdmfDsmComm*, long long);
%rename("ConfigureUniformDef2") XdmfDsm::ConfigureUniform(XdmfDsmComm*, long long, int);
%rename("ReceiveCommandHeaderDef") XdmfDsm::ReceiveCommandHeader(int*, int*, long long*, long long*);
%rename("ReceiveDataDef") XdmfDsm::ReceiveData(int, void*, long long);

%rename("ServiceOnceDef") XdmfDsmBuffer::ServiceOnce();
%rename("ServiceUntilIdleDef") XdmfDsmBuffer::ServiceUntilIdle();
%rename("ServiceLoopDef") XdmfDsmBuffer::ServiceLoop();
%rename("ServiceDef") XdmfDsmBuffer::Service();
%ignore XdmfDsmBufferServiceThread(void*);

%module Xdmf
%{

    /*
#include <XdmfCharArray.h>
    */
#include <XdmfArray.h>
#include <XdmfAttribute.h>
#include <XdmfDOM.h>
#include <XdmfLightData.h>
#include <XdmfInformation.h>
#include <XdmfElement.h>
#include <XdmfDataDesc.h>
#include <XdmfDataStructure.h>
#include <XdmfHDF.h>
#include <XdmfHeavyData.h>
#include <XdmfValues.h>
#include <XdmfValuesXML.h>
#include <XdmfValuesHDF.h>
#include <XdmfExpression.h>
    /*
#include <XdmfFormat.h>
#include <XdmfFormatHDF.h>
#include <XdmfFormatMulti.h>
#include <XdmfFormatXML.h>
    */
#include <XdmfObject.h>
#include <XdmfDomain.h>
#include <XdmfRoot.h>
#include <XdmfTopology.h>
#include <XdmfGeometry.h>
#include <XdmfGrid.h>
#include <XdmfTime.h>
#include <XdmfRegion.h>
#include <XdmfSet.h>
#include <XdmfMap.h>
    /*
#include <XdmfParameter.h>
#include <XdmfTransform.h>
#include <XdmfXNode.h>
#include <XdmfNDGM.h>
    */
#include <XdmfDsm.h>
#include <XdmfDsmMsg.h>
#include <XdmfDsmBuffer.h>
#include <XdmfDsmComm.h>
#ifndef XDMF_NO_MPI
#include <XdmfDsmCommMpi.h>
#endif

#ifndef HAVE_STRTOLL
# define strtoll XDMF_strtoll
inline XDMF_LONG64 XDMF_strtoll(char *str, void*, int)
{
  XDMF_LONG64 result = 0;
  int negative=0;

  while (*str == ' ' || *str == '\t')
    {
    str++;
    }
  if (*str == '+')
    {
    str++;
    }
  else if (*str == '-')
    {
    negative = 1;
    str++;
    }

  while (*str >= '0' && *str <= '9')
    {
    result = (result*10)-(*str++ - '0');
    }
  return negative ? result : -result;
}
#else
# define XDMF_strtoll strtoll
#endif

%}

%include std_string.i

/*
%include XdmfCharArray.h
*/
%include XdmfAttribute.h
%include XdmfArray.h
%include XdmfDOM.h
%include XdmfLightData.h
%include XdmfInformation.h
%include XdmfElement.h
%include XdmfDataDesc.h
%include XdmfDataStructure.h
%include XdmfValues.h
%include XdmfValuesXML.h
%include XdmfValuesHDF.h
%include XdmfExpression.h
/*
%include XdmfFormat.h
%include XdmfFormatHDF.h
%include XdmfFormatMulti.h
%include XdmfFormatXML.h
%include XdmfHDFSupport.h
*/
%include XdmfHeavyData.h
%include XdmfHDF.h
%include XdmfObject.h
%include XdmfDomain.h
%include XdmfRoot.h
%include XdmfTopology.h
%include XdmfGeometry.h
%include XdmfGrid.h
%include XdmfTime.h
%include XdmfRegion.h
%include XdmfSet.h
%include XdmfMap.h
/*
%include XdmfParameter.h
%include XdmfTransform.h
%include XdmfXNode.h
%include XdmfNDGM.h
*/
%include XdmfDsm.h
%include XdmfDsmMsg.h
%include XdmfDsmBuffer.h
%include XdmfDsmComm.h
#ifndef XDMF_NO_MPI
%include XdmfDsmCommMpi.h
#endif

#ifdef SWIGPYTHON
%{
void XdmfSwigException(int code, const char* msg)
{
/*   SWIG_exception(code, msg); */
}
%}
#endif
